"""Trace a photon stream through an optical model (``obdeect ray-tracing trace-stream``).

Ground-frame (G) or telescope-frame (T) photons go in; camera-frame (C) photons
come out. Records pass through unchanged, except that telescopes not selected
are dropped and the efficiency status is updated. One optical model can serve every
telescope, or each telescope id gets its own model; then one process traces a whole array.
The stream component (obdeect-stream) is a required dependency; its frame and chunk
constants are used here. The frames and the efficiency status rules are described in
docs/ray-tracing/photon-stream.md.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import obdeect.stream as obdeect_stream
from obdeect.ray_tracing import __version__
from obdeect.ray_tracing.formats import ModelError
from obdeect.ray_tracing.model import OpticalModel, ThreadPool, _obscurer_support
from obdeect.stream.stage import EFFICIENCY_STATUS_STATES
from obdeect.stream.stream import FRAMES, reject_same_file

FRAME_GROUND = FRAMES["ground"]
FRAME_TELESCOPE = FRAMES["telescope"]
FRAME_CAMERA = FRAMES["camera"]
# Largest photon chunk written (photon stream format limit).
MAX_CHUNK_ROWS = obdeect_stream.max_chunk_rows
SPEED_OF_LIGHT_M_PER_NS = 0.299792458
# Margin above the highest incoming obscurer for the input plane [m] (as simtools).
INPUT_PLANE_MARGIN_M = 1.0
# Largest allowed difference between a cz column and -sqrt(1 - cx^2 - cy^2).
CZ_TOLERANCE = 1e-5
# Efficiency status factors that the optical model applies, in full, in its trace.
OWN_FACTORS = ("mirror_reflectivity", "telescope_transmission")
# Parts of the telescope transmission: all of it, or t(theta) / t0 only.
TRANSMISSION_PARTS = ("full", "angular")
# Efficiency status factor left for a downstream module when only the angular part is applied.
CONSTANT_FACTOR = "telescope_transmission_constant"
# Run-record keys with the telescope pointing [rad]: zenith and azimuth in the array frame.
RUN_POINTING_KEYS = ("telescope_pointing_zenith_rad", "telescope_pointing_azimuth_rad")
# Telescope ids (sim_telarray numbering from 1); the random-number key holds 16 bits.
TELESCOPE_IDS = range(1, 2**16)


def _rotation(alpha: float, beta: float) -> np.ndarray:
    """Matrix M of the row-vector rule v' = v @ M: rotation by alpha about z, then beta about y."""

    def step(angle: float, i: int, j: int) -> np.ndarray:
        matrix = np.eye(3)
        matrix[i, i] = matrix[j, j] = math.cos(angle)
        matrix[i, j], matrix[j, i] = math.sin(angle), -math.sin(angle)
        return matrix

    return step(alpha, 0, 1) @ step(beta, 2, 0)


def ground_to_telescope_rotation(theta_deg: float, phi_deg: float) -> np.ndarray:
    """Rotation R with p_telescope = R p_ground for a telescope at zenith theta, azimuth phi.

    The ground frame has x north, y west and z up (the array frame, as simtools
    ``ground_direction_to_telescope_frame``). phi is the telescope azimuth in that
    frame, not the CORSIKA propagation azimuth of a shower. The telescope frame has
    +z along the optical axis. A shower that propagates along -axis maps to
    direction (0, 0, -1).
    """
    theta, phi = math.radians(theta_deg), math.radians(phi_deg)
    return _rotation(phi, -theta).T


def minimum_input_plane_height_m(
    model: OpticalModel, theta_deg: float, phi_deg: float = 0.0
) -> float:
    """Height [m] of the ground-frame input plane above the telescope sphere centre.

    The plane lies above every incoming obscurer and every mirror in the pointing direction,
    with the axes offset and a margin included (see ``trace_bunches``). A photon starts on
    the plane and moves down to the ground, so a mirror point above the plane would never be
    reached: at grazing pointing the dish reaches far above the obscurers.
    """
    if not 0.0 <= theta_deg < 90.0:
        raise ModelError(f"telescope zenith angle {theta_deg} deg must lie in [0, 90)")
    normal = ground_to_telescope_rotation(theta_deg, phi_deg)[:, 2]
    supports = [
        _obscurer_support(obscurer, normal)
        for obscurer in model.document["optics"]["obscurers"]
        if "incoming" in obscurer.get("legs", [])
    ]
    supports += [support for _, support in model.mirrors_support(normal)]
    support = max(supports, default=0.0)
    return max(0.0, support + float(np.dot(normal, model.axes_offset_m))) + INPUT_PLANE_MARGIN_M


def _run_pointing(metadata: Mapping[str, Any]) -> tuple[float, float] | None:
    """Telescope (zenith, azimuth) [deg] from the run keys (radians), or None if absent."""
    present = [key in metadata for key in RUN_POINTING_KEYS]
    if not any(present):
        return None
    if not all(present):
        raise ModelError(f"run record: give both {' and '.join(RUN_POINTING_KEYS)} or neither")
    zenith, azimuth = (float(metadata[key]) for key in RUN_POINTING_KEYS)
    return (math.degrees(zenith), math.degrees(azimuth))


def _names(text: Any) -> list[str]:
    return [name.strip() for name in str(text).split(",") if name.strip()]


def _provenance(model: OpticalModel) -> str:
    return f"obdeect {__version__} model {model.sha256}"


@dataclass(frozen=True)
class _Transmission:
    """Which part of the telescope transmission the trace applies."""

    part: str
    on_axis: float | None = None

    @property
    def normalise(self) -> bool:
        return self.part == "angular"


def _transmission(model: OpticalModel, part: str) -> _Transmission:
    if part not in TRANSMISSION_PARTS:
        raise ModelError(f"telescope transmission must be one of {TRANSMISSION_PARTS}: {part!r}")
    if part == "full":
        return _Transmission(part)
    block = model.document.get("optics", {}).get("telescope_transmission")
    if block is None:
        raise ModelError(
            "the angular telescope transmission needs a model with telescope_transmission; "
            "this model has none"
        )
    on_axis = float(block["on_axis"])
    if not on_axis > 0.0:
        raise ModelError("the angular telescope transmission needs a positive on-axis value")
    return _Transmission(part, on_axis)


def _own_efficiency_status(
    model: OpticalModel, transmission: _Transmission, metadata: Mapping[str, Any] | None = None
) -> Mapping[str, Any]:
    """The efficiency status record that marks the factors obdeect applies as applied.

    ``metadata`` is an incoming efficiency status: its other factors are kept, and the factors
    obdeect applies are added to it, so that each factor appears once in the record.
    """
    metadata = dict(metadata or {})
    names = _names(metadata.get("factors", ""))
    factors = names + [name for name in _own_factors(model) if name not in names]
    return _efficiency_status_with_applied(
        {**metadata, "factors": ",".join(factors)}, model, transmission
    )


def _own_factors(model: OpticalModel) -> tuple[str, ...]:
    if (
        model.has_structure
        or model.document.get("optics", {}).get("telescope_transmission") is not None
    ):
        return OWN_FACTORS
    return ("mirror_reflectivity",)


def _efficiency_status_with_applied(
    metadata: Mapping[str, Any], model: OpticalModel, transmission: _Transmission
) -> Mapping[str, Any]:
    """The efficiency status record with the factors obdeect applies marked applied."""
    all_names = _names(metadata.get("factors", ""))
    if transmission.normalise and CONSTANT_FACTOR in all_names:
        raise ModelError(f"the efficiency status already names {CONSTANT_FACTOR}")
    names = [name for name in all_names if name in _own_factors(model)]
    if not names:
        return metadata
    updated = dict(metadata)
    for name in names:
        updated[f"{name}.state"] = "applied"
        updated[f"{name}.source"] = _provenance(model)
    if "telescope_transmission" in names:
        updated["telescope_transmission.part"] = transmission.part
        if transmission.normalise:
            updated[f"{CONSTANT_FACTOR}.state"] = "pending"
            updated[f"{CONSTANT_FACTOR}.value"] = transmission.on_axis
            updated["factors"] = ",".join([*_names(updated["factors"]), CONSTANT_FACTOR])
    return updated


@dataclass
class _Factor:
    """One efficiency status entry; for `envelope`, the envelope table of the factor."""

    state: str
    wavelength_nm: Any = None
    value: Any = None


def _efficiency_status_factors(metadata: Mapping[str, Any]) -> dict[str, _Factor]:
    factors = {}
    for name in _names(metadata.get("factors", "")):
        state = str(metadata.get(f"{name}.state", "pending"))
        if state not in EFFICIENCY_STATUS_STATES:
            raise ModelError(f"efficiency_status factor {name!r}: unknown state {state!r}")
        factors[name] = _Factor(
            state, metadata.get(f"{name}.wavelength_nm"), metadata.get(f"{name}.value")
        )
    return factors


def _envelope_table(name: str, factor: _Factor) -> tuple[np.ndarray, np.ndarray]:
    if factor.wavelength_nm is None or factor.value is None:
        raise ModelError(
            f"efficiency status factor {name!r} is 'envelope' but has no "
            f"{name}.wavelength_nm and {name}.value"
        )
    nodes = np.asarray(factor.wavelength_nm, dtype=np.float64).ravel()
    values = np.asarray(factor.value, dtype=np.float64).ravel()
    if (
        nodes.size == 0
        or nodes.shape != values.shape
        or not np.all(np.diff(nodes) > 0)
        or not np.all(values > 0)
    ):
        raise ModelError(
            f"efficiency_status factor {name!r}: the envelope needs increasing wavelengths and "
            "positive values of the same length"
        )
    return nodes, values


@dataclass(frozen=True)
class _Envelopes:
    """Envelope tables of the factors with state `envelope`."""

    tables: tuple[tuple[np.ndarray, np.ndarray], ...] = ()

    def inverse(self, wavelength_nm: np.ndarray) -> np.ndarray:
        """Return 1 / f_env(wavelength): with the traced f, the factor becomes f / f_env."""
        result = np.ones(wavelength_nm.shape)
        for nodes, values in self.tables:
            if wavelength_nm.size and (
                wavelength_nm.min() < nodes[0] or wavelength_nm.max() > nodes[-1]
            ):
                raise ModelError(
                    f"photon wavelengths lie outside the efficiency status envelope "
                    f"[{nodes[0]:g}, {nodes[-1]:g}] nm"
                )
            result /= np.interp(wavelength_nm, nodes, values)
        return result


@dataclass(frozen=True)
class _Optics:
    """Optical model of a telescope, with the part of its transmission that the trace applies."""

    model: OpticalModel
    transmission: _Transmission
    # Input plane of telescope-frame photons (telescope frame, no tilt): above the incoming
    # obscurers and the mirrors, so that every photon starts upstream of them.
    telescope_plane_m: float

    @classmethod
    def of(cls, model: OpticalModel, part: str) -> _Optics:
        tops = [
            z for z in (model.max_incoming_obscurer_z_m(), model.max_mirror_z_m()) if z is not None
        ]
        return cls(model, _transmission(model, part), max(tops, default=0.0) + INPUT_PLANE_MARGIN_M)


@dataclass(frozen=True)
class _Config:
    """How the traced photons of one telescope are weighted."""

    roulette: bool
    envelopes: _Envelopes


@dataclass
class _Telescope:
    """The telescope block that is being read."""

    telescope_id: int
    kept: bool
    factors: dict[str, _Factor] = field(default_factory=dict)
    config: _Config | None = None
    photons_seen: bool = False
    status_written: bool = False
    out_index: int = 0
    weight: float = 0.0


def _telescope_id(value: int) -> int:
    if value not in TELESCOPE_IDS:
        raise ModelError(
            f"telescope id {value} outside {TELESCOPE_IDS.start}..{TELESCOPE_IDS.stop - 1}"
        )
    return value


def _wavelength_law(metadata: Mapping[str, Any]) -> tuple[float, float] | None:
    if "wavelength_min_nm" in metadata and "wavelength_max_nm" in metadata:
        lo, hi = float(metadata["wavelength_min_nm"]), float(metadata["wavelength_max_nm"])
        # One wavelength (monochromatic light sources): the 1/lambda^2 law over
        # [lo, lo + 1) truncates every draw to floor(lo), as sim_telarray does.
        return (lo, lo + 1.0) if hi == lo else (lo, hi)
    return None


def _int_key(metadata: Mapping[str, Any], key: str, record: str) -> int:
    if key not in metadata:
        raise ModelError(f"{record} record without {key}")
    return int(metadata[key])


def _camera_columns(
    traced: dict[str, np.ndarray], input_columns: Mapping[str, np.ndarray], config: _Config
) -> dict[str, np.ndarray]:
    """Columns of the C-frame photons from the output of `trace_bunches`."""
    weight = traced["weight"]
    if config.envelopes.tables:
        weight = weight * config.envelopes.inverse(traced["wavelength_nm"])
    position, direction = traced["position_m"], traced["direction"]
    columns = {
        "x": position[:, 0],
        "y": position[:, 1],
        "cx": direction[:, 0],
        "cy": direction[:, 1],
        "cz": direction[:, 2],
        "t": traced["time_ns"],
        "w": weight,
        "wavelength": traced["wavelength_nm"],
    }
    if "zem" in input_columns:
        columns["zem"] = input_columns["zem"][traced["bunch_index"]]
    return {
        name: np.ascontiguousarray(values, dtype=np.float32) for name, values in columns.items()
    }


class _Tracer:
    """Reads the records of one stream and writes the traced records.

    ``optics`` is one ``_Optics`` for every telescope, or a dict of them by telescope id. With a
    dict, the factors that a model applies are marked in an efficiency status at the start of
    each telescope block (the provenance differs between models), not at run level.
    """

    def __init__(self, optics, writer, *, seed, pool, telescopes, pointing):
        self.optics = optics
        self.per_telescope = isinstance(optics, dict)
        self.writer = writer
        self.seed = seed
        self.pool = pool
        self.telescopes = telescopes
        self.pointing = pointing
        self.run_pointing: tuple[float, float] | None = None
        self.run_factors: dict[str, _Factor] = {}
        self.law: tuple[float, float] | None = None
        self.event: dict[str, Any] | None = None
        # Ground-frame rotation and input plane of the event, by telescope id.
        self.geometry: dict[int, tuple[np.ndarray, float]] = {}
        self.telescope: _Telescope | None = None
        # Row indices continue across repeated blocks of one telescope in an event.
        self.next_index: dict[int, int] = {}
        self.counts = {"events": 0, "photons_in": 0, "photons_out": 0}
        self.run_status_pending = False

    def _optics(self, telescope_id: int) -> _Optics:
        if not self.per_telescope:
            return self.optics
        return self.optics[telescope_id]

    def _all_optics(self) -> list[_Optics]:
        return list(self.optics.values()) if self.per_telescope else [self.optics]

    def handle(self, record: dict[str, Any]) -> None:
        kind = record["type"]
        if self.run_status_pending and kind in ("event", "run_end"):
            self.writer.write(
                "efficiency_status",
                _own_efficiency_status(self.optics.model, self.optics.transmission),
            )
            self.run_status_pending = False
        if kind == "run":
            self._run(record["metadata"])
        elif kind == "efficiency_status":
            self._efficiency_status(record["metadata"])
        elif kind == "event":
            self._event(record["metadata"])
        elif kind == "telescope":
            self._telescope(record["metadata"])
        elif kind == "photons":
            self._photons(record)
        elif kind == "telescope_end":
            self._telescope_end(record["metadata"])
        elif kind == "blob":
            self.writer.write_blob(record["eventio_type"], record["eventio_id"], record["data"])
        elif kind in ("event_end", "run_end"):
            self.writer.write(kind, record["metadata"])
        elif kind != "stream_end":
            raise ModelError(f"unexpected photon-stream record {kind!r}")

    def summary(self) -> dict[str, int]:
        return dict(self.counts)

    def _run(self, metadata: Mapping[str, Any]) -> None:
        self.run_factors = {}
        self.law = _wavelength_law(metadata)
        self.run_pointing = _run_pointing(metadata)
        self.event = None
        self.writer.write("run", metadata)
        self.run_status_pending = not self.per_telescope

    def _efficiency_status(self, metadata: Mapping[str, Any]) -> None:
        if self.telescope is None:
            self.run_factors = {**self.run_factors, **_efficiency_status_factors(metadata)}
            if any(optics.transmission.normalise for optics in self._all_optics()) and (
                CONSTANT_FACTOR in _names(metadata.get("factors", ""))
            ):
                raise ModelError(f"the efficiency status already names {CONSTANT_FACTOR}")
            self.writer.write("efficiency_status", metadata)
            return
        else:
            telescope = self.telescope
            if telescope.photons_seen:
                raise ModelError(
                    f"telescope {telescope.telescope_id}: efficiency_status after its photons"
                )
            if not telescope.kept:
                return
            telescope.factors = {**telescope.factors, **_efficiency_status_factors(metadata)}
            if self.per_telescope:
                if telescope.status_written:
                    raise ModelError(
                        f"telescope {telescope.telescope_id}: more than one efficiency_status "
                        "before its photons"
                    )
                self._write_status(telescope, metadata)
                return
        optics = self._optics(telescope.telescope_id)
        self.writer.write(
            "efficiency_status",
            _efficiency_status_with_applied(metadata, optics.model, optics.transmission),
        )

    def _event(self, metadata: Mapping[str, Any]) -> None:
        _int_key(metadata, "event_number", "event")
        self.event = dict(metadata)
        self.geometry = {}
        self.telescope = None
        self.next_index = {}
        self.counts["events"] += 1
        self.writer.write("event", metadata)

    def _telescope(self, metadata: Mapping[str, Any]) -> None:
        telescope_id = _int_key(metadata, "telescope_id", "telescope")
        kept = self.telescopes is None or telescope_id in self.telescopes
        start = self.next_index.get(telescope_id, 0)
        self.telescope = _Telescope(telescope_id, kept, out_index=start)
        if kept:
            self.writer.write("telescope", {**metadata, "frame": FRAME_CAMERA})

    def _write_status(self, telescope: _Telescope, metadata: Mapping[str, Any] | None) -> None:
        """Writes the telescope's efficiency status, with the factors obdeect applies applied.

        With a mapping of models, this is the one record of the block that marks them applied:
        an incoming status is rewritten, or a record of the own factors is written before the
        first photons, so no factor is applied twice in one telescope scope.
        """
        optics = self._optics(telescope.telescope_id)
        self.writer.write(
            "efficiency_status",
            _own_efficiency_status(optics.model, optics.transmission, metadata),
        )
        telescope.status_written = True

    def _ensure_status(self, telescope: _Telescope) -> None:
        if self.per_telescope and not telescope.status_written:
            self._write_status(telescope, None)

    def _telescope_end(self, metadata: Mapping[str, Any]) -> None:
        telescope = self.telescope
        if telescope is None or _int_key(metadata, "telescope_id", "telescope_end") != (
            telescope.telescope_id
        ):
            raise ModelError("telescope_end record without its telescope record")
        if telescope.kept:
            # Applies the efficiency status checks even for a telescope without photons.
            self._ensure_status(telescope)
            self._config(telescope)
            self.writer.write("telescope_end", {**metadata, "photons": telescope.weight})
            self.next_index[telescope.telescope_id] = telescope.out_index
        self.telescope = None

    def _config(self, telescope: _Telescope) -> _Config:
        if telescope.config is None:
            factors = {**self.run_factors, **telescope.factors}
            tables = []
            for name in _own_factors(self._optics(telescope.telescope_id).model):
                factor = factors.get(name)
                if factor is None or factor.state == "pending":
                    continue
                if factor.state == "applied":
                    raise ModelError(
                        f"telescope {telescope.telescope_id}: the efficiency status says {name} is "
                        "already applied; trace-stream applies it in its trace and cannot skip it"
                    )
                tables.append(_envelope_table(name, factor))
            telescope.config = _Config(roulette=not tables, envelopes=_Envelopes(tuple(tables)))
        return telescope.config

    def _ground_geometry(self, telescope_id: int) -> tuple[np.ndarray, float]:
        if telescope_id not in self.geometry:
            zenith_deg, azimuth_deg = self._pointing()
            model = self._optics(telescope_id).model
            self.geometry[telescope_id] = (
                ground_to_telescope_rotation(zenith_deg, azimuth_deg),
                minimum_input_plane_height_m(model, zenith_deg, azimuth_deg),
            )
        return self.geometry[telescope_id]

    def _pointing(self) -> tuple[float, float]:
        if self.pointing is not None:
            return self.pointing
        if self.run_pointing is not None:
            return self.run_pointing
        raise ModelError(
            "no telescope pointing for ground-frame photons: give --pointing-zenith-deg and "
            "--pointing-azimuth-deg, or the run keys telescope_pointing_zenith_rad and "
            "telescope_pointing_azimuth_rad"
        )

    def _photons(self, record: dict[str, Any]) -> None:
        telescope = self.telescope
        if telescope is None or telescope.telescope_id != record["telescope_id"]:
            raise ModelError(f"photons of telescope {record['telescope_id']} outside their block")
        if not telescope.kept:
            return
        if self.event is None:
            raise ModelError("photons outside an event record")
        telescope.photons_seen = True
        columns = record["columns"]
        if "p" in columns:
            raise ModelError("photon column p (acceptance already applied) is not supported")
        unsupported = set(columns) - {
            "x",
            "y",
            "z",
            "cx",
            "cy",
            "cz",
            "t",
            "w",
            "wavelength",
            "zem",
        }
        if unsupported:
            raise ModelError(f"unsupported photon columns: {', '.join(sorted(unsupported))}")
        rows = len(columns["w"])
        if rows == 0:
            return
        if record["frame"] not in (FRAME_GROUND, FRAME_TELESCOPE):
            raise ModelError(
                f"photon frame {record['frame']} is not supported: trace-stream takes ground "
                "(1) or telescope (2) photons"
            )
        self._ensure_status(telescope)
        config = self._config(telescope)
        optics = self._optics(telescope.telescope_id)
        event_number = _int_key(self.event, "event_number", "event")
        bunches, placement = self._input_bunches(record, event_number, optics)
        traced = optics.model.trace_bunches(
            bunches,
            seed=self.seed,
            event=event_number,
            telescope=telescope.telescope_id,
            length_unit_m=1.0,
            wavelength_law=self.law,
            roulette=config.roulette,
            normalise_transmission=optics.transmission.normalise,
            pool=self.pool,
            **placement,
        )
        self.counts["photons_in"] += rows
        columns = _camera_columns(traced, columns, config)
        self._write(telescope, columns)

    def _input_bunches(
        self, record: dict[str, Any], event_number: int, optics: _Optics
    ) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
        """Bunches on the input plane and the frame that places it (see trace_bunches)."""
        columns = record["columns"]
        rows = len(columns["w"])
        cx = columns["cx"].astype(np.float64)
        cy = columns["cy"].astype(np.float64)
        cz = -np.sqrt(np.maximum(0.0, 1.0 - cx * cx - cy * cy))
        if (
            "cz" in columns
            and np.max(np.abs(columns["cz"].astype(np.float64) - cz), initial=0.0) > CZ_TOLERANCE
        ):
            raise ModelError("photon column cz does not match cx and cy")
        if not np.all(cz < 0.0):
            raise ModelError("photons must move towards -z (downwards)")
        x = columns["x"].astype(np.float64)
        y = columns["y"].astype(np.float64)
        t = columns["t"].astype(np.float64)
        if record["frame"] == FRAME_GROUND:
            if "z" in columns and np.any(columns["z"] != 0):
                raise ModelError("ground-frame photons must lie in the plane z = 0")
            z = np.zeros(rows)
            rotation, plane = self._ground_geometry(record["telescope_id"])
            offset = np.zeros(3)
        else:
            if "z" not in columns:
                raise ModelError("telescope-frame photons need the column z")
            z = columns["z"].astype(np.float64)
            rotation = np.eye(3)
            plane = optics.telescope_plane_m
            # The telescope frame is the frame of trace(): cancel the axes offset that
            # trace_bunches subtracts, so that a photon at p is traced at p.
            offset = np.asarray(optics.model.axes_offset_m, dtype=np.float64)
        # Move each photon along its direction to the input plane; the time follows the path.
        s = (plane - z) / cz
        first_index = int(record["first_index"])
        if event_number >= 2**32 or first_index + rows >= 2**32:
            raise ModelError("bunch ids need event numbers and photon indices below 2^32")
        index = np.arange(first_index, first_index + rows, dtype=np.uint64)
        bunch_id = (np.uint64(event_number) << np.uint64(32)) | index
        float32 = np.float32
        bunches = {
            "x": (x + s * cx).astype(float32),
            "y": (y + s * cy).astype(float32),
            "cx": cx.astype(float32),
            "cy": cy.astype(float32),
            "time_ns": (t + s / SPEED_OF_LIGHT_M_PER_NS).astype(float32),
            "zem": np.asarray(columns.get("zem", np.zeros(rows)), dtype=float32),
            "size": np.asarray(columns["w"], dtype=float32),
            "wavelength_nm": np.asarray(columns.get("wavelength", np.zeros(rows)), dtype=float32),
            "bunch_id": bunch_id,
        }
        placement = {"frame_rotation": rotation, "frame_offset_m": offset, "plane_z": plane}
        return bunches, placement

    def _write(self, telescope: _Telescope, columns: dict[str, np.ndarray]) -> None:
        rows = len(columns["x"])
        for start in range(0, rows, MAX_CHUNK_ROWS):
            stop = min(rows, start + MAX_CHUNK_ROWS)
            chunk = {name: values[start:stop] for name, values in columns.items()}
            self.writer.write_photons(telescope.telescope_id, "camera", telescope.out_index, chunk)
            telescope.out_index += stop - start
        telescope.weight += float(np.sum(columns["w"], dtype=np.float64))
        self.counts["photons_out"] += rows


def _writer_path(output: Any) -> str:
    # The stream writer opens files only; "-" writes to the process's standard output.
    return "/dev/stdout" if str(output) == "-" else str(output)


def _optics_of(
    model: OpticalModel | Mapping[int, OpticalModel],
    telescopes: frozenset[int] | None,
    part: str,
) -> tuple[_Optics | dict[int, _Optics], frozenset[int] | None]:
    """The optics of the tracer and the telescopes to trace (see trace_stream)."""
    if isinstance(model, OpticalModel):
        return _Optics.of(model, part), telescopes
    models = {_telescope_id(int(telescope)): value for telescope, value in model.items()}
    if not models:
        raise ModelError("give at least one optical model")
    if telescopes is None:
        telescopes = frozenset(models)
    without_model = sorted(telescopes - set(models))
    if without_model:
        raise ModelError(f"telescopes without an optical model: {without_model}")
    # One _Optics per distinct model: telescopes of one type share it.
    shared: dict[int, _Optics] = {}
    optics = {}
    for telescope in sorted(telescopes):
        value = models[telescope]
        if id(value) not in shared:
            shared[id(value)] = _Optics.of(value, part)
        optics[telescope] = shared[id(value)]
    return optics, telescopes


def trace_stream(
    model: OpticalModel | Mapping[int, OpticalModel],
    input: Any,
    output: Any,
    *,
    seed: int,
    threads: int = 0,
    telescopes: Iterable[int] | None = None,
    pointing: tuple[float, float] | None = None,
    telescope_transmission: str = "full",
) -> dict[str, int]:
    """Trace the photons of a photon stream and write the camera-frame photons.

    Parameters
    ----------
    model: OpticalModel or mapping of int to OpticalModel
        Optical model of every traced telescope, or the optical model of each telescope id.
        With a mapping, telescopes without a model are dropped, and the efficiency status of
        each telescope block names the model that traced it.
    input, output: str or path-like
        Photon-stream files; "-" is standard input or output.
    seed: int
        Run seed (non-negative). Draws are keyed on (seed, event, telescope, photon), as in the
        other components.
    threads: int
        Worker threads (0 uses every core); results do not depend on it.
    telescopes: iterable of int or None
        Telescope ids to trace; other telescopes are dropped. None traces all (with a mapping
        of models: all telescopes with a model). Each id needs a model.
    pointing: (float, float) or None
        Telescope (zenith, azimuth) [deg]; azimuth in the array frame (x north,
        y west), as simtools' telescope phi. None takes it from the run keys
        telescope_pointing_zenith_rad and telescope_pointing_azimuth_rad [rad].
        Ground-frame photons need one of the two; the event header has no pointing.
    telescope_transmission: "full" or "angular"
        "full" applies the telescope transmission t(theta), including its on-axis
        value t0. "angular" applies t(theta) / t0 only; the efficiency status then leaves
        t0 as the pending factor ``telescope_transmission_constant`` for the
        downstream camera module. Models without telescope transmission reject
        "angular".

    Returns
    -------
    dict
        Counts of ``events``, ``photons_in`` and ``photons_out``.

    Raises
    ------
    ModelError
        For an invalid stream, a missing pointing (ground-frame photons only), an
        unsupported frame or column,
        or an efficiency status that says a factor is already applied, and for an unknown
        telescope_transmission or "angular" with a model that has none, and for a telescope
        id without a model or outside 1..65535.
    """
    if seed < 0:
        raise ModelError("seed must not be negative")
    if telescopes is not None:
        telescopes = frozenset(_telescope_id(int(telescope)) for telescope in telescopes)
        if not telescopes:
            raise ModelError("telescopes: give at least one telescope id")
    optics, telescopes = _optics_of(model, telescopes, telescope_transmission)
    if pointing is not None:
        pointing = (float(pointing[0]), float(pointing[1]))
    try:
        reject_same_file(input, output)
    except ValueError as error:
        raise ModelError(str(error)) from error
    # Photon chunks of other telescopes are skipped by the reader, not decoded.
    kept = None if telescopes is None else sorted(telescopes)
    records = iter(obdeect_stream.Reader(str(input), kept))
    head = next(records, None)
    if head is None or head["type"] != "stream":
        raise ModelError(f"{input}: not a photon stream")
    producer = str(head["metadata"].get("producer", ""))
    writer = obdeect_stream.Writer(
        _writer_path(output),
        f"{producer}; obdeect {__version__}" if producer else f"obdeect {__version__}",
    )
    tracer = _Tracer(
        optics,
        writer,
        seed=seed,
        pool=ThreadPool(threads),
        telescopes=telescopes,
        pointing=pointing,
    )
    for record in records:
        tracer.handle(record)
    # Written only after the whole stream: a failed run leaves no stream_end record.
    writer.close()
    return tracer.summary()
