"""Optical models: build, store and trace."""

from __future__ import annotations

import json
import math
from collections.abc import Callable, Iterator, Mapping, Sequence
from copy import deepcopy
from functools import cached_property
from pathlib import Path
from typing import Any

import numpy as np
from obdeect.ray_tracing import _core
from obdeect.ray_tracing.bounds import UP as _UP
from obdeect.ray_tracing.bounds import mirror_support as _mirror_support
from obdeect.ray_tracing.bounds import obscurer_support as _obscurer_support
from obdeect.ray_tracing.bounds import obscurer_top_z_m as _obscurer_top_z_m
from obdeect.ray_tracing.errors import ModelError
from obdeect.ray_tracing.hashing import model_hash
from obdeect.ray_tracing.importers.simulation_models import compile_model, compile_telescope

STATUS_NAMES: tuple[str, ...] = _core.STATUS_NAMES
DETECTED = STATUS_NAMES.index("detected")
# Path segments of obscurers (obdeect::ray_tracing::Leg bits).
_LEGS = {"incoming": 1, "after_primary": 2, "after_secondary": 4}
ThreadPool = _core.ThreadPool


class OpticalModel:
    """A compiled optical model, ready for tracing.

    The model is immutable; its SHA-256 identifies it in every result.
    """

    def __init__(self, document: dict[str, Any]):
        isolated_document = deepcopy(document)
        if isolated_document.get("sha256") != model_hash(isolated_document):
            raise ModelError("optical model content does not match its sha256")
        try:
            self._native = _core.Model(isolated_document)
        except ValueError as error:
            raise ModelError(str(error)) from error
        self._document = isolated_document
        # Thread pools by thread count, made on first use and shared by the calls on this model.
        self._pools: dict[int, ThreadPool] = {}
        # Results of mirrors_support and optics_bottoms by direction (the model is immutable).
        self._support_cache: dict[tuple[float, ...], list[tuple[str, float]]] = {}
        self._bottoms_cache: dict[tuple[float, ...], list[tuple[str, float]]] = {}

    @property
    def document(self) -> dict[str, Any]:
        """Return a copy of the compiled document (internal code reads `_document` directly)."""
        return deepcopy(self._document)

    @classmethod
    def from_records(
        cls,
        records: Mapping[str, dict[str, Any]],
        read_file: Callable[[str], str],
        *,
        configuration_seed: int,
        propagation_group_index: float,
        zenith_deg: float | None = None,
        provenance: dict[str, Any] | None = None,
    ) -> OpticalModel:
        """Compile structured simulation-model parameter records into an optical model.

        The importer (`obdeect.ray_tracing.importers.simulation_models`) defines the parameter
        names; `propagation_group_index` is required (no default air index).
        """
        document = compile_model(
            records,
            read_file,
            configuration_seed=configuration_seed,
            zenith_deg=zenith_deg,
            propagation_group_index=propagation_group_index,
            provenance=provenance,
        )
        return cls(document)

    @classmethod
    def from_simulation_models(
        cls,
        root: Path,
        model: str,
        version: str,
        *,
        configuration_seed: int,
        zenith_deg: float | None = None,
        propagation_group_index: float | None = None,
    ) -> OpticalModel:
        """Compile an optical model from a simulation-models checkout (see `compile_telescope`)."""
        return cls(
            compile_telescope(
                root,
                model,
                version,
                configuration_seed=configuration_seed,
                zenith_deg=zenith_deg,
                propagation_group_index=propagation_group_index,
            )
        )

    @classmethod
    def load(cls, path: Path) -> OpticalModel:
        return cls(json.loads(Path(path).read_text(encoding="utf-8")))

    def save(self, path: Path) -> None:
        Path(path).write_text(
            json.dumps(self._document, ensure_ascii=False, indent=1) + "\n", encoding="utf-8"
        )

    @property
    def sha256(self) -> str:
        return self._document["sha256"]

    @property
    def has_structure(self) -> bool:
        """True if the model's own geometry shadows instead of ``telescope_transmission``."""
        return self._document["provenance"]["structure"]

    @property
    def sampling_area_m2(self) -> float:
        """Area of the disk on which star photons are sampled."""
        return math.pi * self._document["sampling"]["disk_radius_m"] ** 2

    def _shared_pool(self, threads: int) -> ThreadPool:
        """The thread pool of this model with `threads` threads (0: every core), made once."""
        if threads not in self._pools:
            self._pools[threads] = ThreadPool(threads)
        return self._pools[threads]

    @property
    def focal_surface_z_m(self) -> float:
        """Height of the focal plane, or of the focal-surface vertex, in the telescope frame."""
        focal = self._document["optics"]["focal_surface"]
        return focal["asphere"]["vertex_z_m"] if "asphere" in focal else focal["plane_z_m"]

    @property
    def axes_offset_m(self) -> tuple[float, float, float]:
        """Offset O [m, telescope frame] of the optical origin from the telescope sphere centre.

        A ground-frame point P is traced at R P - O in the telescope frame, with R the input
        rotation of `trace_bunches`, so the dish sits at +O relative to the rotation axes.
        Zero for models without an axes offset.
        """
        x, y, z = self._document["optics"].get("axes_offset_m", (0.0, 0.0, 0.0))
        return (float(x), float(y), float(z))

    @property
    def propagation_group_index(self) -> float:
        """Group index of air the photon times use."""
        return self._document["optics"]["propagation_group_index"]

    def max_incoming_obscurer_z_m(self) -> float | None:
        """Highest z [m, telescope frame] of any obscurer on the incoming leg, or None.

        Each obscurer is bounded from above: a plate by its centre plus its aperture's
        circumradius, a cylinder or frustum by its end points plus the outer end radii
        (wall thickness included), a polyhedron by its highest vertex. Photons
        starting below this height never meet those obscurers.
        """
        top = self._incoming_top
        return None if top is None else top[0]

    def max_mirror_z_m(self) -> float | None:
        """Highest z [m, telescope frame] of the mirrors that photons meet first, or None.

        The primary, and the secondary if it blocks incoming light (see `mirrors_support`).
        """
        tops = [support for _, support in self.mirrors_support(_UP)]
        return max(tops) if tops else None

    def mirrors_support(self, direction: np.ndarray) -> list[tuple[str, float]]:
        """Maximum dot product with a direction over the mirrors that photons meet first [m].

        Photons start on the input plane and meet the primary first, then the secondary
        after it: only the primary and, if it blocks incoming light, the secondary must lie
        below the plane.
        """
        key = tuple(float(x) for x in np.ravel(direction))
        if key not in self._support_cache:
            optics = self._document["optics"]
            roles = ["primary"]
            if optics.get("secondary_blocks_incoming"):
                roles.append("secondary")
            self._support_cache[key] = [
                (role, _mirror_support(optics[role], direction))
                for role in roles
                if optics.get(role) is not None
            ]
        return list(self._support_cache[key])

    def optics_bottoms(self, normal: np.ndarray) -> list[tuple[str, float]]:
        """Lowest point along `normal` of each obscurer (any leg) and of the mirrors [m].

        The mirrors are the primary and the secondary (photons from below meet both). The
        lowest point of an optic is min n.p, i.e. minus its support in direction -n.
        """
        key = tuple(float(x) for x in np.ravel(normal))
        if key not in self._bottoms_cache:
            optics = self._document["optics"]
            flipped = -np.asarray(normal, dtype=float)
            bottoms = [(o["name"], -_obscurer_support(o, flipped)) for o in optics["obscurers"]]
            for role in ("primary", "secondary"):
                if optics.get(role) is not None:
                    bottoms.append((role, -_mirror_support(optics[role], flipped)))
            self._bottoms_cache[key] = bottoms
        return list(self._bottoms_cache[key])

    @cached_property
    def _incoming_top(self) -> tuple[float, str] | None:
        """Highest incoming obscurer top [m] and its name, or None without incoming obscurers."""
        tops = [
            (_obscurer_top_z_m(o), o["name"])
            for o in self._document["optics"]["obscurers"]
            if "incoming" in o.get("legs", [])
        ]
        return max(tops) if tops else None

    def emitter_photons(
        self,
        emitted: float,
        *,
        seed: int,
        position_m: tuple[float, float, float],
        intensity: tuple[np.ndarray, np.ndarray],
        axis: tuple[float, float, float] = (0.0, 0.0, 1.0),
        pulse: tuple[np.ndarray, np.ndarray] | None = None,
        wavelength_nm: float = 400.0,
        disk: tuple[tuple[float, float, float], tuple[float, float, float], float] | None = None,
        first_id: int = 0,
    ) -> dict[str, np.ndarray]:
        """Photons of a point emitter (a calibration light source) that cross a disk.

        The emitter at `position_m` (telescope frame) sends `emitted` photons
        into 4 pi with the relative intensity per solid angle `intensity`
        (angles from `axis` [deg], values; linear in 1 - cos(angle), zero
        beyond the last angle) and emission times from the relative density
        `pulse` (times [ns], values; None: all at 0). Photons are made only
        on their way through `disk` (centre, unit normal, radius [m]; default
        the star-light disk). Their weights are in [0, 1], and their number is
        Poisson distributed, so keeping each photon with the probability of
        its weight gives the photons crossing the disk; the sum of weights
        estimates their mean number. Start positions are the emitter, times
        are emission times.
        """
        source = self._emitter(position_m, intensity, axis, pulse, wavelength_nm, disk)
        return source.sample(first_id, _emitter_count(source, emitted, seed, first_id), seed)

    def trace_emitter(self, emitted: float, **options: Any) -> dict[str, np.ndarray]:
        """Make and trace the photons of a point emitter (see `trace_emitter_blocks`).

        Returns the photons of all blocks reaching the focal surface with a
        positive weight.
        """
        blocks = list(self.trace_emitter_blocks(emitted, **options))
        return {name: np.concatenate([block[name] for block in blocks]) for name in blocks[0]}

    def trace_emitter_blocks(
        self,
        emitted: float,
        *,
        seed: int,
        position_m: tuple[float, float, float],
        intensity: tuple[np.ndarray, np.ndarray],
        axis: tuple[float, float, float] = (0.0, 0.0, 1.0),
        pulse: tuple[np.ndarray, np.ndarray] | None = None,
        wavelength_nm: float = 400.0,
        disk: tuple[tuple[float, float, float], tuple[float, float, float], float] | None = None,
        first_id: int = 0,
        direct: bool = False,
        threads: int = 0,
        block_size: int = 250_000,
    ) -> Iterator[dict[str, np.ndarray]]:
        """Make and trace the photons of a point emitter, `block_size` photons at a time.

        The photons are those of `emitter_photons`, traced as by `trace`;
        memory is bounded by one block. Yields, for each block (at least
        one), the photons reaching the focal surface with a positive weight:
        position_m, direction, time_ns and weight of the trace, and
        wavelength_nm. The photons do not depend on `block_size`.
        """
        if block_size < 1:
            raise ModelError("block_size must be positive")
        source = self._emitter(position_m, intensity, axis, pulse, wavelength_nm, disk)
        count = _emitter_count(source, emitted, seed, first_id)
        for start in range(0, count, block_size) or [0]:
            photons = source.sample(first_id + start, min(block_size, count - start), seed)
            traced = self.trace(photons, seed=seed, direct=direct, threads=threads)
            kept = (traced["status"] == DETECTED) & (traced["weight"] > 0.0)
            yield {name: traced[name][kept] for name in ("position_m", "direction", "time_ns")} | {
                "weight": traced["weight"][kept],
                "wavelength_nm": photons["wavelength_nm"][kept],
            }

    def _emitter(self, position_m, intensity, axis, pulse, wavelength_nm, disk) -> _core.Emitter:
        if disk is None:
            sampling = self._document["sampling"]
            disk = ((0.0, 0.0, sampling["disk_z_m"]), (0.0, 0.0, 1.0), sampling["disk_radius_m"])
        angle_deg, values = (np.asarray(column, dtype=float) for column in intensity)
        pulse_ns, pulse_values = (
            (np.empty(0), np.empty(0))
            if pulse is None
            else (np.asarray(column, dtype=float) for column in pulse)
        )
        centre, normal, radius = disk
        return _core.Emitter(
            np.asarray(position_m, dtype=float),
            np.asarray(axis, dtype=float),
            1.0 - np.cos(np.radians(angle_deg)),
            values,
            pulse_ns,
            pulse_values,
            np.asarray(centre, dtype=float),
            np.asarray(normal, dtype=float),
            float(radius),
            wavelength_nm,
        )

    def star_photons(
        self,
        count: int,
        *,
        seed: int,
        towards_source: tuple[float, float, float] = (0.0, 0.0, 1.0),
        distance_m: float = math.inf,
        wavelength_nm: float = 400.0,
        first_id: int = 0,
    ) -> dict[str, np.ndarray]:
        """Sample photons from a star or a point source.

        `towards_source` is the unit vector to the source in the telescope frame
        (+z along the optical axis, towards the sky; a negative z gives a source below
        the telescope and photons moving up, see docs/ray-tracing/light-sources.md);
        `distance_m` is the source distance from the dish vertex, the centre of the
        sampled disk.
        Times are zero for the wavefront at that point.
        """
        sampling = self._document["sampling"]
        direction = np.asarray(towards_source, dtype=float)
        direction = direction / np.linalg.norm(direction)
        start_z_m = sampling["start_z_m"]
        if direction[2] < 0.0:
            # Photons moving up start on a plane below every optic, not on the plane above them.
            bottom = min((b for _, b in self.optics_bottoms(_UP)), default=0.0)
            start_z_m = bottom - _LOWER_START_MARGIN_M
        return _core.sample_source(
            direction,
            distance_m,
            sampling["disk_radius_m"],
            sampling["disk_z_m"],
            start_z_m,
            wavelength_nm,
            first_id,
            count,
            seed,
        )

    def trace(
        self,
        photons: Mapping[str, np.ndarray],
        *,
        seed: int,
        apply_responses: bool = True,
        threads: int = 0,
        direct: bool = False,
        outputs: Sequence[str] | None = None,
        detected_only: bool = False,
        classify_losses: bool = True,
        pool: ThreadPool | None = None,
    ) -> dict[str, np.ndarray]:
        """Trace photons; returns per-photon arrays (see docs/ray-tracing/python-guide.md).

        `photons` holds position_m (n, 3), direction (n, 3), wavelength_nm,
        time_ns, weight and photon_id. With apply_responses False only geometry
        acts (weights unchanged); `threads` 0 uses every core. `direct`
        photons start inside the telescope on their way to the focal surface
        (e.g. from a flasher): mirrors and telescope transmission do not act,
        the obscurers of that leg do.

        `outputs` lists the columns to return (None: all of them); asking only
        for what is used saves memory traffic. With `detected_only` just the
        photons with status ``detected`` are returned, in input order, and the
        extra column ``index`` (their row in the input) becomes available.
        An unknown column raises ``ValueError``.

        With `classify_losses` False a photon that misses the primary mirror
        is reported as ``missed_primary`` even where it would be ``blocked``
        before missing it; detected photons are the same either way.

        With a `pool` (a `ThreadPool`) its threads run the work and `threads`
        is ignored; the results are the same.
        """
        arrays = _photon_arrays(photons)
        return self._native.trace(
            arrays["position_m"],
            arrays["direction"],
            arrays["wavelength_nm"],
            arrays["time_ns"],
            arrays["weight"],
            arrays["photon_id"],
            seed=seed,
            apply_responses=apply_responses,
            threads=threads,
            direct=direct,
            outputs=None if outputs is None else list(outputs),
            detected_only=detected_only,
            classify_losses=classify_losses,
            pool=pool,
        )

    def obscurer_hits(
        self, photons: Mapping[str, np.ndarray], *, seed: int, threads: int = 0
    ) -> dict[str, Any]:
        """Which obscurers shade the light: a diagnostic, slower than `trace`.

        The photons are traced without obscurers (geometry only, mirror
        scattering from `seed` as in `trace`); along the path of each one that
        then reaches the focal surface, every obscurer it meets is listed. Per
        obscurer, in model order, the result holds ``id``, ``name``, ``group``
        and the summed input weights of the photons whose path ``crossed`` it,
        that it absorbs (``first``: the first obscurer on the path) and that
        only it meets (``sole``: without it they would be detected).
        ``unobscured`` and ``detected`` are the summed weights of the photons
        reaching the focal surface without and with the obscurers. The sequential
        photons (one reflection on the primary of a single mirror, one on each
        mirror of a dual one) are followed along their path; any other photon
        counts its terminal obscurer as crossed.
        """
        free_document = self.document
        obscurers = free_document["optics"]["obscurers"]
        free_document["optics"]["obscurers"] = []
        free_document["sha256"] = model_hash(free_document)
        free = OpticalModel(free_document).trace(
            photons,
            seed=seed,
            apply_responses=False,
            threads=threads,
            outputs=["primary_point_m", "secondary_point_m", "position_m", "index", "reflections"],
            detected_only=True,
        )
        full = self.trace(
            photons, seed=seed, apply_responses=False, threads=threads, outputs=["terminal_surface"]
        )
        index = free["index"]
        n = len(index)
        weight = np.asarray(photons["weight"], dtype=float)[index]
        terminal = full["terminal_surface"][index]
        id_list = [o["id"] for o in obscurers]
        absorbed = np.isin(terminal, id_list)
        # The path of a sequential photon is the start, the primary point, the secondary point and
        # the focal point: one reflection on the primary of a single mirror, one on each mirror of
        # a dual one. The other photons (multiple or direct) count their terminal obscurer only.
        expected = 2 if free_document["optics"].get("secondary") else 1
        sequential = np.flatnonzero(np.asarray(free["reflections"]) == expected)
        start = np.asarray(photons["position_m"], dtype=float)[index[sequential]]
        points = [start, free["primary_point_m"][sequential]]
        legs = [_LEGS["incoming"], _LEGS["after_primary"]]
        if free_document["optics"].get("secondary"):
            points.append(free["secondary_point_m"][sequential])
            legs.append(_LEGS["after_secondary"])
        points.append(free["position_m"][sequential])
        m = len(sequential)
        rows, ids = self._native.obscurer_crossings(
            np.concatenate(points[:-1]),
            np.concatenate(points[1:]),
            np.repeat(np.array(legs, dtype=np.uint8), m),
        )
        crossing_photon = sequential[(rows % m if m else rows).astype(np.int64)]
        crossing_obscurer = _positions(ids, id_list)
        terminal_photon = np.flatnonzero(absorbed)
        terminal_obscurer = _positions(terminal[absorbed], id_list)
        # One entry per photon and obscurer, also if the path meets it twice.
        pairs = (
            np.unique(
                np.column_stack([
                    np.concatenate([crossing_photon, terminal_photon]).astype(np.int64),
                    np.concatenate([crossing_obscurer, terminal_obscurer]),
                ]),
                axis=0,
            )
            if n
            else np.empty((0, 2), dtype=np.int64)
        )
        crossings_per_photon = np.bincount(pairs[:, 0], minlength=n)
        sole = pairs[crossings_per_photon[pairs[:, 0]] == 1]
        size = len(obscurers)
        return {
            "id": np.array(id_list, dtype=np.uint32),
            "name": np.array([o["name"] for o in obscurers]),
            "group": np.array([o.get("group", "") for o in obscurers]),
            "crossed": np.bincount(pairs[:, 1], weights=weight[pairs[:, 0]], minlength=size),
            "first": np.bincount(terminal_obscurer, weights=weight[absorbed], minlength=size),
            "sole": np.bincount(sole[:, 1], weights=weight[sole[:, 0]], minlength=size),
            "unobscured": float(weight.sum()),
            "detected": float(weight[~absorbed].sum()),
        }

    def trace_bunches(
        self,
        bunches: Mapping[str, np.ndarray],
        *,
        seed: int,
        event: int = 0,
        telescope: int = 0,
        frame_rotation: np.ndarray,
        frame_offset_m: Sequence[float],
        length_unit_m: float,
        plane_z: float,
        lower_plane_z: float | None = None,
        acceptance: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
        wavelength_density: tuple[np.ndarray, np.ndarray] | None = None,
        wavelength_law: tuple[float, float] | None = None,
        roulette: bool = True,
        classify_losses: bool = True,
        normalise_transmission: bool = False,
        pool: ThreadPool | None = None,
        threads: int = 0,
    ) -> dict[str, np.ndarray]:
        """Trace photon bunches of an air shower (see docs/ray-tracing/batch-tracing.md).

        `bunches` holds x, y, cx, cy, time_ns, zem, size, wavelength_nm (all
        float32; wavelength 0 means "sample one from `wavelength_density` or
        `wavelength_law`") and bunch_id (uint64, unique per bunch for one
        `seed`). The optional column cz (float32, either sign, never 0) gives the
        direction; without it the bunches move towards -z. With the column z
        (float32, input units; sphere rows, which need cz and may have cz = 0) every
        bunch starts at its own point (x, y, z) in the input frame, with time_ns there,
        and the input planes (`plane_z`, `lower_plane_z`) are not used. The optional column
        path (float64, m) is the path from the emission to each start point; the result then
        gives path_m, the path to each output point (the input path plus the traced path).
        The random draws are keyed
        on (`seed`, `event`, `telescope`, bunch, sub-photon), as in the other components.
        `frame_rotation` (3 x 3, input to telescope frame), `frame_offset_m`,
        `length_unit_m` (metres per input length unit) and `plane_z` (input units) place
        the input plane in the telescope frame; bunches moving up (cz > 0) start on
        `lower_plane_z` (input units), which must lie below every optic. The trace of a
        bunch whose direction in the telescope frame has z >= 0 uses the first-hit path
        (photons from behind the dish). The
        native tracer applies the model's `axes_offset_m` O_t (telescope frame): with zero
        `frame_offset_m`, a ground-frame point P is traced at R P - O_t = R (P - O_g),
        where O_g = R.T O_t is in the ground frame. This assumes
        `frame_rotation` is the ground-to-telescope rotation at the zenith angle the model
        was compiled for.
        `acceptance` is ``(wavelength_nm, zem_m, values)`` with values
        (n_wavelength, n_zem) in [0, 1]; `wavelength_density` is ``(nm,
        values)``. `wavelength_law` is ``(lo_nm, hi_nm)``: it samples
        wavelength-0 bunches, with density proportional to
        1/lambda^2 on [lo_nm, hi_nm] (1 <= lo_nm < hi_nm), continuous in nm. Give at
        most one of `wavelength_density` and `wavelength_law`.

        Each photon of a bunch is thinned by the acceptance and the model's
        response bound, traced, and, with `roulette`, kept with the probability
        min(1, weight / bound) (survivors have weight 1; without roulette they carry
        that capped ratio as weight). Returns the surviving photons as columns position_m,
        direction, time_ns, weight, wavelength_nm, bunch_index (input row) and
        sub_photon, in input order. Invalid input raises ``ValueError`` before
        any tracing. With `normalise_transmission` the roulette weight and
        bound are relative to the telescope transmission's on-axis value (the
        caller applies that value itself). The result does not depend on
        `threads` or `pool`.
        """
        dtypes = {
            "x": np.float32,
            "y": np.float32,
            "cx": np.float32,
            "cy": np.float32,
            "time_ns": np.float32,
            "zem": np.float32,
            "size": np.float32,
            "wavelength_nm": np.float32,
            "bunch_id": np.uint64,
        }
        arrays = {
            name: np.ascontiguousarray(bunches[name], dtype=dtype) for name, dtype in dtypes.items()
        }
        if "cz" in bunches:
            arrays["cz"] = np.ascontiguousarray(bunches["cz"], dtype=np.float32)
            if arrays["cz"].shape != arrays["x"].shape:
                raise ValueError("bunch arrays must have the same length")
        path = None
        if "path" in bunches:
            # Path [m] from the emission to each start point (output: path_m, input plus trace).
            path = np.ascontiguousarray(bunches["path"], dtype=np.float64)
            if path.shape != arrays["x"].shape:
                raise ValueError("bunch arrays must have the same length")
        if "z" in bunches:
            # Sphere rows: each bunch starts at its own height; the input planes do not apply.
            if "cz" not in bunches:
                raise ValueError("start heights z need the direction column cz")
            arrays["z"] = np.ascontiguousarray(bunches["z"], dtype=np.float32)
            if arrays["z"].shape != arrays["x"].shape:
                raise ValueError("bunch arrays must have the same length")
        rotation = np.ascontiguousarray(frame_rotation, dtype=np.float64)
        if rotation.shape != (3, 3):
            raise ValueError("frame_rotation must be a 3 x 3 matrix")
        offset = np.ascontiguousarray(frame_offset_m, dtype=np.float64)
        # The traced input origin is offset - axes_offset_m (applied natively).
        traced_offset = offset - np.asarray(self.axes_offset_m, dtype=np.float64)
        top = self._incoming_top
        if float(length_unit_m) > 0.0 and "z" not in arrays:
            # The input plane is {p : n.p = value} with n = R[:, 2]. Its height on the
            # telescope axis is value / n_z, and an optic with support s = max over its points
            # q of n.q lies entirely below the plane iff s <= value (for n_z > 0).
            normal = rotation[:, 2]
            point = traced_offset + normal * float(plane_z) * float(length_unit_m)
            value = float(np.dot(normal, point))
            if top is not None:
                if normal[2] > 0.0:
                    checks = [
                        (
                            value / float(normal[2]),
                            _obscurer_support(o, normal) / float(normal[2]),
                            o["name"],
                        )
                        for o in self._document["optics"]["obscurers"]
                        if "incoming" in o.get("legs", [])
                    ]
                else:
                    checks = [(-math.inf, top[0], top[1])]
                for plane_height_m, obscurer_top, name in checks:
                    if plane_height_m < obscurer_top:
                        raise ValueError(
                            f"input plane below incoming obscurer {name!r}: the plane is at "
                            f"{plane_height_m:.3f} m along the telescope axis, "
                            "but the obscurer reaches "
                            f"{obscurer_top:.3f} m; raise the input plane above it"
                        )
            # Photons from below start on the lower plane: it must lie below every optic.
            if lower_plane_z is not None and "cz" in arrays and bool(np.any(arrays["cz"] > 0)):
                lower_point = traced_offset + normal * float(lower_plane_z) * float(length_unit_m)
                lower_value = float(np.dot(normal, lower_point))
                for name, bottom in self.optics_bottoms(normal):
                    if lower_value > bottom:
                        raise ValueError(
                            f"input plane above the {name} for photons from below: the plane "
                            f"is at n.p = {lower_value:.3f} m but the optic reaches down to "
                            f"{bottom:.3f} m along the input normal; lower the input plane "
                            "below it"
                        )
            # Photons start on the plane and move on: a mirror point above it is never reached.
            for role, mirror_top in self.mirrors_support(normal):
                if value < mirror_top:
                    raise ValueError(
                        f"input plane below the {role} mirror: the plane is at n.p = {value:.3f} "
                        f"m but the mirror reaches {mirror_top:.3f} m along the input normal; "
                        "raise the input plane above it"
                    )
        if acceptance is not None:
            wavelengths, heights, values = (np.asarray(a, dtype=np.float64) for a in acceptance)
            if wavelengths.ndim != 1 or heights.ndim != 1:
                raise ValueError("acceptance axes must be one-dimensional")
            if values.shape != (wavelengths.size, heights.size):
                raise ValueError("acceptance values must have shape (n_wavelength, n_zem)")
            acceptance = (
                np.ascontiguousarray(wavelengths),
                np.ascontiguousarray(heights),
                np.ascontiguousarray(values).ravel(),
            )
        if wavelength_density is not None:
            nodes, weights = (np.ascontiguousarray(a, dtype=np.float64) for a in wavelength_density)
            if nodes.ndim != 1 or nodes.shape != weights.shape:
                raise ValueError("wavelength_density needs two equally long 1-D arrays")
            wavelength_density = (nodes, weights)
        return self._native.trace_bunches(
            arrays["x"],
            arrays["y"],
            arrays["cx"],
            arrays["cy"],
            arrays["time_ns"],
            arrays["zem"],
            arrays["size"],
            arrays["wavelength_nm"],
            arrays["bunch_id"],
            rotation.ravel(),
            offset,
            float(length_unit_m),
            float(plane_z),
            cz=arrays.get("cz"),
            z=arrays.get("z"),
            path=path,
            lower_plane_z=None if lower_plane_z is None else float(lower_plane_z),
            seed=seed,
            event=event,
            telescope=telescope,
            acceptance=acceptance,
            wavelength_density=wavelength_density,
            wavelength_law=None
            if wavelength_law is None
            else tuple(float(v) for v in wavelength_law),
            roulette=roulette,
            classify_losses=classify_losses,
            normalise_transmission=normalise_transmission,
            pool=pool if pool is not None else self._shared_pool(threads),
        )


# Margin below the lowest optic for the start of photons moving up [m] (as the input plane).
_LOWER_START_MARGIN_M = 1.0


def _positions(ids: Sequence[int] | np.ndarray, id_list: Sequence[int]) -> np.ndarray:
    """Index in `id_list` of each of `ids` (each must be in the list)."""
    keys = np.asarray(id_list, dtype=np.int64)
    wanted = np.asarray(ids, dtype=np.int64).reshape(-1)
    order = np.argsort(keys, kind="stable")
    found = np.searchsorted(keys[order], wanted)
    if wanted.size and (
        np.any(found >= keys.size)
        or np.any(keys[order][np.minimum(found, keys.size - 1)] != wanted)
    ):
        raise ModelError("obscurer id of a crossing is not an obscurer of the model")
    return order[found].astype(np.int64)


def _photon_arrays(photons: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """The photon columns as the contiguous arrays of the native tracer."""
    dtypes = {
        "position_m": np.float64,
        "direction": np.float64,
        "wavelength_nm": np.float64,
        "time_ns": np.float64,
        "weight": np.float64,
        "photon_id": np.uint64,
    }
    return {
        name: np.ascontiguousarray(photons[name], dtype=dtype) for name, dtype in dtypes.items()
    }


def trace_many(
    requests: Sequence[tuple[OpticalModel, Mapping[str, np.ndarray], Mapping[str, Any]]],
    *,
    seed: int,
    pool: ThreadPool,
) -> list[dict[str, np.ndarray]]:
    """Trace the photons of several models in one job on `pool`.

    Each request is ``(model, photons, options)``; `options` may hold
    `apply_responses`, `direct`, `outputs`, `detected_only` and
    `classify_losses` with the meaning and defaults of `OpticalModel.trace`.
    The result is the list of the dicts that one `OpticalModel.trace` call per
    request returns, with the same values; the blocks and models share the
    pool's threads, which keeps them busy when single blocks are small.
    """
    allowed = {"apply_responses", "direct", "outputs", "detected_only", "classify_losses"}
    items = []
    for model, photons, options in requests:
        unknown = set(options) - allowed
        if unknown:
            raise TypeError(f"unknown trace options: {', '.join(sorted(unknown))}")
        outputs = options.get("outputs")
        items.append((
            model._native,  # noqa: SLF001
            _photon_arrays(photons),
            options.get("apply_responses", True),
            options.get("direct", False),
            None if outputs is None else list(outputs),
            options.get("detected_only", False),
            options.get("classify_losses", True),
        ))
    return _core.trace_many(items, seed=seed, pool=pool)


def _emitter_count(source: _core.Emitter, emitted: float, seed: int, first_id: int) -> int:
    """Poisson number of candidate photons of an emitter (see ``OpticalModel.emitter_photons``).

    Drawn by the core from the Philox key of (seed, first_id), like the other photon draws.
    """
    return source.candidate_count(float(emitted), first_id, seed)
