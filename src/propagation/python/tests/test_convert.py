"""Conversion of legacy extinction tables to the dlufd transmission-table format."""

import numpy as np
import pytest
from obdeect.propagation import TransmissionTable, convert_table
from obdeect.propagation.convert import FORMATS, main, read_legacy_table

OBSERVATION_M = 2000.0


def cumulative(wavelength_nm, altitude_km):
    """Smooth, increasing synthetic cumulative optical depth."""
    scale = 1.0 + (500.0 - wavelength_nm) ** 2 / 1.0e5
    return 0.02 * scale * (altitude_km + 0.05 * altitude_km**2)


def expected(wavelength_nm, height_m):
    """tau between the observation level and height_m with c linear in wavelength and height."""
    altitudes = np.arange(0.0, 51.0)
    wavelengths = np.array([180.0 + 5.0 * i for i in range(5)])
    grid = np.array([[cumulative(w, a) for a in altitudes] for w in wavelengths])
    at_height = [np.interp(height_m / 1000.0, altitudes, row) for row in grid]
    at_obs = [np.interp(OBSERVATION_M / 1000.0, altitudes, row) for row in grid]
    return np.interp(wavelength_nm, wavelengths, np.array(at_height) - np.array(at_obs))


def write_atmabs(path, wavelengths=(180, 185, 190, 195, 200), per_line=8):
    lines = ["# EXTINCTION table (synthetic)", " EXTINCTION header 2"]
    for w in wavelengths:
        lines.append(f"  {w}")
        values = [f"{cumulative(w, a):.9e}" for a in range(51)]
        for start in range(0, 51, per_line):
            lines.append(" ".join(values[start : start + per_line]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_matrix(path, first_nm, columns=5, rows=51):
    lines = [
        " ".join(f"{cumulative(first_nm + 5.0 * j, a):.9e}" for j in range(columns))
        for a in range(rows)
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def converted(tmp_path, source, format, **kwargs):
    out = tmp_path / "out.dat"
    convert_table(source, out, format=format, observation_level_m=OBSERVATION_M, **kwargs)
    return out


def heights_of(path):
    header = next(line for line in path.read_text().splitlines() if "H2=" in line)
    return np.array([float(v) for v in header.split("H1=")[1].split()]) * 1000.0


def check_against_expected(path):
    table = TransmissionTable(str(path), OBSERVATION_M)
    assert table.first_wavelength_nm == 180
    assert table.last_wavelength_nm == 200
    heights = heights_of(path)
    # Exact at the grid points, for integer wavelengths between the source columns too.
    for wavelength in (180, 182, 185, 193, 200):
        for height in heights:
            assert table.optical_depth(wavelength, height) == pytest.approx(
                expected(wavelength, height), rel=1e-9, abs=1e-12
            )
    # Between grid points the log-height interpolation follows the linear source.
    for wavelength in (183, 197):
        for height in (2030.0, 2330.0, 3777.0, 7777.0, 21234.0):
            assert table.optical_depth(wavelength, height) == pytest.approx(
                expected(wavelength, height), rel=1e-2, abs=1e-6
            )
    assert table.optical_depth(190, OBSERVATION_M) == 0.0


@pytest.mark.parametrize("per_line", [1, 8, 51])
def test_corsika_atmabs(tmp_path, per_line):
    source = write_atmabs(tmp_path / "atmabs.dat", per_line=per_line)
    check_against_expected(converted(tmp_path, source, "corsika-atmabs"))


@pytest.mark.parametrize(
    ("format", "first_nm"), [("kascade", 180.0), ("modtran4", 200.0), ("modtran5", 205.0)]
)
def test_matrix_formats_default_minimum_wavelength(tmp_path, format, first_nm):
    source = write_matrix(tmp_path / "matrix.ext", first_nm)
    table = read_legacy_table(source, format)
    assert table.wavelength_nm[0] == first_nm
    assert table.depth.shape == (5, 51)
    assert table.depth[2, 10] == pytest.approx(cumulative(first_nm + 10.0, 10.0))


def test_kascade_conversion(tmp_path):
    source = write_matrix(tmp_path / "us76.ext", 180.0)
    check_against_expected(converted(tmp_path, source, "kascade"))


def test_min_wavelength_option(tmp_path):
    source = write_matrix(tmp_path / "us76.ext", 180.0)
    # modtran4 would start at 200 nm; the option overrides it.
    check_against_expected(converted(tmp_path, source, "modtran4", min_wavelength_nm=180.0))
    assert FORMATS["modtran4"] == 200.0 and FORMATS["modtran5"] == 205.0


def test_only_differences_matter(tmp_path):
    source = write_matrix(tmp_path / "a.ext", 180.0)
    shifted = tmp_path / "b.ext"
    rows = [[float(v) + 7.0 for v in line.split()] for line in source.read_text().splitlines()]
    shifted.write_text("\n".join(" ".join(f"{v:.9e}" for v in r) for r in rows) + "\n")
    first = TransmissionTable(str(converted(tmp_path, source, "kascade")), OBSERVATION_M)
    out = tmp_path / "shifted.dat"
    convert_table(shifted, out, format="kascade", observation_level_m=OBSERVATION_M)
    second = TransmissionTable(str(out), OBSERVATION_M)
    assert second.optical_depth(190, 9000.0) == pytest.approx(first.optical_depth(190, 9000.0))


def test_command_line(tmp_path, capsys):
    source = write_atmabs(tmp_path / "atmabs.dat")
    out = tmp_path / "cli.dat"
    status = main([
        str(source),
        "-o",
        str(out),
        "--format",
        "corsika-atmabs",
        "--observation-level",
        "2000",
    ])
    assert status == 0
    check_against_expected(out)
    status = main(
        [str(tmp_path / "missing"), "-o", str(out), "--format", "kascade"]
        + ["--observation-level", "2000"]
    )
    assert status == 1
    assert "dlufd-convert-table" in capsys.readouterr().err


def bad(tmp_path, text, format="kascade", **kwargs):
    source = tmp_path / "bad.dat"
    source.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError):
        convert_table(
            source,
            tmp_path / "out.dat",
            format=format,
            observation_level_m=kwargs.pop("observation_level_m", OBSERVATION_M),
            **kwargs,
        )


def test_malformed_files_rejected(tmp_path):
    good_rows = "\n".join("0.1 0.2 0.3" for _ in range(6)) + "\n"
    bad(tmp_path, "")  # empty
    bad(tmp_path, "0.1 0.2 0.3\n0.1 0.2\n" + good_rows)  # ragged
    bad(tmp_path, "0.1 abc 0.3\n" + good_rows)  # not a number
    bad(tmp_path, "0.1 nan 0.3\n" + good_rows)  # not finite
    bad(tmp_path, "0.1 0.2 0.3\n")  # one altitude row
    bad(tmp_path, "0.1\n0.2\n0.3\n")  # one wavelength column
    # Decreasing with altitude above the site.
    bad(tmp_path, "\n".join(f"{1.0 - 0.1 * i} {1.0 - 0.1 * i}" for i in range(6)) + "\n")
    # Observation level outside the table.
    bad(tmp_path, good_rows, observation_level_m=9000.0)
    bad(tmp_path, good_rows, observation_level_m=-5.0)
    bad(tmp_path, good_rows, observation_level_m=float("nan"))
    bad(tmp_path, good_rows, format="unknown")
    bad(tmp_path, good_rows, min_wavelength_nm=-1.0)


def test_malformed_atmabs_rejected(tmp_path):
    values = " ".join(f"{0.01 * i:.4f}" for i in range(51))
    block = f"180\n{values}\n"
    bad(tmp_path, block + "185\n" + values + " 9.9\n", "corsika-atmabs")  # 52 values
    bad(tmp_path, block + "185\n" + " ".join(values.split()[:50]) + "\n", "corsika-atmabs")
    bad(tmp_path, block + "190\n" + values + "\n", "corsika-atmabs")  # 10 nm step
    bad(tmp_path, block + "185 186\n" + values + "\n", "corsika-atmabs")  # two numbers
    bad(tmp_path, "180\n" + values.replace("0.05", "x") + "\n", "corsika-atmabs")
    bad(tmp_path, block, "corsika-atmabs", min_wavelength_nm=180.0)  # one wavelength only
    bad(tmp_path, "# EXTINCTION only\n", "corsika-atmabs")


def test_local_decrease_above_the_site_rejected(tmp_path):
    # c rises from the site (2 km) to 3 km, then falls at 4 km although it stays above c(2 km).
    rows = ["0 0"] * 2 + ["1 1", "3 3", "2 2"] + ["4 4"] * 3
    bad(tmp_path, "\n".join(rows) + "\n")


def test_step_must_be_exactly_5_nm(tmp_path):
    values = " ".join(f"{0.01 * i:.4f}" for i in range(51))
    bad(tmp_path, f"180\n{values}\n185.00004\n{values}\n", "corsika-atmabs")


def test_too_many_heights_rejected_or_coarsened(tmp_path):
    # 200 km of altitude rows cannot fit 100 heights above a 2 km site.
    rows = "\n".join(f"{0.01 * a} {0.02 * a}" for a in range(201)) + "\n"
    bad(tmp_path, rows)
    # 60 rows still convert, and the output is readable.
    rows = "\n".join(f"{0.01 * a} {0.02 * a}" for a in range(61)) + "\n"
    source = tmp_path / "tall.ext"
    source.write_text(rows, encoding="utf-8")
    out = tmp_path / "tall.dat"
    convert_table(source, out, format="kascade", observation_level_m=OBSERVATION_M)
    assert len(heights_of(out)) < 100  # H2 takes one of the 100 slots
    TransmissionTable(str(out), OBSERVATION_M)


def test_small_local_decrease_in_a_large_table_rejected(tmp_path):
    # A drop of 1e-6 relative beside values of 1e12 is a decrease, not rounding; the tolerance
    # follows the neighbouring values, not the largest value of the table.
    rows = ["0 0"] * 2 + ["1e12 1e12", "2e12 2e12", "1.999999e12 1.999999e12", "3e12 3e12"]
    bad(tmp_path, "\n".join(rows) + "\n" + "\n".join(["4e12 4e12"] * 3) + "\n")
