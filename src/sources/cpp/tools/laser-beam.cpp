// dlamb-laser-beam: primary photons of a laser beam crossing the field of view of a telescope,
// as sim_telarray LightEmission/ls-beam.cc (primary light only; scattering is dlufd's job).
// Writes a gstrickd photon stream (sky frame) to stdout or -o FILE.

#include "dlamb/programs/laser-beam.hpp"

#include "sky_args.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

namespace {

const char *usage_text =
    "Usage: dlamb-laser-beam [options]\n"
    "Primary photons of a laser beam inside the telescope field of view; stream on stdout.\n"
    "Options (original ls-beam.cc option, unit, default):\n"
    "  --events n            events (1)\n"
    "  --step m              step along the beam [m] (0.1)\n"
    "  --bunchsize f         photons per bunch (3; values below 1 are used as 1)\n"
    "  --primary-per-m f     primary photons per m of beam (1e4; NEW, replaces --brightness)\n"
    "  --spectrum wlen|fname wavelength [nm] or table (355)\n"
    "  --lightpulse d|fname  Simple:w, Gauss:s, Gauss-fwhm:f, Exp:tau, Complex:..., or table,\n"
    "                        a bare number is a top-hat width [ns] (Simple:0)\n"
    "  --angular-distribution name  recorded for dlufd (Rayleigh)\n"
    "  --run n               run number (1)\n"
    "  --altitude m          observation level a.s.l. [m] (2150)\n"
    "  --atmosphere n        profile number, recorded only (1)\n"
    "  --laser-position x,y,z  laser w.r.t. the telescope [m] (1000,0,0)\n"
    "  --laser-theta deg     zenith angle of the beam (60)\n"
    "  --laser-phi deg       azimuth of the beam, astronomical sense (180)\n"
    "  --telescope-theta deg zenith angle of the telescope axis (30)\n"
    "  --telescope-phi deg   azimuth of the telescope axis, astronomical sense (0)\n"
    "  --telescope-radius m  fiducial radius [m] (12.5)\n"
    "  --fov deg             field of view, full cone opening angle, 0 < fov < 90 (15)\n"
    "  --seed n              random seed (1, new option)\n"
    "  -o fname              output file (stdout)\n"
    "  --help                this text\n";

// Astronomical azimuth to the mathematical sense, as ls-beam.cc.
double math_azimuth(double astronomical_deg) {
  return astronomical_deg <= 0.0 ? -astronomical_deg : 360.0 - astronomical_deg;
}

int run(int argc, char **argv) {
  dlamb::programs::LsBeamOptions o;
  std::string output;
  std::string spectrum_text = "355";
  std::string pulse_text = "Simple:0";
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    const auto next = [&]() -> std::string {
      if (i + 1 >= argc)
        throw std::invalid_argument(arg + ": missing value");
      return argv[++i];
    };
    if (arg == "--help") {
      std::fputs(usage_text, stdout);
      return 0;
    } else if (arg == "--brightness") {
      throw std::invalid_argument("--brightness is replaced by --primary-per-m "
                                  "(see docs/programs/laser-beam.md)");
    } else if (arg == "--events") {
      o.events = dlamb::cli::to_count(next(), arg);
    } else if (arg == "--step" || arg == "--steplength") {
      o.step_m = dlamb::cli::to_double(next(), arg);
    } else if (arg == "--bunchsize") {
      o.bunch_size = dlamb::cli::to_double(next(), arg);
    } else if (arg == "--primary-per-m") {
      o.primary_photons_per_m = dlamb::cli::to_double(next(), arg);
    } else if (arg == "--spectrum") {
      spectrum_text = next();
    } else if (arg == "--lightpulse" || arg == "--pulse") {
      pulse_text = next();
    } else if (arg == "--angular" || arg == "--angular-distribution") {
      o.angular = next();
    } else if (arg == "--run") {
      o.run_number = dlamb::cli::to_int(next(), arg);
    } else if (arg == "--altitude") {
      o.altitude_m = dlamb::cli::to_double(next(), arg);
    } else if (arg == "--atmosphere") {
      o.atmosphere = dlamb::cli::to_int(next(), arg);
    } else if (arg == "--laser-position") {
      const auto v = dlamb::cli::to_list(next(), arg);
      if (v.size() < 2 || v.size() > 3)
        throw std::invalid_argument("--laser-position needs x,y or x,y,z");
      o.laser_position_m = {v[0], v[1], v.size() == 3 ? v[2] : 0.0};
    } else if (arg == "--laser-theta") {
      o.laser_theta_deg = dlamb::cli::to_double(next(), arg);
    } else if (arg == "--laser-phi") {
      o.laser_phi_deg = math_azimuth(dlamb::cli::to_double(next(), arg));
    } else if (arg == "--telescope-theta") {
      o.telescope_theta_deg = dlamb::cli::to_double(next(), arg);
    } else if (arg == "--telescope-phi") {
      o.telescope_phi_deg = math_azimuth(dlamb::cli::to_double(next(), arg));
    } else if (arg == "--telescope-radius") {
      o.telescope_radius_m = dlamb::cli::to_double(next(), arg);
    } else if (arg == "--fov") {
      o.fov_deg = dlamb::cli::to_double(next(), arg);
    } else if (arg == "--seed") {
      o.seed = static_cast<std::uint64_t>(dlamb::cli::to_int(next(), arg));
    } else if (arg == "-o") {
      output = next();
    } else {
      throw std::invalid_argument("unknown option " + arg + " (see --help)");
    }
  }
  o.spectrum = dlamb::cli::parse_spectrum(spectrum_text, "--spectrum");
  o.pulse = dlamb::cli::parse_pulse(pulse_text, "--lightpulse", false);

  if (output.empty()) {
    std::ios::sync_with_stdio(false);
    dlamb::programs::write_laser_beam(std::cout, o);
    std::cout.flush();
  } else {
    std::ofstream file(output, std::ios::binary);
    if (!file)
      throw std::invalid_argument("cannot write " + output);
    dlamb::programs::write_laser_beam(file, o);
  }
  return 0;
}

} // namespace

int main(int argc, char **argv) {
  try {
    return run(argc, argv);
  } catch (const std::exception &error) {
    std::fprintf(stderr, "dlamb-laser-beam: %s\n", error.what());
    return 1;
  }
}
