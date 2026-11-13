// obdeect sources led: LED in front of selected camera pixels (camera frame).

#include "obdeect/sources/programs/led.hpp"

#include "tool_args.hpp"

#include <fstream>
#include <iostream>
#include <utility>

namespace {

const char *usage =
    R"(usage: obdeect sources led --camera-file FILE (--use LIST | --use-trg-rand N | --use-trg-all) [options]
  --camera-file, -c FILE   camera definition (Pixel and Trigger lines)
  --use, -u LIST           pixels: ids, + required, [a,b] slaves, % comment
  --use-trg-rand N         N random trigger groups      --use-trg-all  all trigger groups
  --require, --required N  pixels per event (0 = all)   --random N     N random pixels per event
  --slaves                 also fire pre-sum slaves     --events N      events per combination
  --photons N  --bunchsize, -b F  --camera-radius, -r R [cm]
  --spectrum W  --lightpulse, --pulse S  --run, -R N  --altitude, -a M  --atmosphere, -A N
  --seed N  --list-only, -L  --quiet, -q  -o FILE  -h, --help
)";

} // namespace

int main(int argc, char **argv) {
  using namespace obdeect::sources;
  return cli::run_tool("obdeect sources led", argc, argv, [&] {
    programs::led::Options o;
    std::string camera_file, output;
    bool list_only = false;
    for (int i = 1; i < argc; ++i) {
      const std::string name = argv[i];
      if (name == "-h" || name == "--help") {
        std::cout << usage;
        return 0;
      }
      if (name == "--slaves") {
        o.slaves = true;
        continue;
      }
      if (name == "--use-trg-all") {
        o.trigger_all = true;
        continue;
      }
      if (name == "--list-only" || name == "-L") {
        list_only = true;
        continue;
      }
      if (name == "--quiet" || name == "-q")
        continue;
      static const char *const with_value[] = {"--camera-file",
                                               "-c",
                                               "--use",
                                               "-u",
                                               "--use-trg-rand",
                                               "--require",
                                               "--required",
                                               "--random",
                                               "--events",
                                               "--photons",
                                               "--bunchsize",
                                               "-b",
                                               "--camera-radius",
                                               "-r",
                                               "--spectrum",
                                               "--lightpulse",
                                               "--pulse",
                                               "--run",
                                               "-R",
                                               "--altitude",
                                               "-a",
                                               "--atmosphere",
                                               "-A",
                                               "--seed",
                                               "-o"};
      bool known = false;
      for (const char *w : with_value)
        known = known || name == w;
      if (!known)
        throw std::invalid_argument("unknown option " + name);
      if (i + 1 >= argc)
        throw std::invalid_argument("option " + name + " needs a value");
      const std::string text = argv[++i];
      const auto integer = [&] {
        return static_cast<int>(programs::parse_integer_value(text, name));
      };
      if (name == "--camera-file" || name == "-c")
        camera_file = text;
      else if (name == "--use" || name == "-u")
        o.use = text;
      else if (name == "--use-trg-rand")
        o.trigger_random = integer();
      else if (name == "--require" || name == "--required")
        o.require = integer();
      else if (name == "--random") {
        o.random = true;
        o.require = integer();
      } else if (name == "--events")
        o.events = integer();
      else if (name == "--photons")
        o.photons = programs::parse_real(text, name);
      else if (name == "--bunchsize" || name == "-b")
        o.bunch_size = programs::parse_real(text, name);
      else if (name == "--camera-radius" || name == "-r")
        o.camera_radius_cm = programs::parse_real(text, name);
      else if (name == "--spectrum")
        o.spectrum = text;
      else if (name == "--lightpulse" || name == "--pulse")
        o.lightpulse = text;
      else if (name == "--run" || name == "-R")
        o.run_number = programs::parse_integer_value(text, name);
      else if (name == "--altitude" || name == "-a")
        o.altitude_m = programs::parse_real(text, name);
      else if (name == "--atmosphere" || name == "-A")
        o.atmosphere = programs::parse_integer_value(text, name);
      else if (name == "--seed")
        o.seed = cli::parse_seed(text, name);
      else if (name == "-o")
        output = text;
    }
    if (camera_file.empty())
      throw std::invalid_argument("no camera definition file defined (--camera-file)");
    std::ifstream file(camera_file);
    if (!file)
      throw std::runtime_error("cannot open camera configuration file " + camera_file);
    o.camera = programs::led::read_camera(file, o.slaves);
    if (list_only) {
      std::cout << programs::led::label_lines(programs::led::build_sources(o));
      return 0;
    }
    cli::Output out(output);
    programs::led::write_stream(out.stream(), o);
    out.close();
    return 0;
  });
}
