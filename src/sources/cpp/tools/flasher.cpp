// dlamb-flasher: flasher photon stream, dish or secondary type.

#include "dlamb/programs/flasher.hpp"

#include "tool_args.hpp"

#include <iostream>
#include <utility>

namespace {

const char *usage = R"(usage: dlamb-flasher [--type dish|secondary] [options] [-o FILE]
  both:      --photons N[,..]  --bunchsize F (--bs, -b)  --spectrum W|Cherenkov:A:B
             --lightpulse S (--pulse)  --angular-distribution S (--angular)  --run N
             --altitude M (-a)  --atmosphere N (-A)  --seed N  --events N[,..]  -o FILE
  dish:      --frame sky|camera  --xy X,Y  --distance D (-z)  --camera-radius R (-r)
  secondary: --flasher-xy X (--xy)  --flasher-depth D (-z, --z)  --flasher-inclination A
             --mirror-camera-distance D  --delays T1,..,T4[,T5]  --fire all|alternate|demo|1..4
)";

} // namespace

int main(int argc, char **argv) {
  using namespace dlamb;
  using programs::flasher::Type;
  const std::string prefix = argc > 0 ? argv[0] : "dlamb-flasher";
  return cli::run_tool(prefix, [&] {
    std::vector<std::pair<std::string, std::string>> pairs;
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "-h" || arg == "--help") {
        std::cout << usage;
        return 0;
      }
      if (i + 1 >= argc)
        throw std::invalid_argument("option " + arg + " needs a value");
      pairs.emplace_back(arg, argv[++i]);
    }
    programs::flasher::Options o;
    for (const auto &[name, text] : pairs)
      if (name == "--type")
        o.type = programs::flasher::parse_type(text);
    const bool dish = o.type == Type::dish;
    const auto only = [&](const std::string &name, Type type) {
      programs::flasher::require_type(o.type, type, name);
    };
    std::string output;
    for (const auto &[name, text] : pairs) {
      if (name == "--type") {
      } else if (name == "--photons") {
        o.photons = programs::parse_real_list(text, name);
      } else if (name == "--bunchsize" || name == "--bs" || name == "-b") {
        o.bunch_size = programs::parse_real(text, name);
      } else if (name == "--spectrum") {
        o.spectrum = text;
      } else if (name == "--lightpulse" || name == "--pulse") {
        o.lightpulse = text;
      } else if (name == "--angular-distribution" || name == "--angular") {
        o.angular = text;
      } else if (name == "--run") {
        o.run_number = programs::parse_integer_value(text, name);
      } else if (name == "--altitude" || name == "-a") {
        o.altitude_m = programs::parse_real(text, name);
      } else if (name == "--atmosphere" || name == "-A") {
        o.atmosphere = programs::parse_integer_value(text, name);
      } else if (name == "--seed") {
        o.seed = cli::parse_seed(text, name);
      } else if (name == "-o") {
        output = text;
      } else if (name == "--events") {
        o.events.clear();
        for (const auto &item : cli::split_items(text, name))
          o.events.push_back(static_cast<int>(programs::parse_integer_value(item, name)));
      } else if (name == "--frame") {
        only(name, Type::dish);
        if (text == "sky")
          o.frame = gstrickd::Frame::sky;
        else if (text == "camera")
          o.frame = gstrickd::Frame::camera;
        else
          throw std::invalid_argument("--frame must be sky or camera");
      } else if (name == "--xy" || name == "--flasher-xy") {
        const auto values = programs::parse_real_list(text, name);
        if (dish) {
          if (values.size() != 2)
            throw std::invalid_argument(name + " needs x,y with --type dish");
          o.x_cm = values[0];
          o.y_cm = values[1];
        } else {
          if (values.size() != 1)
            throw std::invalid_argument(name + " needs one value with --type secondary");
          o.flasher_xy_cm = values[0];
        }
      } else if (name == "--distance") {
        only(name, Type::dish);
        o.distance_cm = programs::parse_real(text, name);
      } else if (name == "-z") {
        (dish ? o.distance_cm : o.flasher_depth_cm) = programs::parse_real(text, name);
      } else if (name == "--camera-radius" || name == "-r") {
        only(name, Type::dish);
        o.camera_radius_cm = programs::parse_real(text, name);
      } else if (name == "--flasher-depth" || name == "--z") {
        only(name, Type::secondary);
        o.flasher_depth_cm = programs::parse_real(text, name);
      } else if (name == "--flasher-inclination" || name == "--inclination") {
        only(name, Type::secondary);
        o.inclination_deg = programs::parse_real(text, name);
      } else if (name == "--mirror-camera-distance" || name == "--mirror-distance") {
        only(name, Type::secondary);
        o.mirror_camera_distance_cm = programs::parse_real(text, name);
      } else if (name == "--delays") {
        only(name, Type::secondary);
        const auto values = programs::parse_real_list(text, name);
        if (values.size() != 4 && values.size() != 5)
          throw std::invalid_argument("--delays needs 4 or 5 values");
        o.delays = {values[0], values[1], values[2], values[3],
                    values.size() == 5 ? values[4] : 0.0};
      } else if (name == "--fire") {
        only(name, Type::secondary);
        o.fire = text;
      } else {
        throw std::invalid_argument("unknown option " + name);
      }
    }
    cli::Output out(output);
    programs::flasher::write_stream(out.stream(), o);
    out.close();
    return 0;
  });
}
