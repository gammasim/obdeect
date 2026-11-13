// obdeect sources beam-scan: parallel downward beam scanned over the camera area (sky frame).

#include "obdeect/sources/programs/beam-scan.hpp"

#include "tool_args.hpp"

#include <iostream>

namespace {

const char *usage = R"(usage: obdeect sources beam-scan [options] [-o FILE]
  -r rcam    camera radius [cm] (123)
  -x xpos    scan start x [cm] (0.3 * 123)    -y ypos  scan start y [cm] (0.2 * 123)
  -h height  source height [cm] (0.05 * 123)
  --seed n   -o FILE   --help
)";

} // namespace

int main(int argc, char **argv) {
  using namespace obdeect::sources;
  return cli::run_tool("obdeect sources beam-scan", argc, argv, [&] {
    cli::Args cursor(argc, argv);
    programs::XylsOptions o;
    std::string output, arg;
    while (cursor.next(arg)) {
      if (arg == "--help") {
        std::cout << usage;
        return 0;
      }
      const auto value = [&] { return cursor.value(arg); };
      if (arg == "-r")
        o.radius_cm = cli::parse_number(value(), arg);
      else if (arg == "-x")
        o.x_cm = cli::parse_number(value(), arg);
      else if (arg == "-y")
        o.y_cm = cli::parse_number(value(), arg);
      else if (arg == "-h")
        o.height_cm = cli::parse_number(value(), arg);
      else if (arg == "--seed")
        o.seed = cli::parse_seed(value(), arg);
      else if (arg == "-o")
        output = value();
      else
        throw std::invalid_argument("unknown option " + arg + " (see --help)");
    }
    cli::Output out(output);
    programs::write_beam_scan(out.stream(), o);
    out.close();
    return 0;
  });
}
