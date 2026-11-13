// obdeect sources muon-track: Cherenkov light of one muon-like track on a telescope sphere (sky
// frame).

#include "obdeect/sources/programs/muon-track.hpp"

#include "tool_args.hpp"

#include <cstdio>
#include <iostream>

int main(int argc, char **argv) {
  using namespace obdeect::sources;
  return cli::run_tool(
      "obdeect sources muon-track", argc, argv,
      [&] {
        cli::Args cursor(argc, argv);
        std::vector<std::string> args;
        std::string arg;
        while (cursor.next(arg))
          args.push_back(arg);
        const auto parsed = muon_track::parse_command_line(args);
        if (parsed.help) {
          std::cout << muon_track::usage();
          return 0;
        }
        muon_track::validate(parsed.options);
        cli::Output out(parsed.output);
        const StreamSummary summary = muon_track::write(out.stream(), parsed.options);
        out.close();
        std::fprintf(
            stderr,
            "obdeect sources muon-track: expected weight %.6g, wrote %llu photons, weight %.6g\n",
            summary.expected, static_cast<unsigned long long>(summary.photons), summary.weight);
        return 0;
      },
      [] { std::cerr << muon_track::usage(); });
}
