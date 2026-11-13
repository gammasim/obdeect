// dlamb-muon-track: Cherenkov light of one muon-like track on a telescope sphere (sky frame).

#include "obdeect/sources/programs/muon-track.hpp"

#include "tool_args.hpp"

#include <cstdio>
#include <iostream>

int main(int argc, char **argv) {
  using namespace obdeect::sources;
  std::vector<std::string> args;
  for (int i = 1; i < argc; ++i)
    args.emplace_back(argv[i]);
  try {
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
        stderr, "dlamb-muon-track: expected weight %.6g, wrote %llu photons, weight %.6g\n",
        summary.expected, static_cast<unsigned long long>(summary.photons), summary.weight);
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "dlamb-muon-track: " << error.what() << '\n' << muon_track::usage();
    return 1;
  }
}
