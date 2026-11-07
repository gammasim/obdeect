// Writes a synthetic CORSIKA IACT EventIO file (normal mode, long-format photon blocks) for
// throughput tests of corsika2ps and ps2simtel.
// Usage: gen_corsika OUTPUT [BUNCHES=10000000] [TELESCOPES=4] [BUNCHES_PER_BLOCK=250000]

#include "gstrickd/eventio.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <string>
#include <vector>

namespace {

std::vector<float> words273() { return std::vector<float>(273, 0.0F); }

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: gen_corsika OUTPUT [BUNCHES] [TELESCOPES] [BUNCHES_PER_BLOCK]\n");
    return 2;
  }
  const auto total = argc > 2 ? std::strtoull(argv[2], nullptr, 10) : 10000000ULL;
  const auto telescopes = argc > 3 ? std::strtoull(argv[3], nullptr, 10) : 4ULL;
  const auto per_block = argc > 4 ? std::strtoull(argv[4], nullptr, 10) : 250000ULL;
  if (total == 0 || telescopes == 0 || per_block == 0) {
    std::fprintf(stderr, "gen_corsika: counts must be positive\n");
    return 2;
  }
  const std::uint64_t per_event = telescopes * per_block;
  const std::uint64_t events = (total + per_event - 1) / per_event;

  std::ofstream out(argv[1], std::ios::binary);
  if (!out) {
    std::fprintf(stderr, "gen_corsika: cannot open %s\n", argv[1]);
    return 1;
  }
  std::vector<char> buffer(std::size_t{1} << 20);
  out.rdbuf()->pubsetbuf(buffer.data(), static_cast<std::streamsize>(buffer.size()));

  std::vector<float> runh = words273();
  runh[1] = 1.0F;
  runh[3] = 7.41F;
  runh[4] = 1.0F;
  runh[5] = 2200000.0F;
  runh[16] = 10.0F;
  runh[17] = 100000.0F;
  gstrickd::write_block(out, gstrickd::make_float_block(gstrickd::iact::runh, 1, runh));

  gstrickd::TelPos pos;
  for (std::uint64_t t = 0; t < telescopes; ++t) {
    pos.x.push_back(10000.0 * static_cast<double>(t));
    pos.y.push_back(0.0);
    pos.z.push_back(0.0);
    pos.r.push_back(1200.0);
  }
  gstrickd::write_block(out, gstrickd::make_telpos_block(pos));

  std::mt19937 rng{20261107};
  std::uniform_real_distribution<float> position{-10000.0F, 10000.0F};
  std::uniform_real_distribution<float> direction{-0.3F, 0.3F};
  std::uniform_real_distribution<float> time{0.0F, 100.0F};
  for (std::uint64_t e = 1; e <= events; ++e) {
    std::vector<float> evth = words273();
    evth[1] = static_cast<float>(e);
    evth[2] = 14.0F;
    evth[3] = 1000.0F;
    evth[6] = 2500000.0F;
    evth[10] = 0.3F;
    evth[11] = 1.2F;
    evth[76] = 5136.0F;
    evth[95] = 250.0F;
    evth[96] = 800.0F;
    gstrickd::write_block(
        out, gstrickd::make_float_block(gstrickd::iact::evth, static_cast<std::uint32_t>(e), evth));
    gstrickd::TelOffset off;
    off.xoff = {0.0};
    off.yoff = {0.0};
    gstrickd::write_block(out, gstrickd::make_teloff_block(off));

    std::vector<gstrickd::EventioBlock> children;
    for (std::uint64_t t = 0; t < telescopes; ++t) {
      gstrickd::TelPhotons p;
      p.array = 0;
      p.tel = static_cast<int>(t);
      p.bunches.reserve(per_block);
      for (std::uint64_t i = 0; i < per_block; ++i) {
        gstrickd::Bunch b;
        b.x = position(rng);
        b.y = position(rng);
        b.cx = direction(rng);
        b.cy = direction(rng);
        b.ctime = time(rng);
        b.zem = 2000000.0F;
        b.photons = (i % 3 == 0) ? 2.0F : 1.0F;
        b.lambda = 400.0F;
        p.photons += b.photons;
        p.bunches.push_back(b);
      }
      children.push_back(gstrickd::make_photons_block(p));
    }
    gstrickd::write_block(out, gstrickd::make_container(gstrickd::iact::telarray, 0, 0, children));
    std::vector<float> evte = words273();
    evte[1] = static_cast<float>(e);
    gstrickd::write_block(
        out, gstrickd::make_float_block(gstrickd::iact::evte, static_cast<std::uint32_t>(e), evte));
  }
  std::vector<float> rune = words273();
  rune[1] = 1.0F;
  gstrickd::write_block(out, gstrickd::make_float_block(gstrickd::iact::rune, 1, rune));
  out.close();
  if (!out) {
    std::fprintf(stderr, "gen_corsika: write failed\n");
    return 1;
  }
  std::printf("%llu events, %llu telescopes, %llu bunches\n",
              static_cast<unsigned long long>(events), static_cast<unsigned long long>(telescopes),
              static_cast<unsigned long long>(events * per_event));
  return 0;
}
