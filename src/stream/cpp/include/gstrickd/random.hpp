#pragma once

#include <array>
#include <cstdint>

namespace gstrickd {

// Philox-4x32-10 (Salmon et al., SC'11), as in obdeect: a counter-based
// generator. A draw depends only on its key and counter, never on the order in
// which photons are processed, so results are identical for any chunk size and
// thread count.
[[nodiscard]] constexpr std::array<std::uint32_t, 4> philox4x32(std::array<std::uint32_t, 4> c,
                                                                std::array<std::uint32_t, 2> k) {
  constexpr std::uint32_t m0 = 0xD2511F53U, m1 = 0xCD9E8D57U;
  constexpr std::uint32_t w0 = 0x9E3779B9U, w1 = 0xBB67AE85U;
  for (int round = 0; round < 10; ++round) {
    const std::uint64_t p0 = std::uint64_t{m0} * c[0];
    const std::uint64_t p1 = std::uint64_t{m1} * c[2];
    c = {static_cast<std::uint32_t>(p1 >> 32) ^ c[1] ^ k[0], static_cast<std::uint32_t>(p1),
         static_cast<std::uint32_t>(p0 >> 32) ^ c[3] ^ k[1], static_cast<std::uint32_t>(p0)};
    k = {k[0] + w0, k[1] + w1};
  }
  return c;
}

// Identifies one random draw: stage (module), event, telescope, photon index
// within (event, telescope), and the stream of draws for that photon. Each
// module uses its own stage id; streams are numbered within the module.
struct DrawKey {
  std::uint64_t seed{};
  std::uint32_t stage{};
  std::uint32_t event{};
  std::uint32_t telescope{};
  std::uint64_t photon{};
  std::uint32_t stream{};
};

// Four uniform deviates in the open interval (0, 1). The 128-bit counter is
// (photon, telescope, event ^ stream << 24); the key mixes seed and stage.
[[nodiscard]] constexpr std::array<double, 4> uniform4(const DrawKey &key,
                                                       std::uint32_t block = 0) {
  const std::uint64_t seed = key.seed ^ (std::uint64_t{key.stage} * 0x9E3779B97F4A7C15ULL);
  const auto r = philox4x32(
      {static_cast<std::uint32_t>(key.photon), static_cast<std::uint32_t>(key.photon >> 32),
       key.telescope, key.event ^ (key.stream << 24) ^ (block << 16)},
      {static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32)});
  constexpr double scale = 1.0 / 4294967296.0;
  return {(r[0] + 0.5) * scale, (r[1] + 0.5) * scale, (r[2] + 0.5) * scale, (r[3] + 0.5) * scale};
}

} // namespace gstrickd
