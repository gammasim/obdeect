#pragma once

#include "obdeect/stream/random.hpp"
#include "obdeect/stream/stages.hpp"

#include <array>
#include <cstdint>
#include <stdexcept>

namespace obdeect::ray_tracing {

// Philox-4x32-10 (Salmon et al., SC'11), shared with the other components: a
// counter-based generator. A draw depends only on its key and counter, never on
// the order in which photons are processed, so results are identical for any
// block size and thread count.
using obdeect::stream::philox4x32;

// Independent random streams of one photon.
// New values are appended; existing ones are never renumbered.
enum class Stream : std::uint32_t {
  source = 0,
  primary_scatter = 1,
  secondary_scatter = 2,
  thinning = 3,
  wavelength = 4
};

// Key of the ray-tracing draws, shared with the other components
// (obdeect::stream::DrawKey with stage_ray_tracing). Outside photon streams the
// event and telescope are 0, so a run seed alone is a key.
struct RandomKey {
  std::uint64_t seed{};
  std::uint32_t event{};
  std::uint32_t telescope{};
  constexpr RandomKey(std::uint64_t run_seed = 0, std::uint32_t event_number = 0,
                      std::uint32_t telescope_id = 0)
      : seed(run_seed), event(event_number), telescope(telescope_id) {}
};

// Photon ids [first, first + count) fit the key of the shared generator (40 bits). Checked
// before tracing: a draw outside the key would throw inside a worker thread.
inline void check_photon_ids(std::uint64_t first, std::uint64_t count) {
  if (count > 0 && (first > obdeect::stream::max_photon_index ||
                    count - 1 > obdeect::stream::max_photon_index - first))
    throw std::out_of_range("photon ids must lie below 2^40 (key of the random numbers)");
}

// Four uniform deviates in the open interval (0, 1).
[[nodiscard]] constexpr std::array<double, 4>
uniform4(const RandomKey &key, std::uint64_t photon_id, Stream stream, std::uint32_t block = 0) {
  return obdeect::stream::uniform4({key.seed, obdeect::stream::stage_ray_tracing, key.event,
                                    key.telescope, photon_id, static_cast<std::uint32_t>(stream)},
                                   block);
}

} // namespace obdeect::ray_tracing
