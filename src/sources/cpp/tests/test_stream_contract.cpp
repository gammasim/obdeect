#include "check.hpp"

#include "obdeect/sources/programs/stream_output.hpp"
#include "obdeect/sources/source.hpp"
#include "obdeect/sources/streams.hpp"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/stages.hpp"

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <string>

// obdeect sources uses obdeect stream's stage ids, frame codes, column names and chunk limit
// directly. These checks fail when obdeect stream changes a value that obdeect sources depends on.
int main() {
  return run_checks([] {
    // Format version written and read by obdeect sources (one named constant, see streams.hpp).
    CHECK(obdeect::stream::format_major == obdeect::sources::photon_stream_format_major);

    // Stage ids are part of the random-number results: never changed.
    CHECK(obdeect::sources::stage_sources == 0x444C0001U);
    CHECK(obdeect::sources::programs::stage_programs == 0x444C0002U);
    CHECK(obdeect::stream::stage_sources == obdeect::sources::stage_sources);
    CHECK(obdeect::stream::stage_programs == obdeect::sources::programs::stage_programs);

    // Frame codes accepted by the Python binding (sky 0, T 2, C 3).
    CHECK(static_cast<std::uint32_t>(obdeect::stream::Frame::sky) == 0U);
    CHECK(static_cast<std::uint32_t>(obdeect::stream::Frame::telescope) == 2U);
    CHECK(static_cast<std::uint32_t>(obdeect::stream::Frame::camera) == 3U);

    // Columns obdeect sources emits, in obdeect stream's names and order.
    const std::string expected[] = {"x", "y", "z", "cx", "cy", "cz", "t", "w", "wavelength"};
    std::size_t k = 0;
    for (std::uint32_t c = 0; c < obdeect::stream::n_columns; ++c) {
      if ((obdeect::sources::source_columns >> c) & 1U) {
        CHECK(k < std::size(expected));
        if (k < std::size(expected))
          CHECK(std::string(obdeect::stream::column_name(
                    static_cast<obdeect::stream::Column>(c))) == expected[k]);
        ++k;
      }
    }
    CHECK(k == std::size(expected));

    CHECK(obdeect::stream::max_chunk_rows == (1U << 16));
  });
}
