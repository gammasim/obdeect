#include "check.hpp"
#include "program_stream.hpp"

#include "obdeect/sources/programs/led.hpp"

#include <cmath>
#include <cstdint>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace obdeect::sources;
using namespace obdeect::sources::programs;

namespace {

// Four pixels at (0,0), (1,0), (0,1), (1,1) cm and two trigger groups.
const char *camera_text = "# test camera\n"
                          "Pixel 0 0 0.0 0.0 1.0\n"
                          "Pixel 1 0 1.0 0.0 1.0\n"
                          "Pixel 2 0 0.0 1.0 1.0\n"
                          "Pixel 3 0 1.0 1.0 1.0\n"
                          "Trigger  * of 0 1 [2]\n"
                          "Trigger  * of 2 3\n";

led::Options base_options(const std::string &use) {
  led::Options options;
  std::istringstream in(camera_text);
  options.camera = led::read_camera(in, options.slaves);
  options.use = use;
  options.photons = 10.0;
  return options;
}

std::string stream_bytes(const led::Options &options) {
  std::ostringstream out;
  led::write_stream(out, options);
  return out.str();
}

// Pixel index of a camera-frame position (metres) of the test camera.
int pixel_at(float x, float y) {
  return static_cast<int>(std::lround(100.0 * x)) + 2 * static_cast<int>(std::lround(100.0 * y));
}

void combinations() {
  // "0 1 2" with two pixels per event: three events, each with 2 x 10 photons.
  auto options = base_options("0 1 2");
  options.require = 2;
  const auto s = read_stream(stream_bytes(options));
  CHECK(count_type(s, obdeect::stream::RecordType::event) == 3);
  CHECK(s.rows.size() == 3 * 20);
  std::set<std::set<int>> sets;
  for (std::size_t e = 0; e < 3; ++e) {
    std::set<int> lit;
    for (std::size_t i = 0; i < 20; ++i) {
      const Row &r = s.rows[e * 20 + i];
      CHECK(r.frame == 3); // camera frame
      CHECK(r.z == 0.0F);
      CHECK(r.cx == 0.0F && r.cy == 0.0F && r.cz == -1.0F);
      CHECK(r.t == 0.0F);
      CHECK(r.w == 1.0F);
      CHECK(r.event == static_cast<std::int64_t>(e + 1));
      lit.insert(pixel_at(r.x, r.y));
    }
    CHECK(lit.size() == 2);
    sets.insert(lit);
  }
  // Combinations in increasing order: {0,1}, {0,2}, {1,2}.
  CHECK(sets.size() == 3);
  CHECK(sets.count({0, 1}) == 1);
  CHECK(sets.count({0, 2}) == 1);
  CHECK(sets.count({1, 2}) == 1);
}

void slaves_and_trigger_groups() {
  // Trigger group 0 is "0 1 [2]": without --slaves the slave is dropped, with it pixel 2
  // fires too (one more emission, same labels).
  auto plain = base_options("");
  plain.use.clear();
  plain.trigger_all = true;
  auto with_slaves = plain;
  with_slaves.slaves = true;
  with_slaves.camera.triggers.clear();
  {
    std::istringstream in(camera_text);
    with_slaves.camera = led::read_camera(in, true);
  }
  const auto a = read_stream(stream_bytes(plain));
  // Groups: {0,1} and {2,3} give two events of 2 x 10 photons.
  CHECK(count_type(a, obdeect::stream::RecordType::event) == 2);
  CHECK(a.rows.size() == 40);

  const auto b = read_stream(stream_bytes(with_slaves));
  // Group 0 with slave 2 is three emissions of 10 photons; group 1 has two.
  CHECK(b.rows.size() == 30 + 20);

  auto labels = led::build_sources(with_slaves);
  CHECK(labels.size() == 2);
  CHECK(labels[0].labels.size() == 2); // slaves are not labels
}

void random_pixels() {
  // Two distinct random pixels of four per event, five events.
  auto options = base_options("0 1 2 3");
  options.require = 2;
  options.random = true;
  options.events = 5;
  const auto s = read_stream(stream_bytes(options));
  CHECK(count_type(s, obdeect::stream::RecordType::event) == 5);
  CHECK(s.rows.size() == 5 * 20);
  for (std::size_t e = 0; e < 5; ++e) {
    std::set<int> lit;
    for (std::size_t i = 0; i < 20; ++i) {
      lit.insert(pixel_at(s.rows[e * 20 + i].x, s.rows[e * 20 + i].y));
    }
    CHECK(lit.size() == 2);
  }
}

void all_pixels_and_list_errors() {
  // All listed pixels fire in one event; a required count above the list is skipped.
  auto options = base_options("+0 1");
  options.photons = 3.0;
  const auto s = read_stream(stream_bytes(options));
  CHECK(count_type(s, obdeect::stream::RecordType::event) == 1);
  CHECK(s.rows.size() == 6);
  auto too_many = base_options("0 1");
  too_many.require = 3;
  CHECK(led::build_sources(too_many).empty());
  bool threw = false;
  try {
    auto bad_pixel = base_options("0 9");
    (void)led::build_sources(bad_pixel);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  CHECK(threw);
}

void non_contiguous_ids() {
  // IDs 20 and 10 are listed out of order and are not 0..N-1: each ID keeps its own
  // position (3 cm and 0 cm), and the labels are the IDs.
  const char *text = "Pixel 20 0 3.0 0.0 1.0\n"
                     "Pixel 10 0 0.0 0.0 1.0\n";
  led::Options options;
  std::istringstream in(text);
  options.camera = led::read_camera(in, options.slaves);
  options.use = "10 20";
  options.photons = 10.0;
  const auto s = read_stream(stream_bytes(options));
  CHECK(s.rows.size() == 20);
  int at_origin = 0;
  int at_3cm = 0;
  for (const Row &r : s.rows) {
    if (std::abs(r.x) < 1e-6) {
      ++at_origin;
    } else if (std::abs(r.x - 0.03) < 1e-6) {
      ++at_3cm;
    }
  }
  CHECK(at_origin == 10);
  CHECK(at_3cm == 10);
  const auto sources = led::build_sources(options);
  CHECK(sources.size() == 1);
  CHECK(sources[0].labels == std::vector<int>({10, 20}));

  bool threw = false;
  try {
    std::istringstream dup("Pixel 5 0 0.0 0.0 1.0\nPixel 5 0 1.0 0.0 1.0\n");
    (void)led::read_camera(dup, false);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  CHECK(threw);
}

void spectrum_wavelengths() {
  // Plain numbers in nm, with a decimal point or an exponent, give a fixed wavelength.
  for (const auto &[text, nm] : std::vector<std::pair<std::string, double>>{
           {"400", 400.0}, {"400.5", 400.5}, {"4e2", 400.0}}) {
    auto options = base_options("0");
    options.spectrum = text;
    const auto s = read_stream(stream_bytes(options));
    CHECK(s.rows.size() == 10);
    for (const Row &r : s.rows) {
      CHECK_NEAR(r.wavelength, nm, 1e-3);
    }
  }
  // Out of range, not a number, or a file name: rejected.
  for (const char *text : {"1e5", "100", "400x", "table.dat"}) {
    auto options = base_options("0");
    options.spectrum = text;
    bool threw = false;
    try {
      (void)led::build_sources(options);
    } catch (const std::invalid_argument &) {
      threw = true;
    }
    CHECK(threw);
  }
}

} // namespace

int main() {
  return run_checks([] {
    combinations();
    slaves_and_trigger_groups();
    random_pixels();
    all_pixels_and_list_errors();
    non_contiguous_ids();
    spectrum_wavelengths();
  });
}
