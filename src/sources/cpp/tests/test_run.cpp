#include "check.hpp"

#include "obdeect/sources/emit.hpp"
#include "obdeect/sources/run.hpp"
#include "obdeect/sources/source.hpp"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/stream.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

using namespace obdeect::sources;
using obdeect::stream::Column;
using obdeect::stream::Frame;
using obdeect::stream::Metadata;
using obdeect::stream::PhotonChunk;
using obdeect::stream::Reader;
using obdeect::stream::Record;
using obdeect::stream::RecordType;

namespace {

// ---- builders

Telescope tel(std::uint32_t id, Vec3 position, double radius_m) {
  return Telescope{id, position, radius_m};
}

Source with_geometry(Geometry geometry, double photons, Spectrum spectrum = Fixed{400.0}) {
  Source s;
  s.geometry = std::move(geometry);
  s.photons = photons;
  s.spectrum = std::move(spectrum);
  return s;
}

Source isotropic_at(Vec3 position, double photons, Spectrum spectrum = Fixed{400.0}) {
  return with_geometry(Point{position, Vec3{0.0, 0.0, 1.0}, Isotropic{}}, photons,
                       std::move(spectrum));
}

Emission make_emission(Source source, double bunch = 1.0, double offset_ns = 0.0,
                       bool to_plane = false) {
  return Emission{std::move(source), bunch, offset_ns, to_plane};
}

Event make_event(std::uint64_t repeat, std::vector<Emission> emissions, Metadata metadata = {}) {
  Event e;
  e.repeat = repeat;
  e.emissions = std::move(emissions);
  e.metadata = std::move(metadata);
  return e;
}

Run make_run(std::vector<Telescope> telescopes, std::vector<Event> events, std::uint64_t seed = 1) {
  Run run;
  run.seed = seed;
  run.observation_level_m = 2150.0;
  run.telescopes = std::move(telescopes);
  run.events = std::move(events);
  return run;
}

// Target of the auto rule for one geometry with the given telescopes.
Target auto_target(Geometry geometry, const std::vector<Telescope> &telescopes) {
  return resolve_auto_target(with_geometry(std::move(geometry), 1.0), telescopes);
}

// ---- writing and reading

std::string write_bytes(const Run &run, std::uint64_t max_cached = max_cached_event_bytes,
                        RunSummary *summary = nullptr) {
  std::ostringstream out;
  const RunSummary s = write_run(out, run, max_cached);
  if (summary)
    *summary = s;
  return out.str();
}

// What the writer throws for the run (empty when it does not throw).
template <class F> std::string error_of(F &&f) {
  try {
    f();
  } catch (const std::invalid_argument &e) {
    return e.what();
  }
  return {};
}

struct TelescopeBlock {
  Metadata meta;
  std::vector<PhotonChunk> chunks;
  Metadata end;
};

struct EventBlock {
  Metadata meta;
  std::vector<TelescopeBlock> telescopes;
};

struct RunBlock {
  std::vector<RecordType> order;
  Metadata meta;
  std::vector<EventBlock> events;
};

RunBlock parse_run(const std::string &bytes) {
  std::istringstream in(bytes);
  Reader reader(in);
  Record r;
  RunBlock run;
  while (reader.next(r)) {
    run.order.push_back(r.type);
    switch (r.type) {
    case RecordType::run:
      run.meta = r.metadata;
      break;
    case RecordType::event: {
      EventBlock e;
      e.meta = r.metadata;
      run.events.push_back(std::move(e));
      break;
    }
    case RecordType::telescope: {
      TelescopeBlock t;
      t.meta = r.metadata;
      if (run.events.empty())
        ++check_failures;
      else
        run.events.back().telescopes.push_back(std::move(t));
      break;
    }
    case RecordType::photons:
      if (!run.events.empty() && !run.events.back().telescopes.empty())
        run.events.back().telescopes.back().chunks.push_back(r.chunk);
      break;
    case RecordType::telescope_end:
      if (!run.events.empty() && !run.events.back().telescopes.empty())
        run.events.back().telescopes.back().end = r.metadata;
      break;
    default:
      break;
    }
  }
  CHECK(reader.ended());
  CHECK(reader.at_clean_end());
  return run;
}

const TelescopeBlock *telescope_of(const RunBlock &r, std::size_t event, std::size_t telescope) {
  if (event < r.events.size() && telescope < r.events[event].telescopes.size())
    return &r.events[event].telescopes[telescope];
  std::fprintf(stderr, "%s:%d: missing telescope block %zu/%zu\n", __FILE__, __LINE__, event,
               telescope);
  ++check_failures;
  return nullptr;
}

template <class T> T field(const Metadata &m, const char *key) {
  const auto it = m.find(key);
  if (it == m.end() || !std::holds_alternative<T>(it->second)) {
    std::fprintf(stderr, "%s:%d: missing or mistyped metadata key %s\n", __FILE__, __LINE__, key);
    ++check_failures;
    return T{};
  }
  return std::get<T>(it->second);
}

double weight_sum(const PhotonChunk &c) {
  double sum = 0.0;
  for (const float w : c[Column::w])
    sum += static_cast<double>(w);
  return sum;
}

// Every source column of two chunks is equal (zem is not compared).
void check_same_columns(const PhotonChunk &got, const PhotonChunk &want) {
  CHECK(got.rows() == want.rows());
  for (const Column c : {Column::x, Column::y, Column::z, Column::cx, Column::cy, Column::cz,
                         Column::t, Column::w, Column::wavelength})
    CHECK(got[c] == want[c]);
}

// The photons one emission gives in a run (event number `event`, first emission of its seed):
// candidates = photons / bunch, time offset added, weights x bunch. Telescope 0 is the key.
PhotonChunk expected_emission(const Source &source, std::uint64_t seed, std::uint32_t event,
                              double bunch, double offset_ns) {
  Source d = source;
  d.photons = source.photons / bunch;
  PhotonChunk c;
  emit(d, seed, event, 0, 0, static_cast<std::size_t>(d.photons), c);
  for (auto &t : c[Column::t])
    t = static_cast<float>(static_cast<double>(t) + offset_ns);
  for (auto &w : c[Column::w])
    w = static_cast<float>(static_cast<double>(w) * bunch);
  return c;
}

// ---- tests

void test_record_order_and_metadata() {
  Metadata tag;
  tag["tag"] = std::string("a");
  const Source first = isotropic_at({0.0, 0.0, 80.0}, 1000.0);
  const Source second = isotropic_at({0.0, 0.0, 60.0}, 500.0, UniformRange{300.0, 600.0});
  Run run = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0), tel(2, {30.0, 0.0, 0.0}, 5.0)},
                     {make_event(1, {make_emission(first), make_emission(second)}, tag)}, 42);
  run.number = 7;
  run.atmosphere_id = 3;
  run.metadata["note"] = std::string("x");

  RunSummary summary;
  const RunBlock r = parse_run(write_bytes(run, max_cached_event_bytes, &summary));
  const std::vector<RecordType> order = {
      RecordType::stream,  RecordType::run,     RecordType::event,         RecordType::telescope,
      RecordType::photons, RecordType::photons, RecordType::telescope_end, RecordType::telescope,
      RecordType::photons, RecordType::photons, RecordType::telescope_end, RecordType::event_end,
      RecordType::run_end};
  CHECK(r.order == order);

  const Metadata &m = r.meta;
  CHECK(field<std::int64_t>(m, "run_number") == 7);
  CHECK(field<std::int64_t>(m, "seed") == 42);
  CHECK(field<std::int64_t>(m, "atmosphere_id") == 3);
  CHECK(field<std::int64_t>(m, "n_telescopes") == 2);
  CHECK(field<std::int64_t>(m, "refraction_applied") == 0);
  CHECK(field<std::int64_t>(m, "curvature_applied") == 0);
  CHECK_NEAR(field<double>(m, "observation_level_m"), 2150.0, 1e-9);
  CHECK_NEAR(field<double>(m, "wavelength_min_nm"), 300.0, 1e-9);
  CHECK_NEAR(field<double>(m, "wavelength_max_nm"), 600.0, 1e-9);
  CHECK(field<std::string>(m, "note") == "x");
  const std::vector<double> xs = {0.0, 30.0};
  const std::vector<double> ys = {0.0, 0.0};
  const std::vector<double> zs = {0.0, 0.0};
  const std::vector<double> rs = {10.0, 5.0};
  CHECK(field<std::vector<double>>(m, "telescope_x_m") == xs);
  CHECK(field<std::vector<double>>(m, "telescope_y_m") == ys);
  CHECK(field<std::vector<double>>(m, "telescope_z_m") == zs);
  CHECK(field<std::vector<double>>(m, "telescope_r_m") == rs);

  CHECK(r.events.size() == 1);
  if (r.events.size() != 1)
    return;
  const EventBlock &ev = r.events[0];
  CHECK(field<std::int64_t>(ev.meta, "event_number") == 1);
  CHECK(field<std::string>(ev.meta, "tag") == "a");
  CHECK(ev.telescopes.size() == 2);
  for (std::size_t t = 0; t < ev.telescopes.size(); ++t) {
    const TelescopeBlock &tb = ev.telescopes[t];
    const auto id = static_cast<std::uint32_t>(t + 1);
    CHECK(field<std::int64_t>(tb.meta, "telescope_id") == static_cast<std::int64_t>(id));
    CHECK(field<std::int64_t>(tb.meta, "frame") == 0);
    CHECK(field<double>(tb.meta, "x_m") == (t == 0 ? 0.0 : 30.0));
    CHECK(field<double>(tb.meta, "r_m") == (t == 0 ? 10.0 : 5.0));
    CHECK(tb.chunks.size() == 2);
    if (tb.chunks.size() != 2)
      continue;
    CHECK(tb.chunks[0].rows() == 1000);
    CHECK(tb.chunks[1].rows() == 500);
    for (const auto &c : tb.chunks) {
      CHECK(c.telescope_id == id);
      CHECK(c.frame == Frame::sky);
      CHECK(c.has(Column::zem));
    }
    for (const float wl : tb.chunks[1][Column::wavelength]) {
      CHECK(wl >= 300.0F);
      CHECK(wl <= 600.0F);
    }
    const double total = weight_sum(tb.chunks[0]) + weight_sum(tb.chunks[1]);
    CHECK_NEAR(total, 1500.0, 1e-9);
    CHECK_NEAR(field<double>(tb.end, "photons"), total, 1e-9);
  }
  CHECK(summary.events == 1);
  CHECK(summary.rows == 3000);
  CHECK_NEAR(summary.weight, 3000.0, 1e-9);
}

void test_event_numbering() {
  Metadata tag;
  tag["tag"] = std::string("third");
  const Source src = isotropic_at({0.0, 0.0, 80.0}, 10.0);
  Run run = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)},
                     {make_event(2, {make_emission(src)}), make_event(0, {make_emission(src)}),
                      make_event(3, {make_emission(src)}, tag)});
  RunSummary summary;
  const RunBlock r = parse_run(write_bytes(run, max_cached_event_bytes, &summary));
  CHECK(r.events.size() == 5);
  CHECK(summary.events == 5);
  CHECK(summary.rows == 50);
  for (std::size_t i = 0; i < r.events.size(); ++i) {
    CHECK(field<std::int64_t>(r.events[i].meta, "event_number") ==
          static_cast<std::int64_t>(i + 1));
    // Events from the third entry carry its metadata; the first two do not.
    CHECK(r.events[i].meta.count("tag") == (i >= 2 ? 1U : 0U));
  }
}

void test_telescopes_share_photons() {
  Source disc = isotropic_at({0.0, 0.0, 50.0}, 2e5);
  disc.target = DiscTarget{Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}, 10.0};
  const Source point = isotropic_at({0.0, 0.0, 30.0}, 300.0);
  Run run = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0), tel(2, {5.0, -3.0, 0.0}, 4.0)},
                     {make_event(2, {make_emission(disc), make_emission(point)})}, 9);
  const RunBlock r = parse_run(write_bytes(run));
  CHECK(r.events.size() == 2);
  if (r.events.size() != 2)
    return;
  for (std::size_t e = 0; e < 2; ++e) {
    const TelescopeBlock *t1 = telescope_of(r, e, 0);
    const TelescopeBlock *t2 = telescope_of(r, e, 1);
    if (!t1 || !t2)
      return;
    CHECK(t1->chunks.size() == 2);
    CHECK(t1->chunks.size() == t2->chunks.size());
    for (std::size_t k = 0; k < t1->chunks.size() && k < t2->chunks.size(); ++k) {
      CHECK(t1->chunks[k].rows() > 0);
      CHECK(t2->chunks[k].telescope_id == 2);
      check_same_columns(t2->chunks[k], t1->chunks[k]);
    }
  }
  // Repeats of an event are independent realisations.
  const TelescopeBlock *a = telescope_of(r, 0, 0);
  const TelescopeBlock *b = telescope_of(r, 1, 0);
  if (a && b && !a->chunks.empty() && !b->chunks.empty())
    CHECK(a->chunks[0][Column::cx] != b->chunks[0][Column::cx]);
}

void test_cache_does_not_change_output() {
  const Source disc = [] {
    Source s = isotropic_at({0.0, 0.0, 50.0}, 2e5);
    s.target = DiscTarget{Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}, 10.0};
    return s;
  }();
  Run run = make_run(
      {tel(1, {0.0, 0.0, 0.0}, 10.0), tel(2, {30.0, 0.0, 0.0}, 5.0), tel(3, {0.0, 40.0, 0.0}, 8.0)},
      {make_event(2, {make_emission(disc),
                      make_emission(isotropic_at({0.0, 0.0, 80.0}, 1000.0), 4.0, 1.5)}),
       make_event(1, {make_emission(isotropic_at({0.0, 0.0, 70.0}, 700.0))})},
      5);
  const std::string cached = write_bytes(run);
  const std::string uncached = write_bytes(run, 0);
  CHECK(!cached.empty());
  CHECK(cached == uncached);
}

void test_seeds_bunch_and_offset() {
  const std::uint64_t seed = 42;
  const Source a = isotropic_at({0.0, 0.0, 80.0}, 1000.0);
  const Source b = isotropic_at({0.0, 0.0, 60.0}, 1000.0);
  // Emission a is k = 0 (event 1, bunch 2, offset 3 ns); emission b is k = 1 (event 2).
  Run run = make_run(
      {tel(1, {0.0, 0.0, 0.0}, 10.0)},
      {make_event(1, {make_emission(a, 2.0, 3.0)}), make_event(1, {make_emission(b)})}, seed);
  const RunBlock r = parse_run(write_bytes(run));
  const TelescopeBlock *ta = telescope_of(r, 0, 0);
  const TelescopeBlock *tb = telescope_of(r, 1, 0);
  if (ta && tb && ta->chunks.size() == 1 && tb->chunks.size() == 1) {
    check_same_columns(ta->chunks[0], expected_emission(a, source_seed(seed, 0), 1, 2.0, 3.0));
    check_same_columns(tb->chunks[0], expected_emission(b, source_seed(seed, 1), 2, 1.0, 0.0));
  } else {
    CHECK(false);
  }

  // Emission numbering follows event order: a repeat-0 event still takes k = 0.
  Run shifted =
      make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)},
               {make_event(0, {make_emission(a)}), make_event(1, {make_emission(b)})}, seed);
  const RunBlock s = parse_run(write_bytes(shifted));
  const TelescopeBlock *ts = telescope_of(s, 0, 0);
  if (ts && ts->chunks.size() == 1)
    check_same_columns(ts->chunks[0], expected_emission(b, source_seed(seed, 1), 1, 1.0, 0.0));
  else
    CHECK(false);

  // Another run seed gives other photons.
  Run other = run;
  other.seed = seed + 1;
  const RunBlock o = parse_run(write_bytes(other));
  const TelescopeBlock *to = telescope_of(o, 0, 0);
  CHECK(to && ta && !to->chunks.empty() && !ta->chunks.empty());
  if (to && ta && !to->chunks.empty() && !ta->chunks.empty())
    CHECK(to->chunks[0][Column::cx] != ta->chunks[0][Column::cx]);
}

void test_bunch_size() {
  Run run = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)},
                     {make_event(1, {make_emission(isotropic_at({0.0, 0.0, 80.0}, 1000.0), 4.0)})});
  RunSummary summary;
  const RunBlock r = parse_run(write_bytes(run, max_cached_event_bytes, &summary));
  const TelescopeBlock *t = telescope_of(r, 0, 0);
  if (t && t->chunks.size() == 1) {
    CHECK(t->chunks[0].rows() == 250);
    for (const float w : t->chunks[0][Column::w])
      CHECK(w == 4.0F);
    CHECK_NEAR(field<double>(t->end, "photons"), 1000.0, 1e-9);
  } else {
    CHECK(false);
  }
  CHECK(summary.rows == 250);
  CHECK_NEAR(summary.weight, 1000.0, 1e-9);

  // 1000 / 3 candidates is not an integer: the count rounds to 333 or 334, weights stay 3.
  Run third =
      make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)},
               {make_event(1, {make_emission(isotropic_at({0.0, 0.0, 80.0}, 1000.0), 3.0)})});
  const RunBlock q = parse_run(write_bytes(third));
  const TelescopeBlock *u = telescope_of(q, 0, 0);
  if (u && u->chunks.size() == 1) {
    const std::size_t rows = u->chunks[0].rows();
    CHECK(rows == 333 || rows == 334);
    CHECK_NEAR(weight_sum(u->chunks[0]), 3.0 * static_cast<double>(rows), 1e-6);
  } else {
    CHECK(false);
  }
}

void test_to_plane_z0() {
  // Photons from (0, 0, 80) along (0.6, 0, -0.8): 100 m to the plane z = 0 at x = 60 m.
  const Source down =
      with_geometry(Point{Vec3{0.0, 0.0, 80.0}, Vec3{0.6, 0.0, -0.8}, FixedAngle{0.0}}, 100.0);
  Run run = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)},
                     {make_event(1, {make_emission(down, 1.0, 0.0, true)})}, 3);
  const RunBlock r = parse_run(write_bytes(run));
  const TelescopeBlock *t = telescope_of(r, 0, 0);
  if (!t)
    return;
  if (t->chunks.size() != 1) {
    CHECK(false);
    return;
  }
  const PhotonChunk &c = t->chunks[0];
  CHECK(c.rows() == 100);
  const double flight_ns = 100.0 / speed_of_light_m_per_ns;
  for (std::size_t i = 0; i < c.rows(); ++i) {
    CHECK(c[Column::z][i] == 0.0F);
    CHECK_NEAR(c[Column::x][i], 60.0, 1e-3);
    CHECK_NEAR(c[Column::y][i], 0.0, 1e-3);
    CHECK_NEAR(c[Column::t][i], flight_ns, 1e-3);
    // zem is the observation level plus z, i.e. the level itself after the move.
    CHECK_NEAR(c[Column::zem][i], 2150.0, 1e-3);
  }

  // Upward photons cannot reach the plane z = 0.
  const Source up =
      with_geometry(Point{Vec3{0.0, 0.0, 80.0}, Vec3{0.0, 0.0, 1.0}, FixedAngle{0.0}}, 100.0);
  Run upward = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)},
                        {make_event(1, {make_emission(up, 1.0, 0.0, true)})});
  const std::string error = error_of([&] { write_bytes(upward); });
  CHECK(error.find("travel down") != std::string::npos);
}

void test_zem_and_frame() {
  Run sky = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)},
                     {make_event(1, {make_emission(isotropic_at({0.0, 0.0, 80.0}, 100.0))})});
  const RunBlock r = parse_run(write_bytes(sky));
  const TelescopeBlock *t = telescope_of(r, 0, 0);
  if (t && t->chunks.size() == 1) {
    const PhotonChunk &c = t->chunks[0];
    CHECK(c.frame == Frame::sky);
    CHECK(c.has(Column::zem));
    for (const float zem : c[Column::zem])
      CHECK_NEAR(zem, 2230.0, 1e-3);
  } else {
    CHECK(false);
  }

  Run telescope_frame = sky;
  telescope_frame.frame = Frame::telescope;
  const RunBlock q = parse_run(write_bytes(telescope_frame));
  const TelescopeBlock *u = telescope_of(q, 0, 0);
  if (u && u->chunks.size() == 1) {
    const PhotonChunk &c = u->chunks[0];
    CHECK(c.frame == Frame::telescope);
    CHECK(!c.has(Column::zem));
    CHECK(c.rows() == 100);
  } else {
    CHECK(false);
  }
}

void test_resolve_auto_target() {
  const std::vector<Telescope> one = {tel(1, {0.0, 0.0, 0.0}, 10.0)};
  const Vec3 up{0.0, 0.0, 1.0};

  // Inside 1.02 R = 10.2 m: no target. Just outside: a sphere of radius R.
  CHECK(std::holds_alternative<std::monostate>(
      auto_target(Point{Vec3{0.0, 0.0, 10.19}, up, Isotropic{}}, one)));
  const Target iso_far = auto_target(Point{Vec3{0.0, 0.0, 50.0}, up, Isotropic{}}, one);
  const auto *sphere = std::get_if<SphereTarget>(&iso_far);
  CHECK(sphere != nullptr);
  if (sphere) {
    CHECK_NEAR(sphere->centre.x, 0.0, 1e-12);
    CHECK_NEAR(sphere->centre.z, 0.0, 1e-12);
    CHECK_NEAR(sphere->radius_m, 10.0, 1e-12);
  }

  // Fixed angle: no target.
  CHECK(std::holds_alternative<std::monostate>(
      auto_target(Point{Vec3{0.0, 0.0, 50.0}, up, FixedAngle{0.1}}, one)));

  // Any other profile: a disc of radius R in the plane z = 0.
  const Target disc_target = auto_target(Point{Vec3{0.0, 0.0, 50.0}, up, FilledCone{0.5}}, one);
  const auto *disc = std::get_if<DiscTarget>(&disc_target);
  CHECK(disc != nullptr);
  if (disc) {
    CHECK_NEAR(disc->centre.z, 0.0, 1e-12);
    CHECK_NEAR(disc->normal.z, 1.0, 1e-12);
    CHECK_NEAR(disc->radius_m, 10.0, 1e-12);
  }

  // Non-isotropic point in the plane z = 0 outside the telescopes: rejected.
  const std::string error =
      error_of([&] { auto_target(Point{Vec3{50.0, 0.0, 0.0}, up, FilledCone{0.5}}, one); });
  CHECK(!error.empty());

  // Lines: a fixed angle or isotropic emission gets the sphere, other profiles none.
  const Target line_fixed = auto_target(
      Line{Vec3{0.0, 0.0, 100.0}, Vec3{0.0, 0.0, -1.0}, 50.0, 1.0, FixedAngle{0.0}}, one);
  CHECK(std::holds_alternative<SphereTarget>(line_fixed));
  const Target line_iso =
      auto_target(Line{Vec3{0.0, 0.0, 100.0}, Vec3{0.0, 0.0, -1.0}, 50.0, 1.0, Isotropic{}}, one);
  CHECK(std::holds_alternative<SphereTarget>(line_iso));
  const Target line_cone = auto_target(
      Line{Vec3{0.0, 0.0, 100.0}, Vec3{0.0, 0.0, -1.0}, 50.0, 1.0, FilledCone{0.5}}, one);
  CHECK(std::holds_alternative<std::monostate>(line_cone));

  // Centre = mean telescope position (30, 0) and (-10, 0) -> (10, 0); R = max(|p - C| + r) =
  // max(20 + 5, 20 + 2) = 25. A single telescope off the origin gets its own sphere.
  const std::vector<Telescope> two = {tel(1, {30.0, 0.0, 0.0}, 5.0),
                                      tel(2, {-10.0, 0.0, 0.0}, 2.0)};
  const Target r25 = auto_target(Point{Vec3{0.0, 0.0, 100.0}, up, Isotropic{}}, two);
  const auto *s25 = std::get_if<SphereTarget>(&r25);
  CHECK(s25 != nullptr);
  if (s25) {
    CHECK_NEAR(s25->centre.x, 10.0, 1e-12);
    CHECK_NEAR(s25->radius_m, 25.0, 1e-12);
  }
  const std::vector<Telescope> off = {tel(1, {100.0, 0.0, 0.0}, 12.5)};
  const Target own = auto_target(
      Line{Vec3{0.0, 0.0, 1000.0}, Vec3{0.0, 0.0, -1.0}, 1000.0, 1.0, FixedAngle{0.01}}, off);
  const auto *s_own = std::get_if<SphereTarget>(&own);
  CHECK(s_own != nullptr);
  if (s_own) {
    CHECK_NEAR(s_own->centre.x, 100.0, 1e-12);
    CHECK_NEAR(s_own->radius_m, 12.5, 1e-12);
  }
}

void test_validation() {
  const Source src = isotropic_at({0.0, 0.0, 80.0}, 10.0);
  const std::vector<Event> events = {make_event(1, {make_emission(src)})};

  // Each invalid run fails before anything is written.
  const auto rejected = [&](const Run &run, const char *needle) {
    std::ostringstream out;
    const std::string message = error_of([&] { write_run(out, run); });
    CHECK(message.find(needle) != std::string::npos);
    CHECK(out.str().empty());
  };

  rejected(make_run({}, events), "at least one telescope");
  rejected(make_run({tel(1, {0.0, 0.0, 0.0}, 10.0), tel(1, {30.0, 0.0, 0.0}, 5.0)}, events),
           "distinct");
  rejected(make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)}, {make_event(1, {make_emission(src, 0.5)})}),
           "bunch size");
  Run ground = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)}, events);
  ground.frame = Frame::ground;
  rejected(ground, "frame");
}

// Writes one emission and checks the weight sum against its expected value. The Monte Carlo
// error is estimated from the rows: n candidates, zero weight for dropped line photons.
void check_weight_matches_expected(const char *label, const Run &run) {
  RunSummary summary;
  const RunBlock r = parse_run(write_bytes(run, max_cached_event_bytes, &summary));
  double m1 = 0.0;
  double m2 = 0.0;
  std::uint64_t rows = 0;
  for (const auto &e : r.events)
    for (const auto &t : e.telescopes)
      for (const auto &c : t.chunks)
        for (const float w : c[Column::w]) {
          const auto wd = static_cast<double>(w);
          m1 += wd;
          m2 += wd * wd;
          ++rows;
        }
  const double n = expected_count(run.events[0].emissions[0].source);
  const double sigma = std::sqrt((m2 - m1 * m1 / n) * n / (n - 1.0));
  std::printf("%s: weight %.6g, expected %.6g, 4 sigma %.4g, rows %llu, candidates %.6g\n", label,
              summary.weight, summary.expected, 4.0 * sigma, static_cast<unsigned long long>(rows),
              n);
  CHECK(summary.events == 1);
  CHECK(summary.rows == rows);
  CHECK_NEAR(summary.weight, m1, 1e-9 * m1);
  CHECK(summary.expected > 0.0);
  CHECK(std::abs(summary.weight - summary.expected) <= 4.0 * sigma);
}

void test_summary_weight_matches_expected() {
  // A point with a disc target: the disc of R = 10 m at 50 m on the axis.
  Source disc = isotropic_at({0.0, 0.0, 50.0}, 2e6);
  disc.target = DiscTarget{Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}, 10.0};
  check_weight_matches_expected(
      "point, disc target",
      make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)}, {make_event(1, {make_emission(disc)})}));

  // A line with a sphere target: isotropic emission, and a 20 deg cone seen from a line offset
  // by 15 m (partial arcs, dropped misses).
  Source iso_line =
      with_geometry(Line{Vec3{0.0, 0.0, 100.0}, Vec3{0.0, 0.0, -1.0}, 80.0, 1.0, Isotropic{}}, 1e5);
  iso_line.target = SphereTarget{Vec3{0.0, 0.0, 0.0}, 10.0};
  check_weight_matches_expected(
      "line isotropic, sphere target",
      make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)}, {make_event(1, {make_emission(iso_line)})}));

  const double cone_u0 = 1.0 - std::cos(20.0 * std::acos(-1.0) / 180.0);
  Source cone_line = with_geometry(
      Line{Vec3{15.0, 0.0, 100.0}, Vec3{0.0, 0.0, -1.0}, 80.0, 1.0, FixedAngle{cone_u0}}, 1e5);
  cone_line.target = SphereTarget{Vec3{0.0, 0.0, 0.0}, 10.0};
  check_weight_matches_expected(
      "line 20 deg cone, sphere target",
      make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)}, {make_event(1, {make_emission(cone_line)})}));
}

} // namespace

int main() {
  return run_checks([] {
    test_record_order_and_metadata();
    test_event_numbering();
    test_telescopes_share_photons();
    test_cache_does_not_change_output();
    test_seeds_bunch_and_offset();
    test_bunch_size();
    test_to_plane_z0();
    test_zem_and_frame();
    test_resolve_auto_target();
    test_validation();
    test_summary_weight_matches_expected();
  });
}
