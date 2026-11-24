#include "check.hpp"

#include "obdeect/sources/emit.hpp"
#include "obdeect/sources/run.hpp"
#include "obdeect/sources/source.hpp"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/stream.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <initializer_list>
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
  Metadata end;
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
    case RecordType::run_end:
      run.end = r.metadata;
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

  // 1000 / 3 candidates is not an integer: the count is 333 for this seed, weights stay 3.
  Run third =
      make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)},
               {make_event(1, {make_emission(isotropic_at({0.0, 0.0, 80.0}, 1000.0), 3.0)})});
  const RunBlock q = parse_run(write_bytes(third));
  const TelescopeBlock *u = telescope_of(q, 0, 0);
  if (u && u->chunks.size() == 1) {
    const std::size_t rows = u->chunks[0].rows();
    // Deterministic for this seed: 1000 / 3 = 333.33 candidates, counted as 333 rows.
    CHECK(rows == 333);
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

  // Upward photons from above the plane move away from it: dropped and counted, not an error.
  const Source up =
      with_geometry(Point{Vec3{0.0, 0.0, 80.0}, Vec3{0.0, 0.0, 1.0}, FixedAngle{0.0}}, 100.0);
  Run upward = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)},
                        {make_event(1, {make_emission(up, 1.0, 0.0, true)})});
  RunSummary away_summary;
  const RunBlock away = parse_run(write_bytes(upward, max_cached_event_bytes, &away_summary));
  CHECK(away_summary.not_on_plane == 100);
  CHECK(away_summary.rows == 0);
  CHECK(away_summary.weight == 0.0);
  CHECK(away_summary.expected > 0.0);
  CHECK(field<std::int64_t>(away.end, "not_on_plane") == 100);
  const TelescopeBlock *none = telescope_of(away, 0, 0);
  if (none) {
    CHECK(none->chunks.empty());
    CHECK(field<double>(none->end, "photons") == 0.0);
  }

  // Upward photons from below the plane reach it: path 25 m along (0.6, 0, 0.8) from z = -20 m.
  const Source from_below =
      with_geometry(Point{Vec3{0.0, 0.0, -20.0}, Vec3{0.6, 0.0, 0.8}, FixedAngle{0.0}}, 100.0);
  Run below = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)},
                       {make_event(1, {make_emission(from_below, 1.0, 0.0, true)})}, 3);
  RunSummary below_summary;
  const RunBlock b = parse_run(write_bytes(below, max_cached_event_bytes, &below_summary));
  CHECK(below_summary.not_on_plane == 0);
  CHECK(below_summary.rows == 100);
  const TelescopeBlock *ub = telescope_of(b, 0, 0);
  if (ub && ub->chunks.size() == 1) {
    const PhotonChunk &c = ub->chunks[0];
    CHECK(c.rows() == 100);
    for (std::size_t i = 0; i < c.rows(); ++i) {
      CHECK(c[Column::z][i] == 0.0F);
      CHECK_NEAR(c[Column::x][i], 15.0, 1e-3);
      CHECK_NEAR(c[Column::y][i], 0.0, 1e-3);
      CHECK_NEAR(c[Column::t][i], 25.0 / speed_of_light_m_per_ns, 1e-3);
    }
  }

  // Horizontal photons (cz = 0) move away from the plane and are dropped as well.
  const Source horizontal =
      with_geometry(Point{Vec3{0.0, 0.0, 80.0}, Vec3{1.0, 0.0, 0.0}, FixedAngle{0.0}}, 10.0);
  RunSummary horizontal_summary;
  write_bytes(make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)},
                       {make_event(1, {make_emission(horizontal, 1.0, 0.0, true)})}),
              max_cached_event_bytes, &horizontal_summary);
  CHECK(horizontal_summary.not_on_plane == 10);
  CHECK(horizontal_summary.rows == 0);
}

// Isotropic light from above the plane: half the photons move down and reach it, the others are
// dropped. The kept rows are the downward candidates in emission order, projected onto the plane.
void test_isotropic_to_plane_keeps_half() {
  const double photons = 1e5;
  Source iso = isotropic_at({0.0, 0.0, 80.0}, photons);
  Run plane = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)},
                       {make_event(1, {make_emission(iso, 1.0, 0.0, true)})}, 11);
  RunSummary summary;
  const RunBlock p = parse_run(write_bytes(plane, max_cached_event_bytes, &summary));
  const double kept = static_cast<double>(summary.rows);
  const double dropped = static_cast<double>(summary.not_on_plane);
  CHECK(kept + dropped == photons);
  const double sigma = std::sqrt(photons * 0.25); // binomial, p = 1/2
  CHECK(std::abs(kept - photons / 2.0) <= 4.0 * sigma);
  CHECK(std::abs(dropped - photons / 2.0) <= 4.0 * sigma);
  CHECK_NEAR(summary.weight, kept, 1e-6 * kept);

  // The same candidates without the plane move: the kept ones, projected by hand, in order.
  Run all = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)},
                     {make_event(1, {make_emission(iso, 1.0, 0.0, false)})}, 11);
  const RunBlock a = parse_run(write_bytes(all));
  std::vector<std::array<double, 3>> expected;
  const TelescopeBlock *ta = telescope_of(a, 0, 0);
  if (ta)
    for (const PhotonChunk &c : ta->chunks)
      for (std::size_t i = 0; i < c.rows(); ++i) {
        const double z = c[Column::z][i], cz = c[Column::cz][i];
        if (!(z * cz < 0.0))
          continue;
        const double path = -z / cz;
        expected.push_back({c[Column::x][i] + path * c[Column::cx][i],
                            c[Column::y][i] + path * c[Column::cy][i],
                            c[Column::t][i] + path / speed_of_light_m_per_ns});
      }
  CHECK(expected.size() == summary.rows);
  const TelescopeBlock *tp = telescope_of(p, 0, 0);
  std::size_t row = 0;
  std::uint64_t next_index = 0;
  if (tp)
    for (const PhotonChunk &c : tp->chunks) {
      CHECK(c.first_index == next_index); // consecutive numbering over the kept rows
      next_index += c.rows();
      for (std::size_t i = 0; i < c.rows() && row < expected.size(); ++i, ++row) {
        CHECK(c[Column::z][i] == 0.0F);
        // float32 rows: relative tolerance (positions and times reach 1e5 m and 1e6 ns here).
        CHECK_NEAR(c[Column::x][i], expected[row][0], 1e-6 * (1.0 + std::abs(expected[row][0])));
        CHECK_NEAR(c[Column::y][i], expected[row][1], 1e-6 * (1.0 + std::abs(expected[row][1])));
        CHECK_NEAR(c[Column::t][i], expected[row][2], 1e-6 * (1.0 + std::abs(expected[row][2])));
      }
    }
  CHECK(row == expected.size());
  CHECK(next_index == summary.rows);

  // Cache on and off give the same bytes.
  CHECK(write_bytes(plane, 0) == write_bytes(plane));
}

// A point source below the telescope sphere emits upward into it: the written weight is the
// solid-angle fraction of the sphere, 1/2 (1 - sqrt(1 - (R/d)^2)), times the photons. The same
// holds for a source above it. Every row of the source below the sphere moves up.
void test_sphere_target_below_and_above() {
  const double photons = 1e6;
  const double radius = 10.0;
  const double distance = 50.0;
  const double fraction = 0.5 * (1.0 - std::sqrt(1.0 - (radius / distance) * (radius / distance)));
  const double mean = photons * fraction;
  for (const double z : {-distance, distance}) {
    Source s = isotropic_at({0.0, 0.0, z}, photons);
    s.target = SphereTarget{Vec3{0.0, 0.0, 0.0}, radius};
    const RunSummary summary = [&] {
      RunSummary out;
      write_bytes(make_run({tel(1, {0.0, 0.0, 0.0}, radius)}, {make_event(1, {make_emission(s)})}),
                  max_cached_event_bytes, &out);
      return out;
    }();
    std::printf("sphere target from z = %g: weight %.6g, analytic %.6g\n", z, summary.weight, mean);
    CHECK(std::abs(summary.weight - mean) <= 4.0 * std::sqrt(mean));
    CHECK(std::abs(summary.expected - mean) <= 1.0 + 1e-9 * mean);
    CHECK(summary.not_on_plane == 0);
    if (z < 0.0) {
      const RunBlock r = parse_run(write_bytes(
          make_run({tel(1, {0.0, 0.0, 0.0}, radius)}, {make_event(1, {make_emission(s)})})));
      const TelescopeBlock *t = telescope_of(r, 0, 0);
      if (t)
        for (const PhotonChunk &c : t->chunks)
          for (std::size_t i = 0; i < c.rows(); ++i)
            CHECK(c[Column::cz][i] > 0.0);
    }
  }
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

  // Inside the sphere (radius R = 10 m): no target. Just outside it: the sphere (no threshold, so
  // the weight sum is continuous across the boundary of the sphere's region).
  CHECK(std::holds_alternative<std::monostate>(
      auto_target(Point{Vec3{0.0, 0.0, 9.99}, up, Isotropic{}}, one)));
  CHECK(std::holds_alternative<SphereTarget>(
      auto_target(Point{Vec3{0.0, 0.0, 10.01}, up, Isotropic{}}, one)));
  CHECK(std::holds_alternative<SphereTarget>(
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
  Run ground = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0), tel(2, {30.0, 0.0, 0.0}, 5.0)}, events);
  ground.frame = Frame::ground;
  rejected(ground, "exactly one telescope");
  Run ground_to_plane = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)},
                                 {make_event(1, {make_emission(src, 1.0, 0.0, true)})});
  ground_to_plane.frame = Frame::ground;
  rejected(ground_to_plane, "to_plane_z0");
  // The seed is stored as a signed 64-bit integer: 2^63 - 1 is the largest accepted seed.
  rejected(make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)}, events, std::uint64_t{1} << 63),
           "seed must be at most 2^63 - 1");
  rejected(make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)}, events, ~std::uint64_t{0}),
           "seed must be at most 2^63 - 1");
  std::ostringstream ok;
  write_run(ok, make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)}, events, (std::uint64_t{1} << 63) - 1));
  CHECK(!ok.str().empty());
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

// ---- ground runs (sphere representation)

// Rows of a ground run: one telescope block, its chunks concatenated and checked by
// validate_sphere_rows. The run must have one event.
std::vector<PhotonChunk> sphere_chunks(const RunBlock &r, double radius_m) {
  std::vector<PhotonChunk> chunks;
  const TelescopeBlock *t = telescope_of(r, 0, 0);
  if (!t)
    return chunks;
  for (const PhotonChunk &c : t->chunks) {
    CHECK(c.frame == Frame::ground);
    CHECK(error_of([&] { obdeect::stream::validate_sphere_rows(c, radius_m); }).empty());
    chunks.push_back(c);
  }
  return chunks;
}

// A line of photons along +x from x = -20 m to 20 m at the origin of the frame, sphere R = 10 m,
// telescope z = 100 m (emission altitude = 2150 + 100 + z of the emission point). Photons
// emitted inside the sphere give rows at their emission points (path 0); photons emitted at
// x < -10 m give entry points at x = -10 m with path -10 - x; photons at x > 10 m move away and
// miss. The muon-like time is s / c along the line, so every entry point has t = 10 m / c.
void test_ground_line_crossing_the_sphere() {
  const double c = speed_of_light_m_per_ns;
  Source track = with_geometry(
      Line{Vec3{-20.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, 40.0, 1.0, FixedAngle{0.0}}, 1000.0);
  Run run =
      make_run({tel(1, {0.0, 0.0, 100.0}, 10.0)}, {make_event(1, {make_emission(track)})}, 21);
  run.frame = Frame::ground;
  RunSummary summary;
  const RunBlock r = parse_run(write_bytes(run, max_cached_event_bytes, &summary));
  CHECK(field<std::string>(r.meta, "ground_representation") == "sphere");
  CHECK(summary.rows + summary.misses_sphere == 1000);
  CHECK(field<std::int64_t>(r.end, "misses_sphere") ==
        static_cast<std::int64_t>(summary.misses_sphere));
  CHECK(summary.not_on_plane == 0);
  // Expected: 500 inside (|x| < 10), 250 entries, 250 misses; binomial with p = 1/4 and 3/4.
  CHECK(std::abs(static_cast<double>(summary.misses_sphere) - 250.0) <= 4.0 * std::sqrt(187.5));
  const auto chunks = sphere_chunks(r, 10.0);
  std::size_t inside = 0;
  std::size_t entries = 0;
  for (const PhotonChunk &ch : chunks) {
    for (std::size_t i = 0; i < ch.rows(); ++i) {
      const double x = ch[Column::x][i];
      const double path = ch[Column::path][i];
      CHECK_NEAR(ch[Column::cx][i], 1.0, 1e-6);
      CHECK_NEAR(ch[Column::zem][i], 2150.0 + 100.0, 1e-3);
      if (path == 0.0) {
        // Inside: the emission point itself.
        ++inside;
        CHECK(std::abs(x) < 10.0);
        CHECK_NEAR(ch[Column::y][i], 0.0, 1e-4);
        CHECK_NEAR(ch[Column::z][i], 0.0, 1e-4);
        CHECK_NEAR(ch[Column::t][i], (x + 20.0) / c, 1e-3);
      } else {
        // Entry point: on the sphere at x = -10 m, reached after a path of -10 - x_emission.
        ++entries;
        CHECK_NEAR(x, -10.0, 1e-4);
        CHECK_NEAR(ch[Column::y][i], 0.0, 1e-4);
        CHECK(path > 0.0 && path < 10.0);
        CHECK_NEAR(ch[Column::t][i], 10.0 / c, 1e-3);
      }
    }
  }
  CHECK(inside + entries == summary.rows);
  CHECK(std::abs(static_cast<double>(inside) - 500.0) <= 4.0 * std::sqrt(187.5));
  CHECK(entries > 0);
}

// Isotropic light from a point at D = 50 m outside the sphere (R = 10 m): the kept fraction is the
// solid-angle fraction (1 - sqrt(1 - (R/D)^2)) / 2. Each entry point lies on the sphere on the
// photon's line, at path s with time s / c; zem is the point's emission altitude.
void test_ground_isotropic_outside_fraction() {
  const double c = speed_of_light_m_per_ns;
  const double photons = 2e5;
  const double radius = 10.0;
  const double distance = 50.0;
  Source iso = isotropic_at({0.0, 0.0, distance}, photons);
  Run run = make_run({tel(1, {0.0, 0.0, 0.0}, radius)}, {make_event(1, {make_emission(iso)})}, 31);
  run.frame = Frame::ground;
  RunSummary summary;
  const RunBlock r = parse_run(write_bytes(run, max_cached_event_bytes, &summary));
  const double fraction = 0.5 * (1.0 - std::sqrt(1.0 - (radius / distance) * (radius / distance)));
  const double expected = photons * fraction;
  const double sigma = std::sqrt(photons * fraction * (1.0 - fraction));
  const double kept = static_cast<double>(summary.rows);
  CHECK(summary.rows + summary.misses_sphere == 200000);
  CHECK(std::abs(kept - expected) <= 4.0 * sigma);
  CHECK_NEAR(summary.weight, kept, 1e-6 * kept);
  CHECK_NEAR(summary.expected, photons, 1e-6 * photons);
  const auto chunks = sphere_chunks(r, radius);
  for (const PhotonChunk &ch : chunks)
    for (std::size_t i = 0; i < ch.rows(); ++i) {
      const Vec3 row{ch[Column::x][i], ch[Column::y][i], ch[Column::z][i]};
      const Vec3 d{ch[Column::cx][i], ch[Column::cy][i], ch[Column::cz][i]};
      const double path = ch[Column::path][i];
      CHECK_NEAR(norm(row), radius, 1e-4);
      CHECK(path > 0.0);
      CHECK_NEAR(ch[Column::t][i], path / c, 1e-3);
      CHECK_NEAR(ch[Column::zem][i], 2150.0 + distance, 1e-3);
      // The entry point lies on the photon's line from the source: source + path d = row.
      const Vec3 source{0.0, 0.0, distance};
      CHECK_NEAR(norm(source + d * path - row), 0.0, 1e-3);
    }
}

// Auto target in a ground run: the isotropic point aims at the sphere at the frame origin (the
// candidates are the photons that reach it), so every candidate hits and none is missed.
void test_ground_auto_target_hits_the_sphere() {
  Source iso = isotropic_at({0.0, 0.0, 50.0}, 1e6);
  iso.target = AutoTarget{};
  Run run = make_run({tel(1, {0.0, 0.0, 0.0}, 10.0)}, {make_event(1, {make_emission(iso)})}, 41);
  run.frame = Frame::ground;
  const double candidates = expected_count(plan(run)[0][0].emission.source);
  RunSummary summary;
  const RunBlock r = parse_run(write_bytes(run, max_cached_event_bytes, &summary));
  CHECK(summary.misses_sphere == 0);
  CHECK(std::abs(static_cast<double>(summary.rows) - candidates) <= 1.0);
  CHECK(std::abs(summary.weight - summary.expected) <= 4.0 * std::sqrt(candidates));
  CHECK(sphere_chunks(r, 10.0).size() > 0);
}

// Rows of a ground run do not depend on the chunk size: emitting the candidates in one chunk or in
// two, and finishing each chunk, gives the same rows in the same order.
// Rows of a ground run lie on the telescope sphere to float32 rounding, also for a photon that
// starts tens of km away (the stored direction cosines have a norm error of about 1e-7, which
// used to move the entry point by millimetres over such a path).
void test_ground_rows_on_the_sphere() {
  const double radius = 30.0;
  // A line 20 km away from the frame origin, crossing the sphere; its direction is not exact in
  // float32 (the cosines are rounded), and the line starts 20 km before the sphere.
  Line track;
  track.start = {-20000.0, -10.0, 5.0};
  track.direction = normalised({1.0, 0.0005, 0.0003});
  track.length_m = 40000.0;
  track.angular = FixedAngle{0.0};
  Source point = with_geometry(track, 1e5);
  point.target = SphereTarget{Vec3{0.0, 0.0, 0.0}, radius};
  Run run =
      make_run({tel(1, {0.0, 0.0, 0.0}, radius)}, {make_event(1, {make_emission(point)})}, 52);
  run.frame = Frame::ground;
  const RunBlock r = parse_run(write_bytes(run));
  double worst_entry = 0.0;
  double worst_outside = 0.0;
  std::size_t entries = 0;
  std::size_t rows = 0;
  for (const PhotonChunk &c : sphere_chunks(r, radius)) {
    for (std::size_t i = 0; i < c.rows(); ++i) {
      const Vec3 row{c[Column::x][i], c[Column::y][i], c[Column::z][i]};
      ++rows;
      worst_outside = std::max(worst_outside, norm(row) - radius);
      if (c[Column::path][i] > 0.0F) {
        // An entry point: on the sphere to float32 rounding (ulp(30 m) = 3.8e-6 m).
        worst_entry = std::max(worst_entry, std::abs(norm(row) - radius));
        ++entries;
      }
    }
  }
  CHECK(rows > 1000 && entries > 1000);
  CHECK(worst_entry <= 4.0 * 3.8e-6);
  CHECK(worst_outside <= 4.0 * 3.8e-6);
}

void test_ground_rows_independent_of_chunk() {
  Source iso = isotropic_at({3.0, -4.0, 40.0}, 1e3);
  const Emission e = make_emission(iso, 2.0, 5.0);
  Run run = make_run({tel(1, {0.0, 0.0, 30.0}, 10.0)}, {make_event(1, {e})}, 51);
  run.frame = Frame::ground;
  const std::size_t n = 1000;
  const auto rows_of = [&](std::initializer_list<std::size_t> sizes) {
    std::vector<PhotonChunk> out;
    std::uint64_t first = 0;
    for (const std::size_t rows : sizes) {
      PhotonChunk chunk;
      Drops drops;
      emit(iso, 51, 1, 1, first, rows, chunk);
      detail::finish_rows(e, run, chunk, drops);
      out.push_back(chunk);
      first += rows;
    }
    return out;
  };
  const auto whole = rows_of({n});
  const auto split = rows_of({300, 1, 699});
  std::vector<float> a_x, a_t, a_w, a_path;
  std::vector<float> b_x, b_t, b_w, b_path;
  for (const auto &c : whole)
    for (std::size_t i = 0; i < c.rows(); ++i) {
      a_x.push_back(c[Column::x][i]);
      a_t.push_back(c[Column::t][i]);
      a_w.push_back(c[Column::w][i]);
      a_path.push_back(c[Column::path][i]);
    }
  for (const auto &c : split)
    for (std::size_t i = 0; i < c.rows(); ++i) {
      b_x.push_back(c[Column::x][i]);
      b_t.push_back(c[Column::t][i]);
      b_w.push_back(c[Column::w][i]);
      b_path.push_back(c[Column::path][i]);
    }
  CHECK(!a_x.empty());
  CHECK(a_x == b_x);
  CHECK(a_t == b_t);
  CHECK(a_w == b_w);
  CHECK(a_path == b_path);
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
    test_isotropic_to_plane_keeps_half();
    test_sphere_target_below_and_above();
    test_zem_and_frame();
    test_resolve_auto_target();
    test_validation();
    test_summary_weight_matches_expected();
    test_ground_line_crossing_the_sphere();
    test_ground_isotropic_outside_fraction();
    test_ground_auto_target_hits_the_sphere();
    test_ground_rows_independent_of_chunk();
    test_ground_rows_on_the_sphere();
  });
}
