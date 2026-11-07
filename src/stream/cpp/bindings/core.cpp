// Python.h must precede all standard headers (it sets feature-test macros such as
// _POSIX_C_SOURCE), so the nanobind headers come first.
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"
#include "gstrickd/simtel.hpp"
#include "gstrickd/stdin_buffer.hpp"
#include "gstrickd/stream.hpp"
#include "gstrickd/version.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

namespace nb = nanobind;

namespace {

// Copies n values into a new numpy array that owns its memory.
template <class T> nb::object to_numpy(const T *data, std::size_t n) {
  auto *owned = new T[n];
  if (n > 0) {
    std::memcpy(owned, data, n * sizeof(T));
  }
  nb::capsule owner(owned, [](void *p) noexcept { delete[] static_cast<T *>(p); });
  return nb::cast(nb::ndarray<nb::numpy, T, nb::ndim<1>>(owned, {n}, owner));
}

struct MetadataValueToPy {
  nb::object operator()(std::int64_t v) const { return nb::int_(v); }
  nb::object operator()(double v) const { return nb::float_(v); }
  nb::object operator()(const std::string &v) const { return nb::str(v.c_str(), v.size()); }
  nb::object operator()(const std::vector<double> &v) const { return to_numpy(v.data(), v.size()); }
};

nb::dict metadata_to_dict(const gstrickd::Metadata &metadata) {
  nb::dict out;
  for (const auto &[key, value] : metadata) {
    out[nb::str(key.c_str(), key.size())] = std::visit(MetadataValueToPy{}, value);
  }
  return out;
}

// Python values: bool/int -> int, float -> float, str -> str, list/tuple -> float array.
gstrickd::Metadata dict_to_metadata(const nb::dict &dict) {
  gstrickd::Metadata metadata;
  for (auto item : dict) {
    const auto key = nb::cast<std::string>(item.first);
    const nb::handle value = item.second;
    if (nb::isinstance<nb::int_>(value)) {
      metadata[key] = nb::cast<std::int64_t>(value);
    } else if (nb::isinstance<nb::float_>(value)) {
      metadata[key] = nb::cast<double>(value);
    } else if (nb::isinstance<nb::str>(value)) {
      metadata[key] = nb::cast<std::string>(value);
    } else if (nb::isinstance<nb::list>(value) || nb::isinstance<nb::tuple>(value)) {
      metadata[key] = nb::cast<std::vector<double>>(value);
    } else {
      throw nb::type_error(("unsupported metadata value for key '" + key + "'").c_str());
    }
  }
  return metadata;
}

std::optional<gstrickd::Column> column_from_name(const std::string &name) {
  for (std::uint32_t i = 0; i < gstrickd::n_columns; ++i) {
    const auto column = static_cast<gstrickd::Column>(i);
    if (name == gstrickd::column_name(column)) {
      return column;
    }
  }
  return std::nullopt;
}

gstrickd::RecordType metadata_type_from_name(const std::string &name) {
  using gstrickd::RecordType;
  for (const auto type :
       {RecordType::run, RecordType::event, RecordType::telescope, RecordType::telescope_end,
        RecordType::event_end, RecordType::run_end, RecordType::ledger}) {
    if (name == gstrickd::record_name(type)) {
      return type;
    }
  }
  throw nb::value_error(("not a metadata record type: '" + name + "'").c_str());
}

class PyReader {
public:
  explicit PyReader(const std::string &path) {
    if (path == "-") {
      input_ = &std::cin;
    } else {
      file_ = std::make_unique<std::ifstream>(path, std::ios::binary);
      if (!*file_) {
        throw std::runtime_error("cannot open " + path);
      }
      input_ = file_.get();
    }
    reader_ = std::make_unique<gstrickd::Reader>(*input_);
  }

  nb::dict next() {
    if (!reader_->next(record_)) {
      if (!reader_->ended()) {
        throw std::runtime_error("gstrickd: truncated stream (no stream_end record)");
      }
      throw nb::stop_iteration();
    }
    nb::dict out;
    const gstrickd::RecordType type = record_.type;
    out["type"] = nb::str(gstrickd::record_name(type));
    switch (type) {
    case gstrickd::RecordType::photons: {
      const auto &chunk = record_.chunk;
      out["telescope_id"] = nb::int_(chunk.telescope_id);
      out["frame"] = nb::int_(static_cast<std::uint32_t>(chunk.frame));
      out["first_index"] = nb::int_(chunk.first_index);
      nb::dict columns;
      for (std::uint32_t i = 0; i < gstrickd::n_columns; ++i) {
        if (chunk.mask() & (1U << i)) {
          const auto &values = chunk.columns[i];
          const auto column = static_cast<gstrickd::Column>(i);
          columns[gstrickd::column_name(column)] = to_numpy(values.data(), values.size());
        }
      }
      out["columns"] = columns;
      break;
    }
    case gstrickd::RecordType::blob:
      out["eventio_type"] = nb::int_(record_.blob_eventio_type);
      out["eventio_id"] = nb::int_(record_.blob_eventio_id);
      out["data"] =
          nb::bytes(reinterpret_cast<const char *>(record_.blob.data()), record_.blob.size());
      break;
    default:
      out["metadata"] = metadata_to_dict(record_.metadata);
      break;
    }
    return out;
  }

private:
  std::unique_ptr<std::ifstream> file_;
  std::istream *input_{nullptr};
  std::unique_ptr<gstrickd::Reader> reader_;
  gstrickd::Record record_;
};

class PyWriter {
public:
  PyWriter(const std::string &path, const std::string &producer) {
    file_ = std::make_unique<std::ofstream>(path, std::ios::binary);
    if (!*file_) {
      throw std::runtime_error("cannot open " + path + " for writing");
    }
    writer_ = std::make_unique<gstrickd::Writer>(*file_, producer);
  }

  void write(const std::string &record_type, const nb::dict &metadata) {
    writer_->write_metadata(metadata_type_from_name(record_type), dict_to_metadata(metadata));
  }

  void write_photons(std::uint32_t telescope_id, std::uint32_t frame, std::uint64_t first_index,
                     const nb::dict &columns) {
    using Array = nb::ndarray<const float, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
    gstrickd::PhotonChunk chunk;
    if (frame > static_cast<std::uint32_t>(gstrickd::Frame::camera)) {
      throw nb::value_error("frame must be 0 (sky), 1 (ground), 2 (telescope) or 3 (camera)");
    }
    chunk.telescope_id = telescope_id;
    chunk.frame = static_cast<gstrickd::Frame>(frame);
    chunk.first_index = first_index;

    std::vector<std::pair<gstrickd::Column, Array>> items;
    std::size_t rows = 0;
    std::uint32_t mask = 0;
    for (auto item : columns) {
      const auto name = nb::cast<std::string>(item.first);
      const auto column = column_from_name(name);
      if (!column) {
        throw nb::value_error(("unknown photon column '" + name + "'").c_str());
      }
      Array array;
      if (!nb::try_cast(nb::handle(item.second), array, false)) {
        throw nb::type_error(
            ("column '" + name + "' must be a 1-D C-contiguous numpy array of float32").c_str());
      }
      if (items.empty()) {
        rows = array.shape(0);
      } else if (array.shape(0) != rows) {
        throw nb::value_error("photon columns differ in length");
      }
      mask |= gstrickd::column_bit(*column);
      items.emplace_back(*column, array);
    }
    chunk.resize(rows, mask);
    for (const auto &[column, array] : items) {
      std::copy_n(array.data(), rows, chunk[column].data());
    }
    writer_->write_chunk(chunk);
  }

  void write_blob(std::uint32_t eventio_type, std::uint32_t eventio_id, const nb::bytes &data) {
    const auto *begin = reinterpret_cast<const std::uint8_t *>(data.c_str());
    writer_->write_blob(eventio_type, eventio_id,
                        std::vector<std::uint8_t>(begin, begin + data.size()));
  }

  void close() {
    writer_->finish();
    if (file_->is_open()) {
      file_->close();
      if (file_->fail()) {
        throw std::runtime_error("gstrickd: writing the photon stream failed (e.g. disk full)");
      }
    }
  }

private:
  std::unique_ptr<std::ofstream> file_;
  std::unique_ptr<gstrickd::Writer> writer_;
};

// Input path or "-" (stdin) as a stream; the file is opened in binary mode.
class InputStream {
public:
  explicit InputStream(const std::string &path) {
    if (path == "-") {
      stream_ = &stdin_;
      return;
    }
    file_.open(path, std::ios::binary);
    if (!file_) {
      throw nb::value_error(("cannot open " + path).c_str());
    }
    stream_ = &file_;
  }
  std::istream &get() { return *stream_; }

private:
  gstrickd::StdinBuffer stdin_buffer_;
  std::istream stdin_{&stdin_buffer_};
  std::ifstream file_;
  std::istream *stream_{nullptr};
};

std::string producer_name() {
  return "corsika2ps " + std::to_string(gstrickd::version_major) + "." +
         std::to_string(gstrickd::version_minor) + "." + std::to_string(gstrickd::version_patch);
}

using CameraEntries = std::vector<std::tuple<std::uint32_t, int, double>>;

} // namespace

NB_MODULE(_core, module) {
  module.doc() = "gstrickd native core";
  module.def("core_version", [] {
    return std::make_tuple(gstrickd::version_major, gstrickd::version_minor,
                           gstrickd::version_patch);
  });
  module.def("column_names", [] {
    std::vector<std::string> names;
    for (std::uint32_t i = 0; i < gstrickd::n_columns; ++i) {
      names.emplace_back(gstrickd::column_name(static_cast<gstrickd::Column>(i)));
    }
    return names;
  });

  module.def(
      "corsika_to_stream",
      [](const std::string &input, const std::string &output) {
        InputStream in{input};
        std::ofstream out(output, std::ios::binary);
        if (!out) {
          throw nb::value_error(("cannot open " + output + " for writing").c_str());
        }
        gstrickd::corsika_to_stream(in.get(), out, producer_name());
      },
      nb::arg("input_path"), nb::arg("output_path"));
  module.def(
      "stream_to_simtel",
      [](const std::string &input, const std::string &output, const CameraEntries &camera) {
        gstrickd::CameraMap map;
        for (const auto &[id, mirror_class, rotation] : camera) {
          static_cast<void>(gstrickd::bypass_optics_level(mirror_class)); // validates
          if (!map.emplace(id, gstrickd::CameraTelescope{mirror_class, rotation}).second) {
            throw nb::value_error(("telescope " + std::to_string(id) + " listed twice").c_str());
          }
        }
        InputStream in{input};
        std::ofstream out(output, std::ios::binary);
        if (!out) {
          throw nb::value_error(("cannot open " + output + " for writing").c_str());
        }
        gstrickd::stream_to_simtel(in.get(), out, std::move(map));
      },
      nb::arg("input_path"), nb::arg("output_path"), nb::arg("camera"));
  module.def(
      "parse_camera_map",
      [](const std::string &text) {
        std::istringstream in(text);
        CameraEntries entries;
        for (const auto &[id, telescope] : gstrickd::read_camera_map(in)) {
          entries.emplace_back(id, telescope.mirror_class, telescope.rotation_deg);
        }
        return entries;
      },
      nb::arg("text"));

  nb::class_<PyReader>(module, "Reader")
      .def(nb::init<const std::string &>(), nb::arg("path"))
      .def("__iter__", [](nb::object self) { return self; })
      .def("__next__", &PyReader::next);

  nb::class_<PyWriter>(module, "Writer")
      .def(nb::init<const std::string &, const std::string &>(), nb::arg("path"),
           nb::arg("producer"))
      .def("write", &PyWriter::write, nb::arg("record_type"), nb::arg("metadata"))
      .def("write_photons", &PyWriter::write_photons, nb::arg("telescope_id"), nb::arg("frame"),
           nb::arg("first_index"), nb::arg("columns"))
      .def("write_blob", &PyWriter::write_blob, nb::arg("eventio_type"), nb::arg("eventio_id"),
           nb::arg("data"))
      .def("close", &PyWriter::close);
}
