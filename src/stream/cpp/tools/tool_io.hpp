#pragma once

// Input and output of the stream tools. "-" is standard input or output (buffered); any other name
// is a file. Output is finished with finish(), which flushes and reports every write error,
// including those of the final flush and of close (a deferred error must not end with exit status
// 0). Standard output is also checked with fflush and ferror.

#include "obdeect/stream/stdin_buffer.hpp"

#include <cstddef>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace obdeect::stream::tool {

inline constexpr std::size_t io_buffer_size = std::size_t{1} << 20;

// The one standard input of the process. Read errors throw (the stream sets its exceptions mask,
// so they are not taken for the end of the input). Its buffer is kept between uses, so a tool
// that reads "-" twice does not lose the bytes buffered by the first read.
inline std::istream &standard_input() {
  static StdinBuffer buffer;
  static std::istream stream{&buffer};
  stream.exceptions(std::ios::badbit);
  return stream;
}

// Input from a file, or from standard input for "-".
class Input {
public:
  explicit Input(const std::string &name) {
    if (name == "-") {
      in_ = &standard_input();
      return;
    }
    buffer_.resize(io_buffer_size);
    file_.rdbuf()->pubsetbuf(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
    file_.open(name, std::ios::binary);
    if (!file_) {
      throw std::runtime_error("cannot open " + name);
    }
    file_.exceptions(std::ios::badbit);
    in_ = &file_;
  }

  [[nodiscard]] std::istream &stream() { return *in_; }

private:
  std::vector<char> buffer_;
  std::ifstream file_;
  std::istream *in_{nullptr};
};

// Output to a file, or to standard output for "-". Create it, write to stream(), then finish().
class Output {
public:
  explicit Output(const std::string &name) : name_{name} {
    if (name == "-") {
      buffer_.resize(io_buffer_size);
      std::cout.rdbuf()->pubsetbuf(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
      out_ = &std::cout;
      return;
    }
    buffer_.resize(io_buffer_size);
    file_.rdbuf()->pubsetbuf(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
    file_.open(name, std::ios::binary);
    if (!file_) {
      throw std::runtime_error("cannot open " + name + " for writing");
    }
    out_ = &file_;
  }

  [[nodiscard]] std::ostream &stream() { return *out_; }

  // Flushes and closes the output. Throws std::runtime_error when any write failed.
  void finish() {
    out_->flush();
    if (out_ == &std::cout) {
      if (std::fflush(stdout) != 0 || std::ferror(stdout) != 0 || !std::cout) {
        throw std::runtime_error("write error on standard output");
      }
      return;
    }
    file_.close();
    if (!file_) {
      throw std::runtime_error("write error on " + name_);
    }
  }

private:
  std::string name_;
  std::vector<char> buffer_;
  std::ofstream file_;
  std::ostream *out_{nullptr};
};

} // namespace obdeect::stream::tool
