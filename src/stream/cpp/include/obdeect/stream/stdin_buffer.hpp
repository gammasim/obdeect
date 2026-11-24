#pragma once

// Buffered std::streambuf over a file descriptor (default: standard input). Each read returns what
// the pipe has (partial reads), so a pipeline streams. std::cin is much slower for binary input on
// some standard libraries (libc++ reads it one character at a time after sync_with_stdio(false)),
// so the "-" input of the tools uses this buffer.

#include <cerrno>
#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <vector>

#include <unistd.h>

namespace obdeect::stream {

// A read error of the input (not the end of the input): the C API maps it to an I/O error.
class InputReadError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class StdinBuffer : public std::streambuf {
public:
  explicit StdinBuffer(std::size_t size = std::size_t{1} << 20, int fd = 0)
      : buffer_(size), fd_(fd) {
    char *end = buffer_.data() + buffer_.size();
    setg(end, end, end);
  }

protected:
  // Returns eof at the end of the input. A read error throws: the std::istream reading from this
  // buffer turns the exception into its bad state, which the readers report as an error (not as
  // a truncated or complete stream). Only EINTR is retried.
  int_type underflow() override {
    if (gptr() < egptr()) {
      return traits_type::to_int_type(*gptr());
    }
    for (;;) {
      const auto n = static_cast<long long>(::read(fd_, buffer_.data(), buffer_.size()));
      if (n > 0) {
        setg(buffer_.data(), buffer_.data(), buffer_.data() + n);
        return traits_type::to_int_type(*gptr());
      }
      if (n == 0) {
        return traits_type::eof();
      }
      if (errno != EINTR) {
        throw InputReadError(std::string("read error on the input: ") + std::strerror(errno));
      }
    }
  }

private:
  std::vector<char> buffer_;
  int fd_;
};

} // namespace obdeect::stream
