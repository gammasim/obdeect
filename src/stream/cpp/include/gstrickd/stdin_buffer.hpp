#pragma once

// Buffered std::streambuf over a file descriptor (default: standard input). Each read returns what
// the pipe has (partial reads), so a pipeline streams. std::cin is much slower for binary input on
// some standard libraries (libc++ reads it one character at a time after sync_with_stdio(false)),
// so the "-" input of the tools uses this buffer.

#include <cerrno>
#include <cstddef>
#include <streambuf>
#include <vector>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

namespace gstrickd {

class StdinBuffer : public std::streambuf {
public:
  explicit StdinBuffer(std::size_t size = std::size_t{1} << 20, int fd = 0)
      : buffer_(size), fd_(fd) {
    char *end = buffer_.data() + buffer_.size();
    setg(end, end, end);
  }

protected:
  int_type underflow() override {
    if (gptr() < egptr()) {
      return traits_type::to_int_type(*gptr());
    }
    for (;;) {
#ifdef _WIN32
      const auto n = static_cast<long long>(
          _read(fd_, buffer_.data(), static_cast<unsigned int>(buffer_.size())));
#else
      const auto n = static_cast<long long>(::read(fd_, buffer_.data(), buffer_.size()));
#endif
      if (n > 0) {
        setg(buffer_.data(), buffer_.data(), buffer_.data() + n);
        return traits_type::to_int_type(*gptr());
      }
      // End of input or a read error: the caller sees the end, and a truncated stream fails
      // its stream_end check.
      if (n == 0 || errno != EINTR) {
        return traits_type::eof();
      }
    }
  }

private:
  std::vector<char> buffer_;
  int fd_;
};

} // namespace gstrickd
