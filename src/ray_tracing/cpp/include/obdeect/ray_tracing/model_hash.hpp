#pragma once

// SHA-256 of an optical-model document, as compiler.model_hash in Python: the hash of the
// canonical JSON of the top-level object without its "sha256" entry. Canonical JSON is what
// Python's json.dumps(..., ensure_ascii=False, sort_keys=True, separators=(",", ":")) writes for
// the parsed document: keys sorted by code point, no spaces, integers as written, floats as
// Python's repr, strings with Python's escapes. Standard library only.

#include "obdeect/ray_tracing/model.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace obdeect::ray_tracing {

namespace hash_detail {

// FIPS 180-4 SHA-256.
class Sha256 {
public:
  void update(std::string_view data) {
    for (const char c : data) {
      block_[used_++] = static_cast<std::uint8_t>(c);
      if (used_ == 64) {
        compress();
        used_ = 0;
      }
    }
    length_ += data.size();
  }

  [[nodiscard]] std::string hex() {
    const std::uint64_t bits = length_ * 8U;
    block_[used_++] = 0x80;
    if (used_ > 56) {
      std::fill(block_.begin() + static_cast<std::ptrdiff_t>(used_), block_.end(), 0);
      compress();
      used_ = 0;
    }
    std::fill(block_.begin() + static_cast<std::ptrdiff_t>(used_), block_.begin() + 56, 0);
    for (int i = 0; i < 8; ++i)
      block_[56 + static_cast<std::size_t>(i)] = static_cast<std::uint8_t>(bits >> (56 - 8 * i));
    compress();
    std::string out;
    char digits[9];
    for (const std::uint32_t word : state_) {
      std::snprintf(digits, sizeof digits, "%08x", static_cast<unsigned>(word));
      out += digits;
    }
    return out;
  }

private:
  static constexpr std::uint32_t rotr(std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

  void compress() {
    static constexpr std::array<std::uint32_t, 64> k = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2};
    std::array<std::uint32_t, 64> w{};
    for (std::size_t i = 0; i < 16; ++i)
      w[i] = (std::uint32_t{block_[4 * i]} << 24) | (std::uint32_t{block_[4 * i + 1]} << 16) |
             (std::uint32_t{block_[4 * i + 2]} << 8) | std::uint32_t{block_[4 * i + 3]};
    for (std::size_t i = 16; i < 64; ++i) {
      const std::uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
      const std::uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    auto [a, b, c, d, e, f, g, h] = state_;
    for (std::size_t i = 0; i < 64; ++i) {
      const std::uint32_t t1 =
          h + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ (~e & g)) + k[i] + w[i];
      const std::uint32_t t2 =
          (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
      h = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
    }
    const std::array<std::uint32_t, 8> add = {a, b, c, d, e, f, g, h};
    for (std::size_t i = 0; i < 8; ++i)
      state_[i] += add[i];
  }

  std::array<std::uint32_t, 8> state_ = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                                         0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
  std::array<std::uint8_t, 64> block_{};
  std::size_t used_{};
  std::uint64_t length_{};
};

// Python's repr of a finite float: shortest round-trip digits, exponent form when the decimal
// point would lie more than 16 digits right of or 4 digits left of the first digit.
inline std::string python_float_repr(double value) {
  if (value == 0.0)
    return std::signbit(value) ? "-0.0" : "0.0";
  char buffer[64];
  const auto result =
      std::to_chars(buffer, buffer + sizeof buffer, value, std::chars_format::scientific);
  const std::string text(buffer, result.ptr);
  const std::size_t e = text.find('e');
  std::string sign, digits;
  for (std::size_t i = 0; i < e; ++i) {
    if (text[i] == '-')
      sign = "-";
    else if (text[i] != '.')
      digits += text[i];
  }
  const int exponent = std::stoi(text.substr(e + 1));
  const int point = exponent + 1; // Digits before the decimal point.
  const int count = static_cast<int>(digits.size());
  if (point <= -4 || point > 16) {
    std::string mantissa = digits.substr(0, 1);
    if (count > 1)
      mantissa += "." + digits.substr(1);
    char tail[16];
    std::snprintf(tail, sizeof tail, "e%c%02d", exponent < 0 ? '-' : '+', std::abs(exponent));
    return sign + mantissa + tail;
  }
  if (point <= 0)
    return sign + "0." + std::string(static_cast<std::size_t>(-point), '0') + digits;
  if (point >= count)
    return sign + digits + std::string(static_cast<std::size_t>(point - count), '0') + ".0";
  return sign + digits.substr(0, static_cast<std::size_t>(point)) + "." +
         digits.substr(static_cast<std::size_t>(point));
}

// Writes a string as json.dumps(..., ensure_ascii=False) does.
inline void append_python_string(std::string &out, std::string_view text) {
  out += '"';
  for (const char c : text) {
    switch (c) {
    case '"':
      out += "\\\"";
      break;
    case '\\':
      out += "\\\\";
      break;
    case '\n':
      out += "\\n";
      break;
    case '\r':
      out += "\\r";
      break;
    case '\t':
      out += "\\t";
      break;
    case '\b':
      out += "\\b";
      break;
    case '\f':
      out += "\\f";
      break;
    default:
      if (static_cast<unsigned char>(c) < 0x20) {
        char escape[8];
        std::snprintf(escape, sizeof escape, "\\u%04x", static_cast<unsigned>(c));
        out += escape;
      } else {
        out += c;
      }
    }
  }
  out += '"';
}

// Parses JSON text and writes its canonical form (see the file comment).
class Canonical {
public:
  explicit Canonical(std::string_view text) : text_(text) {}

  // Canonical text of the document; with `drop_key`, that top-level key is left out.
  std::string document(std::string_view drop_key) {
    space();
    std::string out = value(0, drop_key);
    space();
    if (at_ != text_.size())
      fail("trailing characters");
    return out;
  }

private:
  [[noreturn]] void fail(const char *what) const {
    throw ModelError("optical model: JSON syntax error at byte " + std::to_string(at_) + ": " +
                     what);
  }

  void space() {
    while (at_ < text_.size() &&
           (text_[at_] == ' ' || text_[at_] == '\n' || text_[at_] == '\r' || text_[at_] == '\t'))
      ++at_;
  }

  [[nodiscard]] char peek() const { return at_ < text_.size() ? text_[at_] : '\0'; }

  void expect(char c) {
    if (peek() != c)
      fail("unexpected character");
    ++at_;
  }

  std::string value(int depth, std::string_view drop_key = {}) {
    if (depth > 256)
      fail("nested too deeply");
    const char c = peek();
    if (c == '{')
      return object(depth, drop_key);
    if (c == '[')
      return array(depth);
    if (c == '"') {
      std::string out;
      append_python_string(out, string());
      return out;
    }
    if (c == '-' || (c >= '0' && c <= '9'))
      return number();
    for (const std::string_view word : {"true", "false", "null"})
      if (text_.substr(at_, word.size()) == word) {
        at_ += word.size();
        return std::string(word);
      }
    fail("unexpected value");
  }

  std::string object(int depth, std::string_view drop_key) {
    expect('{');
    // Later duplicates replace earlier ones, as in Python's json.loads.
    std::vector<std::pair<std::string, std::string>> members;
    space();
    if (peek() == '}') {
      ++at_;
      return "{}";
    }
    while (true) {
      space();
      if (peek() != '"')
        fail("expected a key");
      std::string key = string();
      space();
      expect(':');
      space();
      std::string item = value(depth + 1);
      const auto found = std::find_if(members.begin(), members.end(),
                                      [&](const auto &member) { return member.first == key; });
      if (found != members.end())
        found->second = std::move(item);
      else
        members.emplace_back(std::move(key), std::move(item));
      space();
      if (peek() == ',') {
        ++at_;
        continue;
      }
      expect('}');
      break;
    }
    // UTF-8 byte order is code point order.
    std::sort(members.begin(), members.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });
    std::string out = "{";
    bool first = true;
    for (const auto &[key, item] : members) {
      if (!drop_key.empty() && key == drop_key)
        continue;
      if (!first)
        out += ',';
      first = false;
      append_python_string(out, key);
      out += ':';
      out += item;
    }
    return out + "}";
  }

  std::string array(int depth) {
    expect('[');
    std::string out = "[";
    space();
    if (peek() == ']') {
      ++at_;
      return "[]";
    }
    while (true) {
      space();
      out += value(depth + 1);
      space();
      if (peek() == ',') {
        ++at_;
        out += ',';
        continue;
      }
      expect(']');
      break;
    }
    return out + "]";
  }

  std::string number() {
    const std::size_t start = at_;
    if (peek() == '-')
      ++at_;
    if (peek() == '0') {
      ++at_;
    } else if (peek() >= '1' && peek() <= '9') {
      while (peek() >= '0' && peek() <= '9')
        ++at_;
    } else {
      fail("invalid number");
    }
    bool integer = true;
    if (peek() == '.') {
      integer = false;
      ++at_;
      if (!(peek() >= '0' && peek() <= '9'))
        fail("invalid number");
      while (peek() >= '0' && peek() <= '9')
        ++at_;
    }
    if (peek() == 'e' || peek() == 'E') {
      integer = false;
      ++at_;
      if (peek() == '+' || peek() == '-')
        ++at_;
      if (!(peek() >= '0' && peek() <= '9'))
        fail("invalid number");
      while (peek() >= '0' && peek() <= '9')
        ++at_;
    }
    const std::string_view token = text_.substr(start, at_ - start);
    if (integer)
      return token == "-0" ? "0" : std::string(token);
    double value = 0.0;
    const auto result = std::from_chars(token.data(), token.data() + token.size(), value);
    if (result.ec != std::errc{} || !std::isfinite(value))
      fail("number out of range");
    return python_float_repr(value);
  }

  void append_utf8(std::string &out, std::uint32_t code) {
    if (code < 0x80) {
      out += static_cast<char>(code);
    } else if (code < 0x800) {
      out += static_cast<char>(0xC0 | (code >> 6));
      out += static_cast<char>(0x80 | (code & 0x3F));
    } else if (code < 0x10000) {
      out += static_cast<char>(0xE0 | (code >> 12));
      out += static_cast<char>(0x80 | ((code >> 6) & 0x3F));
      out += static_cast<char>(0x80 | (code & 0x3F));
    } else {
      out += static_cast<char>(0xF0 | (code >> 18));
      out += static_cast<char>(0x80 | ((code >> 12) & 0x3F));
      out += static_cast<char>(0x80 | ((code >> 6) & 0x3F));
      out += static_cast<char>(0x80 | (code & 0x3F));
    }
  }

  std::uint32_t hex4() {
    if (at_ + 4 > text_.size())
      fail("truncated \\u escape");
    std::uint32_t code = 0;
    const auto result = std::from_chars(text_.data() + at_, text_.data() + at_ + 4, code, 16);
    if (result.ec != std::errc{} || result.ptr != text_.data() + at_ + 4)
      fail("invalid \\u escape");
    at_ += 4;
    return code;
  }

  // The decoded (UTF-8) string at the cursor.
  std::string string() {
    expect('"');
    std::string out;
    while (at_ < text_.size()) {
      const char c = text_[at_++];
      if (c == '"')
        return out;
      if (static_cast<unsigned char>(c) < 0x20)
        fail("control character in string");
      if (c != '\\') {
        out += c;
        continue;
      }
      const char e = peek();
      ++at_;
      switch (e) {
      case '"':
      case '\\':
      case '/':
        out += e;
        break;
      case 'b':
        out += '\b';
        break;
      case 'f':
        out += '\f';
        break;
      case 'n':
        out += '\n';
        break;
      case 'r':
        out += '\r';
        break;
      case 't':
        out += '\t';
        break;
      case 'u': {
        std::uint32_t code = hex4();
        if (code >= 0xD800 && code < 0xDC00) {
          if (text_.substr(at_, 2) != "\\u")
            fail("unpaired surrogate");
          at_ += 2;
          const std::uint32_t low = hex4();
          if (!(low >= 0xDC00 && low < 0xE000))
            fail("unpaired surrogate");
          code = 0x10000 + ((code - 0xD800) << 10) + (low - 0xDC00);
        } else if (code >= 0xDC00 && code < 0xE000) {
          fail("unpaired surrogate");
        }
        append_utf8(out, code);
        break;
      }
      default:
        fail("invalid escape");
      }
    }
    fail("unterminated string");
  }

  std::string_view text_;
  std::size_t at_{};
};

} // namespace hash_detail

// SHA-256 (hex) of the canonical JSON of an optical-model document without its "sha256" entry.
[[nodiscard]] inline std::string model_hash(std::string_view text) {
  hash_detail::Sha256 sha;
  sha.update(hash_detail::Canonical(text).document("sha256"));
  return sha.hex();
}

} // namespace obdeect::ray_tracing
