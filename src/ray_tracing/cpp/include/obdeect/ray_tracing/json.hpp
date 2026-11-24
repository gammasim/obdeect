#pragma once

#include "obdeect/stream/parse.hpp"
#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <variant>
#include <vector>
#include <version>

namespace obdeect::ray_tracing::json {

// Nesting limit of optical-model documents: the parser and the canonical hash reject deeper input.
inline constexpr int max_depth = 64;

// JSON white space: exactly these four characters (std::isspace also takes \f and \v).
[[nodiscard]] constexpr bool is_space(char c) {
  return c == ' ' || c == '\t' || c == '\n' || c == '\r';
}

// Appends the UTF-8 encoding of a code point.
inline void append_utf8(std::string &out, std::uint32_t code) {
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

// Minimal strict JSON document model for optical-model files.
struct Value {
  enum class Kind { null, boolean, number, string, array, object };
  using Array = std::vector<Value>;
  // Parallel vectors: std::vector accepts the incomplete Value, std::pair does not.
  struct Object {
    std::vector<std::string> names;
    std::vector<Value> values;
  };

  std::variant<std::monostate, bool, double, std::string, Array, Object> data;

  [[nodiscard]] Kind kind() const { return static_cast<Kind>(data.index()); }
  [[nodiscard]] const Value *find(std::string_view key) const {
    if (kind() != Kind::object)
      return nullptr;
    const auto &object = std::get<Object>(data);
    for (std::size_t i = 0; i < object.names.size(); ++i)
      if (object.names[i] == key)
        return &object.values[i];
    return nullptr;
  }
};

class Parser {
public:
  // Throws std::runtime_error with the byte offset of the first syntax error.
  [[nodiscard]] static Value parse(std::string_view text) {
    Parser parser{text};
    parser.space();
    Value value = parser.value(0);
    parser.space();
    if (parser.at_ != text.size())
      parser.fail("trailing characters");
    return value;
  }

private:
  explicit Parser(std::string_view text) : text_(text) {}

  [[noreturn]] void fail(const char *what) const {
    throw std::runtime_error("JSON syntax error at byte " + std::to_string(at_) + ": " + what);
  }
  void space() {
    while (at_ < text_.size() && is_space(text_[at_]))
      ++at_;
  }
  bool take(char c) {
    if (at_ < text_.size() && text_[at_] == c) {
      ++at_;
      return true;
    }
    return false;
  }
  void expect(char c) {
    if (!take(c))
      fail("unexpected character");
  }

  Value value(int depth) {
    if (depth > max_depth)
      fail("nesting too deep");
    if (at_ == text_.size())
      fail("unexpected end of input");
    switch (text_[at_]) {
    case '{':
      return object(depth);
    case '[':
      return array(depth);
    case '"':
      return Value{string()};
    case 't':
      return literal("true", Value{true});
    case 'f':
      return literal("false", Value{false});
    case 'n':
      return literal("null", Value{});
    default:
      return Value{number()};
    }
  }

  Value literal(std::string_view word, Value result) {
    if (text_.substr(at_, word.size()) != word)
      fail("invalid literal");
    at_ += word.size();
    return result;
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

  // The decoded UTF-8 string at the cursor. Surrogate pairs decode to one code point; a lone
  // surrogate is rejected (as the canonical form of the hash does).
  std::string string() {
    expect('"');
    std::string result;
    while (at_ < text_.size()) {
      const char c = text_[at_++];
      if (c == '"')
        return result;
      if (static_cast<unsigned char>(c) < 0x20)
        fail("control character in string");
      if (c != '\\') {
        result.push_back(c);
        continue;
      }
      if (at_ == text_.size())
        break;
      const char e = text_[at_++];
      const std::string_view plain = "\"\\/", escaped = "bfnrt", decoded = "\b\f\n\r\t";
      if (plain.find(e) != std::string_view::npos) {
        result.push_back(e);
      } else if (const auto i = escaped.find(e); i != std::string_view::npos) {
        result.push_back(decoded[i]);
      } else if (e == 'u') {
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
        append_utf8(result, code);
      } else {
        fail("unsupported escape sequence");
      }
    }
    fail("unterminated string");
  }

  double number() {
    const std::size_t start = at_;
    take('-');
    if (!take('0')) {
      if (at_ == text_.size() || !std::isdigit(static_cast<unsigned char>(text_[at_])))
        fail("invalid number");
      digits();
    }
    if (take('.') && !digits())
      fail("invalid fraction");
    if (take('e') || take('E')) {
      if (!take('+'))
        take('-');
      if (!digits())
        fail("invalid exponent");
    }
    // The text was checked against the JSON grammar above. Parse it independent of
    // the process locale (strtod would read "1.5" as 1 under a decimal-comma locale).
    const std::string_view literal = text_.substr(start, at_ - start);
    double result = 0.0;
    if (!stream::parse::finite_number(literal, result))
      fail("number out of range");
    return result;
  }
  bool digits() {
    const std::size_t start = at_;
    while (at_ < text_.size() && std::isdigit(static_cast<unsigned char>(text_[at_])))
      ++at_;
    return at_ > start;
  }

  Value array(int depth) {
    expect('[');
    Value::Array items;
    space();
    if (!take(']')) {
      do {
        space();
        items.push_back(value(depth + 1));
        space();
      } while (take(','));
      expect(']');
    }
    return Value{std::move(items)};
  }

  Value object(int depth) {
    expect('{');
    Value::Object members;
    space();
    if (!take('}')) {
      do {
        space();
        std::string key = string();
        space();
        expect(':');
        space();
        Value item = value(depth + 1);
        // Later duplicates replace earlier ones, as in Python's json.loads.
        const auto found = std::find(members.names.begin(), members.names.end(), key);
        if (found != members.names.end()) {
          members.values[static_cast<std::size_t>(found - members.names.begin())] = std::move(item);
        } else {
          members.names.push_back(std::move(key));
          members.values.push_back(std::move(item));
        }
        space();
      } while (take(','));
      expect('}');
    }
    return Value{std::move(members)};
  }

  std::string_view text_;
  std::size_t at_{};
};

} // namespace obdeect::ray_tracing::json
