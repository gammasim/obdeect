#pragma once

// Locale-independent number parsing. The whole token must be the number; a leading '+' is
// accepted. Callers add their own context to the error message.

#include <cctype>
#include <charconv>
#include <cmath>
#include <locale>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>

namespace obdeect::stream::parse {

// A finite real number filling the whole token.
[[nodiscard]] inline bool finite_number(std::string_view token, double &value) {
  if (token.size() > 1 && token.front() == '+' && token[1] != '-' && token[1] != '+')
    token.remove_prefix(1);
  if (token.empty())
    return false;
#if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
  double parsed = 0.0;
  const auto result = std::from_chars(token.data(), token.data() + token.size(), parsed);
  if (result.ec != std::errc{} || result.ptr != token.data() + token.size())
    return false;
#else
  if (std::isspace(static_cast<unsigned char>(token.front())))
    return false;
  std::istringstream stream{std::string(token)};
  stream.imbue(std::locale::classic());
  double parsed = 0.0;
  if (!(stream >> parsed) || stream.peek() != std::char_traits<char>::eof())
    return false;
#endif
  if (!std::isfinite(parsed))
    return false;
  value = parsed;
  return true;
}

// A decimal integer (optional sign) filling the whole token and fitting in T.
template <class T> [[nodiscard]] inline bool integer(std::string_view token, T &value) {
  if (token.size() > 1 && token.front() == '+' && token[1] != '-' && token[1] != '+')
    token.remove_prefix(1);
  if (token.empty())
    return false;
  T parsed{};
  const auto result = std::from_chars(token.data(), token.data() + token.size(), parsed, 10);
  if (result.ec != std::errc{} || result.ptr != token.data() + token.size())
    return false;
  value = parsed;
  return true;
}

} // namespace obdeect::stream::parse
