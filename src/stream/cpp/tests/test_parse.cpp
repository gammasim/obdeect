#include "check.hpp"

#include "obdeect/stream/parse.hpp"

#include <clocale>
#include <cstdint>

using namespace obdeect::stream::parse;

int main() {
  return run_checks([] {
    // A decimal-comma locale must not change the result (best effort: the locale may be absent).
    (void)std::setlocale(LC_ALL, "de_DE.UTF-8");

    double v = 0.0;
    CHECK(finite_number("1.5", v) && v == 1.5);
    CHECK(finite_number("-4e2", v) && v == -400.0);
    CHECK(finite_number("+2", v) && v == 2.0);
    CHECK(!finite_number("", v));
    CHECK(!finite_number("+", v));
    CHECK(!finite_number("4e2x", v));
    CHECK(!finite_number(" 1", v));
    CHECK(!finite_number("1 ", v));
    CHECK(!finite_number("nan", v));
    CHECK(!finite_number("inf", v));
    CHECK(!finite_number("1e999", v));
    CHECK(!finite_number("++1", v));

    long n = 0;
    CHECK(integer("42", n) && n == 42);
    CHECK(integer("-7", n) && n == -7);
    CHECK(integer("+7", n) && n == 7);
    CHECK(!integer("", n));
    CHECK(!integer("1.5", n));
    CHECK(!integer("12x", n));
    std::uint8_t small = 0;
    CHECK(integer("255", small) && small == 255);
    CHECK(!integer("256", small));
    CHECK(!integer("-1", small));
  });
}
