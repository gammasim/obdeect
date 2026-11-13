#include "check.hpp"

#include "obdeect/stream/cli.hpp"

#include <functional>
#include <string>
#include <vector>

using namespace obdeect::stream;

namespace {

bool usage_printed = false;

int run(std::vector<std::string> words, const std::function<int()> &body) {
  std::vector<std::string> storage{"tool"};
  storage.insert(storage.end(), words.begin(), words.end());
  std::vector<char *> argv;
  for (auto &word : storage)
    argv.push_back(word.data());
  usage_printed = false;
  return cli::run_tool(
      "tool", "1.2.3", static_cast<int>(argv.size()), argv.data(), [] { usage_printed = true; },
      body);
}

} // namespace

int main() {
  return run_checks([] {
    // Option cursor: value() fails with "<option>: missing value".
    std::vector<std::string> storage{"tool", "-n", "5", "-x"};
    std::vector<char *> argv;
    for (auto &word : storage)
      argv.push_back(word.data());
    cli::Args args(static_cast<int>(argv.size()), argv.data());
    std::string arg;
    CHECK(args.next(arg) && arg == "-n");
    CHECK(args.value(arg) == "5");
    CHECK(args.next(arg) && arg == "-x");
    bool missing = false;
    try {
      (void)args.value(arg);
    } catch (const cli::UsageError &error) {
      missing = std::string(error.what()) == "-x: missing value";
    }
    CHECK(missing);
    CHECK(!args.next(arg));

    CHECK(cli::number("1.5", "-a") == 1.5);
    CHECK(cli::integer("-7", "-n") == -7);
    CHECK(cli::integer("7", "-n", 1, 7) == 7);
    const auto fails = [](const std::function<void()> &call) {
      try {
        call();
      } catch (const cli::UsageError &) {
        return true;
      }
      return false;
    };
    CHECK(fails([] { (void)cli::number("1.5x", "-a"); }));
    CHECK(fails([] { (void)cli::number("nan", "-a"); }));
    CHECK(fails([] { (void)cli::integer("8", "-n", 1, 7); }));
    CHECK(fails([] { (void)cli::integer("7.5", "-n"); }));
    try {
      (void)cli::number("abc", "--qe");
      CHECK(false);
    } catch (const cli::UsageError &error) {
      CHECK(std::string(error.what()) == "--qe: not a finite number: abc");
    }

    // Exit status: 0 ok, 2 usage error (usage printed), 1 anything else (no usage).
    CHECK(run({}, [] { return 0; }) == 0 && !usage_printed);
    CHECK(run({"--version"}, [] { return 5; }) == 0 && !usage_printed);
    CHECK(run({}, []() -> int { throw cli::UsageError("bad"); }) == 2 && usage_printed);
    CHECK(run({}, []() -> int { throw std::runtime_error("io"); }) == 1 && !usage_printed);
    CHECK(run({}, []() -> int { throw std::invalid_argument("data"); }) == 1 && !usage_printed);

    CHECK(cli::version_line("tool", "1.2.3").rfind("tool 1.2.3 (photon-stream format ", 0) == 0);
  });
}
