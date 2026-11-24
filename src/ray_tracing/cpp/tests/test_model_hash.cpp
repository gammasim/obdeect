// Optical-model hash: SHA-256, Python float repr and canonical JSON match Python's
// compiler.model_hash (expected values computed with Python 3.14).
#include "check.hpp"

#include "obdeect/ray_tracing/model_hash.hpp"

#include <string>

using namespace obdeect::ray_tracing;

namespace {

std::string sha256(const std::string &text) {
  hash_detail::Sha256 sha;
  sha.update(text);
  return sha.hex();
}

void test_sha256_vectors() {
  CHECK(sha256("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
  CHECK(sha256("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
  // Two blocks: the padding does not fit after 56 bytes.
  CHECK(sha256("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq") ==
        "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
}

void test_python_float_repr() {
  const std::pair<double, const char *> cases[] = {{1.0, "1.0"},
                                                   {0.1, "0.1"},
                                                   {1e-5, "1e-05"},
                                                   {1e16, "1e+16"},
                                                   {1234567890123456.0, "1234567890123456.0"},
                                                   {0.0001, "0.0001"},
                                                   {1.5e300, "1.5e+300"},
                                                   {5e-324, "5e-324"},
                                                   {-2.5, "-2.5"},
                                                   {100.0, "100.0"},
                                                   {1e22, "1e+22"},
                                                   {0.30000000000000004, "0.30000000000000004"},
                                                   {-0.0, "-0.0"}};
  for (const auto &[value, text] : cases)
    CHECK(hash_detail::python_float_repr(value) == text);
}

void test_canonical_json_and_hash() {
  // Keys sorted, the last duplicate kept, "sha256" dropped, numbers and strings as Python writes.
  const std::string doc =
      R"({"z": [1, 2.50, -0, 1e-05, 1E16, 0.0001, 123456789012345678], )"
      R"("a": {"éé": "tab\t\"q\" \u0001 /\/", "b": true, "c": null, "d": -0.0}, )"
      R"("sha256": "x", "m": 0.1, "m": 3})";
  CHECK(model_hash(doc) == "914cdbdebd92fb3effffbf44399f18073e4dccc9c25b49e302d1ae31d141d0ae");
  bool rejected = false;
  try {
    (void)model_hash(R"({"a": 1,})");
  } catch (const ModelError &) {
    rejected = true;
  }
  CHECK(rejected);
}

// True if hashing `text` throws a ModelError.
bool rejects(const std::string &text) {
  try {
    (void)model_hash(text);
  } catch (const ModelError &) {
    return true;
  }
  return false;
}

// The parser and the canonical hash accept the same documents: escapes with surrogate pairs,
// the four JSON white-space characters, and the depth limit.
void test_same_documents_as_the_parser() {
  // \u escapes decode to UTF-8 (U+00E9 and U+1F600), as the raw UTF-8 text does.
  const std::string raw = "{\"a\": \"\xc3\xa9\xf0\x9f\x98\x80\"}";
  CHECK(model_hash(R"({"a": "é😀"})") == model_hash(raw));
  CHECK(rejects(R"({"a": "\ud800"})"));
  CHECK(rejects(R"({"a": "\udc00"})"));
  CHECK(rejects("{\"a\":\f1}"));
  CHECK(rejects("{\"a\":\v1}"));
  CHECK(model_hash("{\n\t\"a\" :\r1}") == model_hash(R"({"a": 1})"));
  std::string deep;
  for (int i = 0; i < json::max_depth + 1; ++i)
    deep += "[";
  deep += std::string(static_cast<std::size_t>(json::max_depth + 1), ']');
  CHECK(rejects("{\"a\": " + deep + "}"));
  std::string limit;
  for (int i = 0; i < json::max_depth - 1; ++i)
    limit += "[";
  limit += std::string(static_cast<std::size_t>(json::max_depth - 1), ']');
  CHECK(!rejects("{\"a\": " + limit + "}"));
}

void all_tests() {
  test_sha256_vectors();
  test_python_float_repr();
  test_canonical_json_and_hash();
  test_same_documents_as_the_parser();
}

} // namespace

int main() { return run_checks(all_tests); }
