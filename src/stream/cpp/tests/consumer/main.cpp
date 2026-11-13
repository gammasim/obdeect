// Uses the installed obdeect stream headers: a stream round trip and the stage ids.

#include "obdeect/stream/stages.hpp"
#include "obdeect/stream/stream.hpp"

#include <cstdio>
#include <exception>
#include <sstream>

static_assert(obdeect::stream::stage_atmosphere == 0x444C5546U);

static int run() {
  std::stringstream buffer;
  obdeect::stream::Writer writer(buffer, "consumer");
  writer.finish();

  obdeect::stream::Reader reader(buffer);
  obdeect::stream::Record record;
  while (reader.next(record)) {
  }
  return reader.ended() ? 0 : 1;
}

int main() {
  try {
    return run();
  } catch (const std::exception &error) {
    std::fprintf(stderr, "consumer: %s\n", error.what());
    return 1;
  }
}
