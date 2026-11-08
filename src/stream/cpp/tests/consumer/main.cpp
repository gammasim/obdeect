// Uses the installed gstrickd headers: a stream round trip and the stage ids.

#include "gstrickd/stages.hpp"
#include "gstrickd/stream.hpp"

#include <sstream>

static_assert(gstrickd::stage_atmosphere == 0x444C5546U);

int main() {
  std::stringstream buffer;
  gstrickd::Writer writer(buffer, "consumer");
  writer.finish();

  gstrickd::Reader reader(buffer);
  gstrickd::Record record;
  while (reader.next(record)) {
  }
  return reader.ended() ? 0 : 1;
}
