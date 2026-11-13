/* Uses the installed C library through obdeect::stream_c. */
#include "obdeect/stream/obdeect_stream.h"

int main(void) {
  return obdeect_stream_abi_version() == OBDEECT_STREAM_C_ABI_VERSION && obdeect_stream_format_major() >= 1 ? 0 : 1;
}
