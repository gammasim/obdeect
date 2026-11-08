/* Uses the installed C library through gstrickd::gstrickd_c. */
#include "gstrickd/gstrickd.h"

int main(void) {
  return gstrickd_abi_version() == GSTRICKD_C_ABI_VERSION && gstrickd_format_major() >= 1 ? 0 : 1;
}
