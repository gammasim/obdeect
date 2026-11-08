/* C11 test of the C API: writes a stream (by path and by file descriptor), reads it back,
 * checks errors, and reads the conformance streams. Usage: test_capi <output.ps> <conformance dir>.
 * The output file is compared with the C++ reader by test_capi_cpp. */
#define _POSIX_C_SOURCE 200809L

#include "gstrickd/gstrickd.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fcntl.h>
#include <unistd.h>

static int failures = 0;

#define CHECK(condition)                                                                           \
  do {                                                                                             \
    if (!(condition)) {                                                                            \
      fprintf(stderr, "%s:%d: CHECK failed: %s (%s)\n", __FILE__, __LINE__, #condition,            \
              gstrickd_last_error());                                                              \
      ++failures;                                                                                  \
    }                                                                                              \
  } while (0)

enum { ROWS = 5 };
static const int mask_columns[] = {GSTRICKD_COLUMN_X,         GSTRICKD_COLUMN_Y, GSTRICKD_COLUMN_CX,
                                   GSTRICKD_COLUMN_CY,        GSTRICKD_COLUMN_T, GSTRICKD_COLUMN_W,
                                   GSTRICKD_COLUMN_WAVELENGTH};

static float value_of(int column, int row) { return (float)column * 100.0F + (float)row * 0.5F; }

static void write_stream(gstrickd_writer *writer) {
  float data[GSTRICKD_N_COLUMNS][ROWS];
  const float *columns[GSTRICKD_N_COLUMNS] = {0};
  const double positions[] = {0.0, 12.5, -3.0};
  const uint8_t block[] = {1, 2, 3, 250};
  size_t i;
  int row;

  CHECK(gstrickd_writer_set_int(writer, "run_number", 42) == GSTRICKD_OK);
  CHECK(gstrickd_writer_set_float(writer, "energy_min_gev", 0.25) == GSTRICKD_OK);
  CHECK(gstrickd_writer_set_string(writer, "site", "La Palma") == GSTRICKD_OK);
  CHECK(gstrickd_writer_set_array(writer, "telescope_x_m", positions, 3) == GSTRICKD_OK);
  CHECK(gstrickd_writer_write_record(writer, GSTRICKD_RECORD_RUN) == GSTRICKD_OK);
  /* The pending metadata was cleared: this record is empty. */
  CHECK(gstrickd_writer_write_record(writer, GSTRICKD_RECORD_EVENT) == GSTRICKD_OK);
  CHECK(gstrickd_writer_set_int(writer, "telescope_id", 3) == GSTRICKD_OK);
  CHECK(gstrickd_writer_set_int(writer, "frame", GSTRICKD_FRAME_GROUND) == GSTRICKD_OK);
  CHECK(gstrickd_writer_write_record(writer, GSTRICKD_RECORD_TELESCOPE) == GSTRICKD_OK);

  for (i = 0; i < sizeof mask_columns / sizeof mask_columns[0]; ++i) {
    for (row = 0; row < ROWS; ++row) {
      data[mask_columns[i]][row] = value_of(mask_columns[i], row);
    }
    columns[mask_columns[i]] = data[mask_columns[i]];
  }
  CHECK(gstrickd_writer_write_photons(writer, 3, GSTRICKD_FRAME_GROUND, 100, ROWS, columns) ==
        GSTRICKD_OK);
  CHECK(gstrickd_writer_write_photons(writer, 3, GSTRICKD_FRAME_CAMERA, 105, 0, NULL) ==
        GSTRICKD_OK);
  CHECK(gstrickd_writer_write_blob(writer, 1204, 7, block, sizeof block) == GSTRICKD_OK);
  CHECK(gstrickd_writer_set_float(writer, "photons", 5.0) == GSTRICKD_OK);
  CHECK(gstrickd_writer_write_record(writer, GSTRICKD_RECORD_TELESCOPE_END) == GSTRICKD_OK);
  CHECK(gstrickd_writer_write_record(writer, GSTRICKD_RECORD_RUN_END) == GSTRICKD_OK);
}

static void expect_next(gstrickd_reader *reader, int type) {
  int got = 0;
  CHECK(gstrickd_reader_next(reader, &got) == GSTRICKD_OK);
  CHECK(got == type);
}

static void check_stream(gstrickd_reader *reader) {
  int type = 0;
  int value_type = 0;
  int64_t integer = 0;
  double real = 0.0;
  const char *text = NULL;
  size_t length = 0, n = 0, i;
  const double *array = NULL;
  uint32_t telescope = 0, frame = 0, rows = 0, mask = 0, eventio_type = 0, eventio_id = 0;
  uint64_t first = 0;
  const uint8_t *bytes = NULL;
  int row;

  expect_next(reader, GSTRICKD_RECORD_STREAM);
  CHECK(gstrickd_reader_meta_int(reader, "format_major", &integer) == GSTRICKD_OK);
  CHECK(integer == gstrickd_format_major());
  CHECK(gstrickd_reader_meta_string(reader, "producer", &text, &length) == GSTRICKD_OK);
  CHECK(strcmp(text, "test_capi") == 0 && length == 9);

  expect_next(reader, GSTRICKD_RECORD_RUN);
  CHECK(gstrickd_reader_n_keys(reader, &n) == GSTRICKD_OK && n == 4);
  CHECK(strcmp(gstrickd_reader_key(reader, 0), "energy_min_gev") == 0); /* keys are sorted */
  CHECK(gstrickd_reader_key(reader, 4) == NULL);
  CHECK(gstrickd_reader_meta_type(reader, "run_number", &value_type) == GSTRICKD_OK);
  CHECK(value_type == GSTRICKD_VALUE_INT);
  CHECK(gstrickd_reader_meta_type(reader, "telescope_x_m", &value_type) == GSTRICKD_OK);
  CHECK(value_type == GSTRICKD_VALUE_ARRAY);
  CHECK(gstrickd_reader_meta_int(reader, "run_number", &integer) == GSTRICKD_OK && integer == 42);
  CHECK(gstrickd_reader_meta_float(reader, "energy_min_gev", &real) == GSTRICKD_OK);
  CHECK(real == 0.25);
  CHECK(gstrickd_reader_meta_float(reader, "run_number", &real) == GSTRICKD_OK && real == 42.0);
  CHECK(gstrickd_reader_meta_string(reader, "site", &text, NULL) == GSTRICKD_OK);
  CHECK(strcmp(text, "La Palma") == 0);
  CHECK(gstrickd_reader_meta_array(reader, "telescope_x_m", &array, &n) == GSTRICKD_OK);
  CHECK(n == 3 && array[1] == 12.5 && array[2] == -3.0);
  CHECK(gstrickd_reader_meta_int(reader, "site", &integer) == GSTRICKD_ERROR_TYPE);
  CHECK(gstrickd_reader_meta_int(reader, "missing", &integer) == GSTRICKD_ERROR_NOT_FOUND);
  CHECK(gstrickd_reader_photons(reader, NULL, NULL, NULL, NULL, NULL) == GSTRICKD_ERROR_STATE);
  CHECK(gstrickd_reader_column(reader, GSTRICKD_COLUMN_X) == NULL);

  expect_next(reader, GSTRICKD_RECORD_EVENT);
  CHECK(gstrickd_reader_n_keys(reader, &n) == GSTRICKD_OK && n == 0);
  expect_next(reader, GSTRICKD_RECORD_TELESCOPE);
  CHECK(gstrickd_reader_meta_int(reader, "telescope_id", &integer) == GSTRICKD_OK && integer == 3);

  expect_next(reader, GSTRICKD_RECORD_PHOTONS);
  CHECK(gstrickd_reader_photons(reader, &telescope, &frame, &first, &rows, &mask) == GSTRICKD_OK);
  CHECK(telescope == 3 && frame == GSTRICKD_FRAME_GROUND && first == 100 && rows == ROWS);
  for (i = 0; i < GSTRICKD_N_COLUMNS; ++i) {
    const float *column = gstrickd_reader_column(reader, (int)i);
    const int present = (mask >> i) & 1U;
    CHECK((column != NULL) == present);
    for (row = 0; present && row < ROWS; ++row) {
      CHECK(column[row] == value_of((int)i, row));
    }
  }
  CHECK(mask ==
        ((1U << GSTRICKD_COLUMN_X) | (1U << GSTRICKD_COLUMN_Y) | (1U << GSTRICKD_COLUMN_CX) |
         (1U << GSTRICKD_COLUMN_CY) | (1U << GSTRICKD_COLUMN_T) | (1U << GSTRICKD_COLUMN_W) |
         (1U << GSTRICKD_COLUMN_WAVELENGTH)));
  CHECK(gstrickd_reader_n_keys(reader, &n) == GSTRICKD_ERROR_STATE);

  expect_next(reader, GSTRICKD_RECORD_PHOTONS);
  CHECK(gstrickd_reader_photons(reader, &telescope, &frame, &first, &rows, &mask) == GSTRICKD_OK);
  CHECK(frame == GSTRICKD_FRAME_CAMERA && first == 105 && rows == 0 && mask == 0);
  CHECK(gstrickd_reader_column(reader, GSTRICKD_COLUMN_X) == NULL);

  expect_next(reader, GSTRICKD_RECORD_BLOB);
  CHECK(gstrickd_reader_blob(reader, &eventio_type, &eventio_id, &bytes, &n) == GSTRICKD_OK);
  CHECK(eventio_type == 1204 && eventio_id == 7 && n == 4 && bytes[3] == 250);

  expect_next(reader, GSTRICKD_RECORD_TELESCOPE_END);
  CHECK(gstrickd_reader_meta_float(reader, "photons", &real) == GSTRICKD_OK && real == 5.0);
  expect_next(reader, GSTRICKD_RECORD_RUN_END);
  CHECK(gstrickd_reader_next(reader, &type) == GSTRICKD_END && type == 0);
}

static void test_arguments(const char *directory) {
  gstrickd_writer *writer = NULL;
  gstrickd_reader *reader = NULL;
  char path[1024];
  float column[2] = {0.0F, 0.0F};
  const float *columns[GSTRICKD_N_COLUMNS] = {0};
  float *many = NULL;
  int type = 0;

  CHECK(gstrickd_abi_version() == GSTRICKD_C_ABI_VERSION);
  CHECK(gstrickd_format_major() == 1);
  CHECK(strlen(gstrickd_version()) > 0);
  CHECK(strcmp(gstrickd_record_type_name(GSTRICKD_RECORD_PHOTONS), "photons") == 0);
  CHECK(strcmp(gstrickd_column_name(GSTRICKD_COLUMN_WAVELENGTH), "wavelength") == 0);
  CHECK(gstrickd_column_name(11) == NULL && gstrickd_record_type_name(77) == NULL);

  CHECK(gstrickd_reader_open("/nonexistent/x.ps", &reader) == GSTRICKD_ERROR_IO);
  CHECK(reader == NULL && strlen(gstrickd_last_error()) > 0);
  CHECK(gstrickd_reader_open(NULL, &reader) == GSTRICKD_ERROR_ARGUMENT);
  CHECK(gstrickd_writer_open("/nonexistent/x.ps", "t", &writer) == GSTRICKD_ERROR_IO);
  gstrickd_reader_close(NULL);
  CHECK(gstrickd_writer_close(NULL) == GSTRICKD_OK);

  snprintf(path, sizeof path, "%.500s/arguments.ps", directory);
  CHECK(gstrickd_writer_open(path, "t", &writer) == GSTRICKD_OK);
  CHECK(gstrickd_writer_write_record(writer, GSTRICKD_RECORD_PHOTONS) == GSTRICKD_ERROR_ARGUMENT);
  CHECK(gstrickd_writer_write_record(writer, GSTRICKD_RECORD_STREAM) == GSTRICKD_ERROR_ARGUMENT);
  CHECK(gstrickd_writer_write_photons(writer, 1, 4, 0, 2, columns) == GSTRICKD_ERROR_ARGUMENT);
  /* Required columns missing. */
  columns[GSTRICKD_COLUMN_X] = column;
  CHECK(gstrickd_writer_write_photons(writer, 1, 1, 0, 2, columns) == GSTRICKD_ERROR_ARGUMENT);
  CHECK(gstrickd_writer_write_photons(writer, 1, 1, 0, 2, NULL) == GSTRICKD_ERROR_ARGUMENT);
  /* More than 2^16 rows. */
  many = calloc((size_t)GSTRICKD_MAX_CHUNK_ROWS + 1, sizeof(float));
  CHECK(many != NULL);
  if (many != NULL) {
    int c;
    for (c = 0; c < GSTRICKD_N_COLUMNS; ++c) {
      columns[c] = many;
    }
    CHECK(gstrickd_writer_write_photons(writer, 1, 1, 0, GSTRICKD_MAX_CHUNK_ROWS + 1, columns) ==
          GSTRICKD_ERROR_ARGUMENT);
    CHECK(gstrickd_writer_write_photons(writer, 1, 1, 0, GSTRICKD_MAX_CHUNK_ROWS, columns) ==
          GSTRICKD_OK);
    free(many);
  }
  /* Empty array and empty blob with NULL pointers. */
  CHECK(gstrickd_writer_set_array(writer, "empty", NULL, 0) == GSTRICKD_OK);
  CHECK(gstrickd_writer_write_record(writer, GSTRICKD_RECORD_RUN) == GSTRICKD_OK);
  CHECK(gstrickd_writer_write_blob(writer, 1, 2, NULL, 0) == GSTRICKD_OK);
  CHECK(gstrickd_writer_close(writer) == GSTRICKD_OK);

  CHECK(gstrickd_reader_open(path, &reader) == GSTRICKD_OK);
  expect_next(reader, GSTRICKD_RECORD_STREAM);
  expect_next(reader, GSTRICKD_RECORD_PHOTONS);
  expect_next(reader, GSTRICKD_RECORD_RUN);
  {
    const double *values = NULL;
    size_t n = 7;
    CHECK(gstrickd_reader_meta_array(reader, "empty", &values, &n) == GSTRICKD_OK && n == 0);
  }
  expect_next(reader, GSTRICKD_RECORD_BLOB);
  CHECK(gstrickd_reader_next(reader, &type) == GSTRICKD_END);
  gstrickd_reader_close(reader);
}

/* Every stream in the conformance set: the valid one ends cleanly after 11 records, the
 * invalid ones fail with the given code. */
static int count_records(const char *path, int *last) {
  gstrickd_reader *reader = NULL;
  int type = 0, code, count = 0;
  if (gstrickd_reader_open(path, &reader) != GSTRICKD_OK) {
    *last = GSTRICKD_ERROR_IO;
    return -1;
  }
  while ((code = gstrickd_reader_next(reader, &type)) == GSTRICKD_OK) {
    ++count;
  }
  *last = code;
  gstrickd_reader_close(reader);
  return count;
}

static void test_conformance(const char *directory) {
  const char *invalid[] = {"bad_magic.ps", "unknown_major.ps", "truncated_chunk.ps",
                           "too_many_rows.ps"};
  const int codes[] = {GSTRICKD_ERROR_FORMAT, GSTRICKD_ERROR_VERSION, GSTRICKD_ERROR_FORMAT,
                       GSTRICKD_ERROR_FORMAT};
  char path[1024];
  int last = 0, i;

  snprintf(path, sizeof path, "%.500s/valid.ps", directory);
  CHECK(count_records(path, &last) == 11 && last == GSTRICKD_END);
  for (i = 0; i < 4; ++i) {
    snprintf(path, sizeof path, "%.500s/%.64s", directory, invalid[i]);
    count_records(path, &last);
    CHECK(last == codes[i]);
  }
}

int main(int argc, char **argv) {
  gstrickd_writer *writer = NULL;
  gstrickd_reader *reader = NULL;
  char directory[512];
  const char *slash;

  if (argc != 3) {
    fprintf(stderr, "usage: test_capi <output.ps> <conformance dir>\n");
    return 2;
  }
  slash = strrchr(argv[1], '/');
  snprintf(directory, sizeof directory, "%.*s", slash ? (int)(slash - argv[1]) : 1,
           slash ? argv[1] : ".");

  CHECK(gstrickd_writer_open(argv[1], "test_capi", &writer) == GSTRICKD_OK);
  if (writer != NULL) {
    write_stream(writer);
    CHECK(gstrickd_writer_close(writer) == GSTRICKD_OK);
  }
  CHECK(gstrickd_reader_open(argv[1], &reader) == GSTRICKD_OK);
  if (reader != NULL) {
    check_stream(reader);
    gstrickd_reader_close(reader);
  }

  { /* The same through file descriptors. */
    char path[1100];
    int fd;
    snprintf(path, sizeof path, "%.1000s.fd", argv[1]);
    fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    CHECK(fd >= 0);
    CHECK(gstrickd_writer_open_fd(fd, "test_capi", &writer) == GSTRICKD_OK);
    write_stream(writer);
    CHECK(gstrickd_writer_close(writer) == GSTRICKD_OK);
    CHECK(close(fd) == 0);
    fd = open(path, O_RDONLY);
    CHECK(fd >= 0);
    CHECK(gstrickd_reader_open_fd(fd, &reader) == GSTRICKD_OK);
    check_stream(reader);
    gstrickd_reader_close(reader);
    close(fd);
    CHECK(gstrickd_reader_open_fd(-1, &reader) == GSTRICKD_ERROR_ARGUMENT);
  }

  test_arguments(directory);
  test_conformance(argv[2]);

  if (failures != 0) {
    fprintf(stderr, "%d check(s) failed\n", failures);
    return 1;
  }
  return 0;
}
