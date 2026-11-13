/* C11 test of the C API: writes a stream (by path and by file descriptor), reads it back,
 * checks errors, and reads the conformance streams. Usage: test_capi <output.ps> <conformance dir>.
 * The output file is compared with the C++ reader by test_capi_cpp. */
#define _POSIX_C_SOURCE 200809L

#include "obdeect/stream/obdeect_stream.h"

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
              obdeect_stream_last_error());                                                              \
      ++failures;                                                                                  \
    }                                                                                              \
  } while (0)

enum { ROWS = 5 };
static const int mask_columns[] = {OBDEECT_STREAM_COLUMN_X,         OBDEECT_STREAM_COLUMN_Y, OBDEECT_STREAM_COLUMN_CX,
                                   OBDEECT_STREAM_COLUMN_CY,        OBDEECT_STREAM_COLUMN_T, OBDEECT_STREAM_COLUMN_W,
                                   OBDEECT_STREAM_COLUMN_WAVELENGTH};

static float value_of(int column, int row) { return (float)column * 100.0F + (float)row * 0.5F; }

static void write_stream(obdeect_stream_writer *writer) {
  float data[OBDEECT_STREAM_N_COLUMNS][ROWS];
  const float *columns[OBDEECT_STREAM_N_COLUMNS] = {0};
  const double positions[] = {0.0, 12.5, -3.0};
  const uint8_t block[] = {1, 2, 3, 250};
  size_t i;
  int row;

  CHECK(obdeect_stream_writer_set_int(writer, "run_number", 42) == OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_set_float(writer, "energy_min_gev", 0.25) == OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_set_string(writer, "site", "La Palma") == OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_set_array(writer, "telescope_x_m", positions, 3) == OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_write_record(writer, OBDEECT_STREAM_RECORD_RUN) == OBDEECT_STREAM_OK);
  /* The pending metadata was cleared: this record is empty. */
  CHECK(obdeect_stream_writer_write_record(writer, OBDEECT_STREAM_RECORD_EVENT) == OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_set_int(writer, "telescope_id", 3) == OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_set_int(writer, "frame", OBDEECT_STREAM_FRAME_GROUND) == OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_write_record(writer, OBDEECT_STREAM_RECORD_TELESCOPE) == OBDEECT_STREAM_OK);

  for (i = 0; i < sizeof mask_columns / sizeof mask_columns[0]; ++i) {
    for (row = 0; row < ROWS; ++row) {
      data[mask_columns[i]][row] = value_of(mask_columns[i], row);
    }
    columns[mask_columns[i]] = data[mask_columns[i]];
  }
  CHECK(obdeect_stream_writer_write_photons(writer, 3, OBDEECT_STREAM_FRAME_GROUND, 100, ROWS, columns) ==
        OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_write_photons(writer, 3, OBDEECT_STREAM_FRAME_CAMERA, 105, 0, NULL) ==
        OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_write_blob(writer, 1204, 7, block, sizeof block) == OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_set_float(writer, "photons", 5.0) == OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_write_record(writer, OBDEECT_STREAM_RECORD_TELESCOPE_END) == OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_write_record(writer, OBDEECT_STREAM_RECORD_RUN_END) == OBDEECT_STREAM_OK);
}

static void expect_next(obdeect_stream_reader *reader, int type) {
  int got = 0;
  CHECK(obdeect_stream_reader_next(reader, &got) == OBDEECT_STREAM_OK);
  CHECK(got == type);
}

static void check_stream(obdeect_stream_reader *reader) {
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

  expect_next(reader, OBDEECT_STREAM_RECORD_STREAM);
  CHECK(obdeect_stream_reader_meta_int(reader, "format_major", &integer) == OBDEECT_STREAM_OK);
  CHECK(integer == obdeect_stream_format_major());
  CHECK(obdeect_stream_reader_meta_string(reader, "producer", &text, &length) == OBDEECT_STREAM_OK);
  CHECK(strcmp(text, "test_capi") == 0 && length == 9);

  expect_next(reader, OBDEECT_STREAM_RECORD_RUN);
  CHECK(obdeect_stream_reader_n_keys(reader, &n) == OBDEECT_STREAM_OK && n == 4);
  CHECK(strcmp(obdeect_stream_reader_key(reader, 0), "energy_min_gev") == 0); /* keys are sorted */
  CHECK(obdeect_stream_reader_key(reader, 4) == NULL);
  CHECK(obdeect_stream_reader_meta_type(reader, "run_number", &value_type) == OBDEECT_STREAM_OK);
  CHECK(value_type == OBDEECT_STREAM_VALUE_INT);
  CHECK(obdeect_stream_reader_meta_type(reader, "telescope_x_m", &value_type) == OBDEECT_STREAM_OK);
  CHECK(value_type == OBDEECT_STREAM_VALUE_ARRAY);
  CHECK(obdeect_stream_reader_meta_int(reader, "run_number", &integer) == OBDEECT_STREAM_OK && integer == 42);
  CHECK(obdeect_stream_reader_meta_float(reader, "energy_min_gev", &real) == OBDEECT_STREAM_OK);
  CHECK(real == 0.25);
  CHECK(obdeect_stream_reader_meta_float(reader, "run_number", &real) == OBDEECT_STREAM_OK && real == 42.0);
  CHECK(obdeect_stream_reader_meta_string(reader, "site", &text, NULL) == OBDEECT_STREAM_OK);
  CHECK(strcmp(text, "La Palma") == 0);
  CHECK(obdeect_stream_reader_meta_array(reader, "telescope_x_m", &array, &n) == OBDEECT_STREAM_OK);
  CHECK(n == 3 && array[1] == 12.5 && array[2] == -3.0);
  CHECK(obdeect_stream_reader_meta_int(reader, "site", &integer) == OBDEECT_STREAM_ERROR_TYPE);
  CHECK(obdeect_stream_reader_meta_int(reader, "missing", &integer) == OBDEECT_STREAM_ERROR_NOT_FOUND);
  CHECK(obdeect_stream_reader_photons(reader, NULL, NULL, NULL, NULL, NULL) == OBDEECT_STREAM_ERROR_STATE);
  CHECK(obdeect_stream_reader_column(reader, OBDEECT_STREAM_COLUMN_X) == NULL);

  expect_next(reader, OBDEECT_STREAM_RECORD_EVENT);
  CHECK(obdeect_stream_reader_n_keys(reader, &n) == OBDEECT_STREAM_OK && n == 0);
  expect_next(reader, OBDEECT_STREAM_RECORD_TELESCOPE);
  CHECK(obdeect_stream_reader_meta_int(reader, "telescope_id", &integer) == OBDEECT_STREAM_OK && integer == 3);

  expect_next(reader, OBDEECT_STREAM_RECORD_PHOTONS);
  CHECK(obdeect_stream_reader_photons(reader, &telescope, &frame, &first, &rows, &mask) == OBDEECT_STREAM_OK);
  CHECK(telescope == 3 && frame == OBDEECT_STREAM_FRAME_GROUND && first == 100 && rows == ROWS);
  for (i = 0; i < OBDEECT_STREAM_N_COLUMNS; ++i) {
    const float *column = obdeect_stream_reader_column(reader, (int)i);
    const int present = (mask >> i) & 1U;
    CHECK((column != NULL) == present);
    for (row = 0; present && row < ROWS; ++row) {
      CHECK(column[row] == value_of((int)i, row));
    }
  }
  CHECK(mask ==
        ((1U << OBDEECT_STREAM_COLUMN_X) | (1U << OBDEECT_STREAM_COLUMN_Y) | (1U << OBDEECT_STREAM_COLUMN_CX) |
         (1U << OBDEECT_STREAM_COLUMN_CY) | (1U << OBDEECT_STREAM_COLUMN_T) | (1U << OBDEECT_STREAM_COLUMN_W) |
         (1U << OBDEECT_STREAM_COLUMN_WAVELENGTH)));
  CHECK(obdeect_stream_reader_n_keys(reader, &n) == OBDEECT_STREAM_ERROR_STATE);

  expect_next(reader, OBDEECT_STREAM_RECORD_PHOTONS);
  CHECK(obdeect_stream_reader_photons(reader, &telescope, &frame, &first, &rows, &mask) == OBDEECT_STREAM_OK);
  CHECK(frame == OBDEECT_STREAM_FRAME_CAMERA && first == 105 && rows == 0 && mask == 0);
  CHECK(obdeect_stream_reader_column(reader, OBDEECT_STREAM_COLUMN_X) == NULL);

  expect_next(reader, OBDEECT_STREAM_RECORD_BLOB);
  CHECK(obdeect_stream_reader_blob(reader, &eventio_type, &eventio_id, &bytes, &n) == OBDEECT_STREAM_OK);
  CHECK(eventio_type == 1204 && eventio_id == 7 && n == 4 && bytes[3] == 250);

  expect_next(reader, OBDEECT_STREAM_RECORD_TELESCOPE_END);
  CHECK(obdeect_stream_reader_meta_float(reader, "photons", &real) == OBDEECT_STREAM_OK && real == 5.0);
  expect_next(reader, OBDEECT_STREAM_RECORD_RUN_END);
  CHECK(obdeect_stream_reader_next(reader, &type) == OBDEECT_STREAM_END && type == 0);
}

static void test_arguments(const char *directory) {
  obdeect_stream_writer *writer = NULL;
  obdeect_stream_reader *reader = NULL;
  char path[1024];
  float column[2] = {0.0F, 0.0F};
  const float *columns[OBDEECT_STREAM_N_COLUMNS] = {0};
  float *many = NULL;
  int type = 0;

  CHECK(obdeect_stream_abi_version() == OBDEECT_STREAM_C_ABI_VERSION);
  CHECK(obdeect_stream_format_major() == 1);
  CHECK(strlen(obdeect_stream_version()) > 0);
  CHECK(strcmp(obdeect_stream_record_type_name(OBDEECT_STREAM_RECORD_PHOTONS), "photons") == 0);
  CHECK(strcmp(obdeect_stream_column_name(OBDEECT_STREAM_COLUMN_WAVELENGTH), "wavelength") == 0);
  CHECK(obdeect_stream_column_name(11) == NULL && obdeect_stream_record_type_name(77) == NULL);

  CHECK(obdeect_stream_reader_open("/nonexistent/x.ps", &reader) == OBDEECT_STREAM_ERROR_IO);
  CHECK(reader == NULL && strlen(obdeect_stream_last_error()) > 0);
  CHECK(obdeect_stream_reader_open(NULL, &reader) == OBDEECT_STREAM_ERROR_ARGUMENT);
  CHECK(obdeect_stream_writer_open("/nonexistent/x.ps", "t", &writer) == OBDEECT_STREAM_ERROR_IO);
  obdeect_stream_reader_close(NULL);
  CHECK(obdeect_stream_writer_close(NULL) == OBDEECT_STREAM_OK);

  snprintf(path, sizeof path, "%.500s/arguments.ps", directory);
  CHECK(obdeect_stream_writer_open(path, "t", &writer) == OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_write_record(writer, OBDEECT_STREAM_RECORD_PHOTONS) == OBDEECT_STREAM_ERROR_ARGUMENT);
  CHECK(obdeect_stream_writer_write_record(writer, OBDEECT_STREAM_RECORD_STREAM) == OBDEECT_STREAM_ERROR_ARGUMENT);
  CHECK(obdeect_stream_writer_write_photons(writer, 1, 4, 0, 2, columns) == OBDEECT_STREAM_ERROR_ARGUMENT);
  /* Required columns missing. */
  columns[OBDEECT_STREAM_COLUMN_X] = column;
  CHECK(obdeect_stream_writer_write_photons(writer, 1, 1, 0, 2, columns) == OBDEECT_STREAM_ERROR_ARGUMENT);
  CHECK(obdeect_stream_writer_write_photons(writer, 1, 1, 0, 2, NULL) == OBDEECT_STREAM_ERROR_ARGUMENT);
  /* More than 2^16 rows. */
  many = calloc((size_t)OBDEECT_STREAM_MAX_CHUNK_ROWS + 1, sizeof(float));
  CHECK(many != NULL);
  if (many != NULL) {
    int c;
    for (c = 0; c < OBDEECT_STREAM_N_COLUMNS; ++c) {
      columns[c] = many;
    }
    CHECK(obdeect_stream_writer_write_photons(writer, 1, 1, 0, OBDEECT_STREAM_MAX_CHUNK_ROWS + 1, columns) ==
          OBDEECT_STREAM_ERROR_ARGUMENT);
    CHECK(obdeect_stream_writer_write_photons(writer, 1, 1, 0, OBDEECT_STREAM_MAX_CHUNK_ROWS, columns) ==
          OBDEECT_STREAM_OK);
    free(many);
  }
  /* Empty array and empty blob with NULL pointers. */
  CHECK(obdeect_stream_writer_set_array(writer, "empty", NULL, 0) == OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_write_record(writer, OBDEECT_STREAM_RECORD_RUN) == OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_write_blob(writer, 1, 2, NULL, 0) == OBDEECT_STREAM_OK);
  CHECK(obdeect_stream_writer_close(writer) == OBDEECT_STREAM_OK);

  CHECK(obdeect_stream_reader_open(path, &reader) == OBDEECT_STREAM_OK);
  expect_next(reader, OBDEECT_STREAM_RECORD_STREAM);
  expect_next(reader, OBDEECT_STREAM_RECORD_PHOTONS);
  expect_next(reader, OBDEECT_STREAM_RECORD_RUN);
  {
    const double *values = NULL;
    size_t n = 7;
    CHECK(obdeect_stream_reader_meta_array(reader, "empty", &values, &n) == OBDEECT_STREAM_OK && n == 0);
  }
  expect_next(reader, OBDEECT_STREAM_RECORD_BLOB);
  CHECK(obdeect_stream_reader_next(reader, &type) == OBDEECT_STREAM_END);
  obdeect_stream_reader_close(reader);
}

/* Every stream in the conformance set: the valid one ends cleanly after 11 records, the
 * invalid ones fail with the given code. */
static int count_records(const char *path, int *last) {
  obdeect_stream_reader *reader = NULL;
  int type = 0, code, count = 0;
  if (obdeect_stream_reader_open(path, &reader) != OBDEECT_STREAM_OK) {
    *last = OBDEECT_STREAM_ERROR_IO;
    return -1;
  }
  while ((code = obdeect_stream_reader_next(reader, &type)) == OBDEECT_STREAM_OK) {
    ++count;
  }
  *last = code;
  obdeect_stream_reader_close(reader);
  return count;
}

static void test_conformance(const char *directory) {
  const char *invalid[] = {"bad_magic.ps", "unknown_major.ps", "truncated_chunk.ps",
                           "too_many_rows.ps"};
  const int codes[] = {OBDEECT_STREAM_ERROR_FORMAT, OBDEECT_STREAM_ERROR_VERSION, OBDEECT_STREAM_ERROR_FORMAT,
                       OBDEECT_STREAM_ERROR_FORMAT};
  char path[1024];
  int last = 0, i;

  snprintf(path, sizeof path, "%.500s/valid.ps", directory);
  CHECK(count_records(path, &last) == 11 && last == OBDEECT_STREAM_END);
  for (i = 0; i < 4; ++i) {
    snprintf(path, sizeof path, "%.500s/%.64s", directory, invalid[i]);
    count_records(path, &last);
    CHECK(last == codes[i]);
  }
}

int main(int argc, char **argv) {
  obdeect_stream_writer *writer = NULL;
  obdeect_stream_reader *reader = NULL;
  char directory[512];
  const char *slash;

  if (argc != 3) {
    fprintf(stderr, "usage: test_capi <output.ps> <conformance dir>\n");
    return 2;
  }
  slash = strrchr(argv[1], '/');
  snprintf(directory, sizeof directory, "%.*s", slash ? (int)(slash - argv[1]) : 1,
           slash ? argv[1] : ".");

  CHECK(obdeect_stream_writer_open(argv[1], "test_capi", &writer) == OBDEECT_STREAM_OK);
  if (writer != NULL) {
    write_stream(writer);
    CHECK(obdeect_stream_writer_close(writer) == OBDEECT_STREAM_OK);
  }
  CHECK(obdeect_stream_reader_open(argv[1], &reader) == OBDEECT_STREAM_OK);
  if (reader != NULL) {
    check_stream(reader);
    obdeect_stream_reader_close(reader);
  }

  { /* The same through file descriptors. */
    char path[1100];
    int fd;
    snprintf(path, sizeof path, "%.1000s.fd", argv[1]);
    fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    CHECK(fd >= 0);
    CHECK(obdeect_stream_writer_open_fd(fd, "test_capi", &writer) == OBDEECT_STREAM_OK);
    write_stream(writer);
    CHECK(obdeect_stream_writer_close(writer) == OBDEECT_STREAM_OK);
    CHECK(close(fd) == 0);
    fd = open(path, O_RDONLY);
    CHECK(fd >= 0);
    CHECK(obdeect_stream_reader_open_fd(fd, &reader) == OBDEECT_STREAM_OK);
    check_stream(reader);
    obdeect_stream_reader_close(reader);
    close(fd);
    CHECK(obdeect_stream_reader_open_fd(-1, &reader) == OBDEECT_STREAM_ERROR_ARGUMENT);
  }

  test_arguments(directory);
  test_conformance(argv[2]);

  if (failures != 0) {
    fprintf(stderr, "%d check(s) failed\n", failures);
    return 1;
  }
  return 0;
}
