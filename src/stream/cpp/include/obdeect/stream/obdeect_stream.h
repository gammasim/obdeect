/*
 * gstrickd C API: read and write photon streams (docs/format.md) from C and from any
 * language with a C FFI (Rust, Julia, Python ctypes, ...). Pure C, no exceptions cross
 * the boundary.
 *
 * Conventions
 * - Functions return OBDEECT_STREAM_OK (0), OBDEECT_STREAM_END (1, only obdeect_stream_reader_next) or a
 *   negative OBDEECT_STREAM_ERROR_* code. After a negative code obdeect_stream_last_error() returns a
 *   message for the calling thread.
 * - Handles are opaque. A handle may be used from one thread at a time; different handles
 *   may be used from different threads.
 * - Ownership: every pointer returned by the library (strings, arrays, photon columns,
 *   blob bytes, keys) belongs to the library. It stays valid until the next call of
 *   obdeect_stream_reader_next() or obdeect_stream_reader_close() on the same reader, whichever comes
 *   first. Copy what you need to keep. Pointers you pass in are only read during the call.
 * - Photon columns are float32, little-endian host values, one array of `rows` values per
 *   column; no copy is made. An absent column is NULL.
 * - All integers and floats of the stream are converted to host byte order by the library.
 */
#ifndef OBDEECT_STREAM_GSTRICKD_H
#define OBDEECT_STREAM_GSTRICKD_H

#include <stddef.h>
#include <stdint.h>

#define OBDEECT_STREAM_C_API __attribute__((visibility("default")))

#ifdef __cplusplus
extern "C" {
#endif

/* Version of this C ABI. Incremented when an existing function changes; new functions
 * do not change it. Compare obdeect_stream_abi_version() with it at start-up. */
#define OBDEECT_STREAM_C_ABI_VERSION 1

/* Maximum number of photons in one chunk. */
#define OBDEECT_STREAM_MAX_CHUNK_ROWS 65536

/* Return codes. */
enum {
  OBDEECT_STREAM_OK = 0,
  OBDEECT_STREAM_END = 1,              /* obdeect_stream_reader_next: stream_end reached */
  OBDEECT_STREAM_ERROR_ARGUMENT = -1,  /* NULL pointer, bad enum value, invalid chunk or metadata */
  OBDEECT_STREAM_ERROR_IO = -2,        /* cannot open, read or write */
  OBDEECT_STREAM_ERROR_FORMAT = -3,    /* malformed or truncated stream */
  OBDEECT_STREAM_ERROR_VERSION = -4,   /* stream has an unsupported format major version */
  OBDEECT_STREAM_ERROR_NOT_FOUND = -5, /* metadata key does not exist */
  OBDEECT_STREAM_ERROR_TYPE = -6,      /* metadata value has another type */
  OBDEECT_STREAM_ERROR_STATE = -7,     /* no current record of the required kind */
  OBDEECT_STREAM_ERROR_MEMORY = -8,    /* out of memory */
  OBDEECT_STREAM_ERROR_INTERNAL = -9   /* unexpected failure */
};

/* Record types (values as in docs/format.md). */
enum {
  OBDEECT_STREAM_RECORD_STREAM = 1,
  OBDEECT_STREAM_RECORD_RUN = 2,
  OBDEECT_STREAM_RECORD_EVENT = 3,
  OBDEECT_STREAM_RECORD_TELESCOPE = 4,
  OBDEECT_STREAM_RECORD_PHOTONS = 5,
  OBDEECT_STREAM_RECORD_TELESCOPE_END = 6,
  OBDEECT_STREAM_RECORD_EVENT_END = 7,
  OBDEECT_STREAM_RECORD_RUN_END = 8,
  OBDEECT_STREAM_RECORD_BLOB = 9,
  OBDEECT_STREAM_RECORD_EFFICIENCY_STATUS = 10
};

/* Photon frames. */
enum {
  OBDEECT_STREAM_FRAME_SKY = 0,
  OBDEECT_STREAM_FRAME_GROUND = 1,
  OBDEECT_STREAM_FRAME_TELESCOPE = 2,
  OBDEECT_STREAM_FRAME_CAMERA = 3
};

/* Photon columns: index into the column arrays, bit (1 << column) in the column mask.
 * x y cx cy t w are required in every chunk with rows. */
enum {
  OBDEECT_STREAM_COLUMN_X = 0,
  OBDEECT_STREAM_COLUMN_Y = 1,
  OBDEECT_STREAM_COLUMN_Z = 2,
  OBDEECT_STREAM_COLUMN_CX = 3,
  OBDEECT_STREAM_COLUMN_CY = 4,
  OBDEECT_STREAM_COLUMN_CZ = 5,
  OBDEECT_STREAM_COLUMN_T = 6,
  OBDEECT_STREAM_COLUMN_ZEM = 7,
  OBDEECT_STREAM_COLUMN_W = 8,
  OBDEECT_STREAM_COLUMN_WAVELENGTH = 9,
  OBDEECT_STREAM_COLUMN_P = 10,
  OBDEECT_STREAM_N_COLUMNS = 11
};

/* Metadata value types. */
enum {
  OBDEECT_STREAM_VALUE_INT = 1,
  OBDEECT_STREAM_VALUE_FLOAT = 2,
  OBDEECT_STREAM_VALUE_STRING = 3,
  OBDEECT_STREAM_VALUE_ARRAY = 4 /* array of float64 */
};

typedef struct obdeect_stream_reader obdeect_stream_reader;
typedef struct obdeect_stream_writer obdeect_stream_writer;

/* ---- Versions and errors ---- */

/* Package version, e.g. "0.1.0". Static string. */
OBDEECT_STREAM_C_API const char *obdeect_stream_version(void);
/* Value of OBDEECT_STREAM_C_ABI_VERSION the library was built with. */
OBDEECT_STREAM_C_API int obdeect_stream_abi_version(void);
/* Photon-stream format version the library reads and writes (readers accept any minor). */
OBDEECT_STREAM_C_API int obdeect_stream_format_major(void);
OBDEECT_STREAM_C_API int obdeect_stream_format_minor(void);
/* Message of the last failed call in the calling thread; never NULL. Valid until the next
 * failing call in this thread. Not cleared by successful calls. */
OBDEECT_STREAM_C_API const char *obdeect_stream_last_error(void);
/* Names ("run", "photons", ...; "x", "cx", ...), or NULL for an unknown value. Static. */
OBDEECT_STREAM_C_API const char *obdeect_stream_record_type_name(int record_type);
OBDEECT_STREAM_C_API const char *obdeect_stream_column_name(int column);

/* ---- Reader ---- */

/* Opens a stream file. On success *out is a new reader that the caller closes with
 * obdeect_stream_reader_close(). */
OBDEECT_STREAM_C_API int obdeect_stream_reader_open(const char *path, obdeect_stream_reader **out);
/* Reads from an open file descriptor (e.g. 0 for stdin; a pipe works). The reader never
 * closes the descriptor; the caller does, after obdeect_stream_reader_close(). Reading is
 * buffered, so do not read from fd yourself while the reader is in use. */
OBDEECT_STREAM_C_API int obdeect_stream_reader_open_fd(int fd, obdeect_stream_reader **out);
/* Frees the reader and everything it returned. NULL is allowed. */
OBDEECT_STREAM_C_API void obdeect_stream_reader_close(obdeect_stream_reader *reader);

/* Advances to the next record and stores its type in *record_type. The first record is
 * OBDEECT_STREAM_RECORD_STREAM. Returns OBDEECT_STREAM_END (and *record_type = 0) after the
 * stream_end record; a stream that ends without it, or has data after it, is
 * OBDEECT_STREAM_ERROR_FORMAT. Unknown record types are skipped. A failed reader stays failed.
 * Invalidates all pointers returned for the previous record. */
OBDEECT_STREAM_C_API int obdeect_stream_reader_next(obdeect_stream_reader *reader, int *record_type);

/* Metadata of the current record (stream, run, event, telescope, telescope_end,
 * event_end, run_end, efficiency_status); OBDEECT_STREAM_ERROR_STATE for other records. Keys are sorted
 * and unique. */
OBDEECT_STREAM_C_API int obdeect_stream_reader_n_keys(const obdeect_stream_reader *reader, size_t *n);
/* Key i (0 <= i < n) as a NUL-terminated UTF-8 string, or NULL if out of range. */
OBDEECT_STREAM_C_API const char *obdeect_stream_reader_key(const obdeect_stream_reader *reader, size_t i);
/* Type of the value of `key`: OBDEECT_STREAM_VALUE_* in *value_type. */
OBDEECT_STREAM_C_API int obdeect_stream_reader_meta_type(const obdeect_stream_reader *reader, const char *key,
                                             int *value_type);
OBDEECT_STREAM_C_API int obdeect_stream_reader_meta_int(const obdeect_stream_reader *reader, const char *key,
                                            int64_t *value);
/* Accepts float values and, converted, int values. */
OBDEECT_STREAM_C_API int obdeect_stream_reader_meta_float(const obdeect_stream_reader *reader, const char *key,
                                              double *value);
/* *value is NUL-terminated; *length (optional, may be NULL) is the length without the NUL. */
OBDEECT_STREAM_C_API int obdeect_stream_reader_meta_string(const obdeect_stream_reader *reader, const char *key,
                                               const char **value, size_t *length);
OBDEECT_STREAM_C_API int obdeect_stream_reader_meta_array(const obdeect_stream_reader *reader, const char *key,
                                              const double **values, size_t *n);

/* Current photons record (OBDEECT_STREAM_ERROR_STATE otherwise). Every output pointer may be
 * NULL. column_mask has bit (1 << OBDEECT_STREAM_COLUMN_*) set for each present column. */
OBDEECT_STREAM_C_API int obdeect_stream_reader_photons(const obdeect_stream_reader *reader, uint32_t *telescope_id,
                                           uint32_t *frame, uint64_t *first_index, uint32_t *rows,
                                           uint32_t *column_mask);
/* Pointer to `rows` float32 values of a column of the current photons record, or NULL if
 * the column is absent, the chunk is empty or the current record is not photons. Valid
 * until the next call of obdeect_stream_reader_next() or obdeect_stream_reader_close(). */
OBDEECT_STREAM_C_API const float *obdeect_stream_reader_column(const obdeect_stream_reader *reader, int column);

/* Current blob record: raw EventIO block (header and body). */
OBDEECT_STREAM_C_API int obdeect_stream_reader_blob(const obdeect_stream_reader *reader, uint32_t *eventio_type,
                                        uint32_t *eventio_id, const uint8_t **data, size_t *size);

/* ---- Writer ---- */

/* Creates a stream file and writes the stream record (producer: non-NULL name of the
 * writing program). The caller finishes with obdeect_stream_writer_close(). */
OBDEECT_STREAM_C_API int obdeect_stream_writer_open(const char *path, const char *producer,
                                        obdeect_stream_writer **out);
/* Writes to an open file descriptor; the writer never closes it. */
OBDEECT_STREAM_C_API int obdeect_stream_writer_open_fd(int fd, const char *producer, obdeect_stream_writer **out);

/* Pending metadata of the next obdeect_stream_writer_write_record(). Setting a key again
 * replaces its value. Strings are copied (NUL-terminated); arrays are copied. */
OBDEECT_STREAM_C_API int obdeect_stream_writer_set_int(obdeect_stream_writer *writer, const char *key, int64_t value);
OBDEECT_STREAM_C_API int obdeect_stream_writer_set_float(obdeect_stream_writer *writer, const char *key,
                                             double value);
OBDEECT_STREAM_C_API int obdeect_stream_writer_set_string(obdeect_stream_writer *writer, const char *key,
                                              const char *value);
OBDEECT_STREAM_C_API int obdeect_stream_writer_set_array(obdeect_stream_writer *writer, const char *key,
                                             const double *values, size_t n);

/* Writes a metadata record (RUN, EVENT, TELESCOPE, TELESCOPE_END, EVENT_END, RUN_END or
 * EFFICIENCY_STATUS) with the pending metadata and clears it, also when the call fails. */
OBDEECT_STREAM_C_API int obdeect_stream_writer_write_record(obdeect_stream_writer *writer, int record_type);

/* Writes one photon chunk. columns points to OBDEECT_STREAM_N_COLUMNS pointers, each NULL
 * (absent) or `rows` float32 values; x y cx cy t w are required if rows > 0. For rows == 0
 * an empty chunk is written and columns may be NULL. At most OBDEECT_STREAM_MAX_CHUNK_ROWS rows.
 * The values are copied; the arrays may be reused after the call. */
OBDEECT_STREAM_C_API int obdeect_stream_writer_write_photons(obdeect_stream_writer *writer, uint32_t telescope_id,
                                                 uint32_t frame, uint64_t first_index,
                                                 uint32_t rows, const float *const *columns);

/* Writes a blob record: a raw EventIO block (header and body). data may be NULL if size is 0. */
OBDEECT_STREAM_C_API int obdeect_stream_writer_write_blob(obdeect_stream_writer *writer, uint32_t eventio_type,
                                              uint32_t eventio_id, const uint8_t *data,
                                              size_t size);

/* Writes stream_end, flushes, closes a file opened by path and frees the writer, also
 * when it fails; the handle is invalid afterwards. NULL is allowed. Returns an error if any
 * earlier or this write failed (e.g. disk full). */
OBDEECT_STREAM_C_API int obdeect_stream_writer_close(obdeect_stream_writer *writer);

#ifdef __cplusplus
}
#endif

#endif /* OBDEECT_STREAM_GSTRICKD_H */
