/*
 * gstrickd C API: read and write photon streams (docs/format.md) from C and from any
 * language with a C FFI (Rust, Julia, Python ctypes, ...). Pure C, no exceptions cross
 * the boundary.
 *
 * Conventions
 * - Functions return GSTRICKD_OK (0), GSTRICKD_END (1, only gstrickd_reader_next) or a
 *   negative GSTRICKD_ERROR_* code. After a negative code gstrickd_last_error() returns a
 *   message for the calling thread.
 * - Handles are opaque. A handle may be used from one thread at a time; different handles
 *   may be used from different threads.
 * - Ownership: every pointer returned by the library (strings, arrays, photon columns,
 *   blob bytes, keys) belongs to the library. It stays valid until the next call of
 *   gstrickd_reader_next() or gstrickd_reader_close() on the same reader, whichever comes
 *   first. Copy what you need to keep. Pointers you pass in are only read during the call.
 * - Photon columns are float32, little-endian host values, one array of `rows` values per
 *   column; no copy is made. An absent column is NULL.
 * - All integers and floats of the stream are converted to host byte order by the library.
 */
#ifndef GSTRICKD_GSTRICKD_H
#define GSTRICKD_GSTRICKD_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32)
#if defined(GSTRICKD_C_BUILDING)
#define GSTRICKD_C_API __declspec(dllexport)
#elif defined(GSTRICKD_C_STATIC)
#define GSTRICKD_C_API
#else
#define GSTRICKD_C_API __declspec(dllimport)
#endif
#else
#define GSTRICKD_C_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Version of this C ABI. Incremented when an existing function changes; new functions
 * do not change it. Compare gstrickd_abi_version() with it at start-up. */
#define GSTRICKD_C_ABI_VERSION 1

/* Maximum number of photons in one chunk. */
#define GSTRICKD_MAX_CHUNK_ROWS 65536

/* Return codes. */
enum {
  GSTRICKD_OK = 0,
  GSTRICKD_END = 1,              /* gstrickd_reader_next: stream_end reached */
  GSTRICKD_ERROR_ARGUMENT = -1,  /* NULL pointer, bad enum value, invalid chunk or metadata */
  GSTRICKD_ERROR_IO = -2,        /* cannot open, read or write */
  GSTRICKD_ERROR_FORMAT = -3,    /* malformed or truncated stream */
  GSTRICKD_ERROR_VERSION = -4,   /* stream has an unsupported format major version */
  GSTRICKD_ERROR_NOT_FOUND = -5, /* metadata key does not exist */
  GSTRICKD_ERROR_TYPE = -6,      /* metadata value has another type */
  GSTRICKD_ERROR_STATE = -7,     /* no current record of the required kind */
  GSTRICKD_ERROR_MEMORY = -8,    /* out of memory */
  GSTRICKD_ERROR_INTERNAL = -9   /* unexpected failure */
};

/* Record types (values as in docs/format.md). */
enum {
  GSTRICKD_RECORD_STREAM = 1,
  GSTRICKD_RECORD_RUN = 2,
  GSTRICKD_RECORD_EVENT = 3,
  GSTRICKD_RECORD_TELESCOPE = 4,
  GSTRICKD_RECORD_PHOTONS = 5,
  GSTRICKD_RECORD_TELESCOPE_END = 6,
  GSTRICKD_RECORD_EVENT_END = 7,
  GSTRICKD_RECORD_RUN_END = 8,
  GSTRICKD_RECORD_BLOB = 9,
  GSTRICKD_RECORD_LEDGER = 10
};

/* Photon frames. */
enum {
  GSTRICKD_FRAME_SKY = 0,
  GSTRICKD_FRAME_GROUND = 1,
  GSTRICKD_FRAME_TELESCOPE = 2,
  GSTRICKD_FRAME_CAMERA = 3
};

/* Photon columns: index into the column arrays, bit (1 << column) in the column mask.
 * x y cx cy t w are required in every chunk with rows. */
enum {
  GSTRICKD_COLUMN_X = 0,
  GSTRICKD_COLUMN_Y = 1,
  GSTRICKD_COLUMN_Z = 2,
  GSTRICKD_COLUMN_CX = 3,
  GSTRICKD_COLUMN_CY = 4,
  GSTRICKD_COLUMN_CZ = 5,
  GSTRICKD_COLUMN_T = 6,
  GSTRICKD_COLUMN_ZEM = 7,
  GSTRICKD_COLUMN_W = 8,
  GSTRICKD_COLUMN_WAVELENGTH = 9,
  GSTRICKD_COLUMN_P = 10,
  GSTRICKD_N_COLUMNS = 11
};

/* Metadata value types. */
enum {
  GSTRICKD_VALUE_INT = 1,
  GSTRICKD_VALUE_FLOAT = 2,
  GSTRICKD_VALUE_STRING = 3,
  GSTRICKD_VALUE_ARRAY = 4 /* array of float64 */
};

typedef struct gstrickd_reader gstrickd_reader;
typedef struct gstrickd_writer gstrickd_writer;

/* ---- Versions and errors ---- */

/* Package version, e.g. "0.1.0". Static string. */
GSTRICKD_C_API const char *gstrickd_version(void);
/* Value of GSTRICKD_C_ABI_VERSION the library was built with. */
GSTRICKD_C_API int gstrickd_abi_version(void);
/* Photon-stream format version the library reads and writes (readers accept any minor). */
GSTRICKD_C_API int gstrickd_format_major(void);
GSTRICKD_C_API int gstrickd_format_minor(void);
/* Message of the last failed call in the calling thread; never NULL. Valid until the next
 * failing call in this thread. Not cleared by successful calls. */
GSTRICKD_C_API const char *gstrickd_last_error(void);
/* Names ("run", "photons", ...; "x", "cx", ...), or NULL for an unknown value. Static. */
GSTRICKD_C_API const char *gstrickd_record_type_name(int record_type);
GSTRICKD_C_API const char *gstrickd_column_name(int column);

/* ---- Reader ---- */

/* Opens a stream file. On success *out is a new reader that the caller closes with
 * gstrickd_reader_close(). */
GSTRICKD_C_API int gstrickd_reader_open(const char *path, gstrickd_reader **out);
/* Reads from an open file descriptor (e.g. 0 for stdin; a pipe works). The reader never
 * closes the descriptor; the caller does, after gstrickd_reader_close(). Reading is
 * buffered, so do not read from fd yourself while the reader is in use. */
GSTRICKD_C_API int gstrickd_reader_open_fd(int fd, gstrickd_reader **out);
/* Frees the reader and everything it returned. NULL is allowed. */
GSTRICKD_C_API void gstrickd_reader_close(gstrickd_reader *reader);

/* Advances to the next record and stores its type in *record_type. The first record is
 * GSTRICKD_RECORD_STREAM. Returns GSTRICKD_END (and *record_type = 0) after the
 * stream_end record; a stream that ends without it, or has data after it, is
 * GSTRICKD_ERROR_FORMAT. Unknown record types are skipped. A failed reader stays failed.
 * Invalidates all pointers returned for the previous record. */
GSTRICKD_C_API int gstrickd_reader_next(gstrickd_reader *reader, int *record_type);

/* Metadata of the current record (stream, run, event, telescope, telescope_end,
 * event_end, run_end, ledger); GSTRICKD_ERROR_STATE for other records. Keys are sorted
 * and unique. */
GSTRICKD_C_API int gstrickd_reader_n_keys(const gstrickd_reader *reader, size_t *n);
/* Key i (0 <= i < n) as a NUL-terminated UTF-8 string, or NULL if out of range. */
GSTRICKD_C_API const char *gstrickd_reader_key(const gstrickd_reader *reader, size_t i);
/* Type of the value of `key`: GSTRICKD_VALUE_* in *value_type. */
GSTRICKD_C_API int gstrickd_reader_meta_type(const gstrickd_reader *reader, const char *key,
                                             int *value_type);
GSTRICKD_C_API int gstrickd_reader_meta_int(const gstrickd_reader *reader, const char *key,
                                            int64_t *value);
/* Accepts float values and, converted, int values. */
GSTRICKD_C_API int gstrickd_reader_meta_float(const gstrickd_reader *reader, const char *key,
                                              double *value);
/* *value is NUL-terminated; *length (optional, may be NULL) is the length without the NUL. */
GSTRICKD_C_API int gstrickd_reader_meta_string(const gstrickd_reader *reader, const char *key,
                                               const char **value, size_t *length);
GSTRICKD_C_API int gstrickd_reader_meta_array(const gstrickd_reader *reader, const char *key,
                                              const double **values, size_t *n);

/* Current photons record (GSTRICKD_ERROR_STATE otherwise). Every output pointer may be
 * NULL. column_mask has bit (1 << GSTRICKD_COLUMN_*) set for each present column. */
GSTRICKD_C_API int gstrickd_reader_photons(const gstrickd_reader *reader, uint32_t *telescope_id,
                                           uint32_t *frame, uint64_t *first_index, uint32_t *rows,
                                           uint32_t *column_mask);
/* Pointer to `rows` float32 values of a column of the current photons record, or NULL if
 * the column is absent, the chunk is empty or the current record is not photons. Valid
 * until the next call of gstrickd_reader_next() or gstrickd_reader_close(). */
GSTRICKD_C_API const float *gstrickd_reader_column(const gstrickd_reader *reader, int column);

/* Current blob record: raw EventIO block (header and body). */
GSTRICKD_C_API int gstrickd_reader_blob(const gstrickd_reader *reader, uint32_t *eventio_type,
                                        uint32_t *eventio_id, const uint8_t **data, size_t *size);

/* ---- Writer ---- */

/* Creates a stream file and writes the stream record (producer: non-NULL name of the
 * writing program). The caller finishes with gstrickd_writer_close(). */
GSTRICKD_C_API int gstrickd_writer_open(const char *path, const char *producer,
                                        gstrickd_writer **out);
/* Writes to an open file descriptor; the writer never closes it. */
GSTRICKD_C_API int gstrickd_writer_open_fd(int fd, const char *producer, gstrickd_writer **out);

/* Pending metadata of the next gstrickd_writer_write_record(). Setting a key again
 * replaces its value. Strings are copied (NUL-terminated); arrays are copied. */
GSTRICKD_C_API int gstrickd_writer_set_int(gstrickd_writer *writer, const char *key, int64_t value);
GSTRICKD_C_API int gstrickd_writer_set_float(gstrickd_writer *writer, const char *key,
                                             double value);
GSTRICKD_C_API int gstrickd_writer_set_string(gstrickd_writer *writer, const char *key,
                                              const char *value);
GSTRICKD_C_API int gstrickd_writer_set_array(gstrickd_writer *writer, const char *key,
                                             const double *values, size_t n);

/* Writes a metadata record (RUN, EVENT, TELESCOPE, TELESCOPE_END, EVENT_END, RUN_END or
 * LEDGER) with the pending metadata and clears it, also when the call fails. */
GSTRICKD_C_API int gstrickd_writer_write_record(gstrickd_writer *writer, int record_type);

/* Writes one photon chunk. columns points to GSTRICKD_N_COLUMNS pointers, each NULL
 * (absent) or `rows` float32 values; x y cx cy t w are required if rows > 0. For rows == 0
 * an empty chunk is written and columns may be NULL. At most GSTRICKD_MAX_CHUNK_ROWS rows.
 * The values are copied; the arrays may be reused after the call. */
GSTRICKD_C_API int gstrickd_writer_write_photons(gstrickd_writer *writer, uint32_t telescope_id,
                                                 uint32_t frame, uint64_t first_index,
                                                 uint32_t rows, const float *const *columns);

/* Writes a blob record: a raw EventIO block (header and body). data may be NULL if size is 0. */
GSTRICKD_C_API int gstrickd_writer_write_blob(gstrickd_writer *writer, uint32_t eventio_type,
                                              uint32_t eventio_id, const uint8_t *data,
                                              size_t size);

/* Writes stream_end, flushes, closes a file opened by path and frees the writer, also
 * when it fails; the handle is invalid afterwards. NULL is allowed. Returns an error if any
 * earlier or this write failed (e.g. disk full). */
GSTRICKD_C_API int gstrickd_writer_close(gstrickd_writer *writer);

#ifdef __cplusplus
}
#endif

#endif /* GSTRICKD_GSTRICKD_H */
