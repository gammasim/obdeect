/*
 * Writes photons3d.iact: a small CORSIKA IACT EventIO file with 3D photon bunches (1217).
 *
 * Test fixture for obdeect stream from-corsika (cpp/tests/test_corsika.cpp, python tests).
 * NOT part of the obdeect build. It is compiled against the CORSIKA IACT/ATMO library
 * (Konrad Bernloehr, IACT/ATMO package 1.69, corsika7/bernlohr), so the 1217 layout is the one
 * that write_tel_photons3d() of io_simtel.c writes. See photons3d_provenance.txt.
 *
 * Build and run from the obdeect repository root:
 *
 *   cc -std=c99 -O1 -w -I../corsika7/bernlohr -o /tmp/make_photons3d_fixture \
 *      src/stream/cpp/tests/data/make_photons3d_fixture.c \
 *      ../corsika7/bernlohr/eventio.c ../corsika7/bernlohr/io_simtel.c \
 *      ../corsika7/bernlohr/fileopen.c ../corsika7/bernlohr/warning.c \
 *      ../corsika7/bernlohr/straux.c ../corsika7/bernlohr/mc_atmprof.c -lm
 *   /tmp/make_photons3d_fixture src/stream/cpp/tests/data/photons3d.iact
 * (atmo.c is not needed: it depends on CORSIKA Fortran routines.)
 *
 * All bunch values are CORSIKA units (cm, ns, nm). The expected stream values of every bunch
 * are derived by hand in test_corsika.cpp (photons3d_fixture_rows); keep them in sync.
 *
 * Telescope 1 (array 0, tel 0): position (0, 0, 1800) cm, sphere radius 1200 cm.
 * Telescope 2 (array 0, tel 1): position (1000, -500, 2000) cm, sphere radius 800 cm.
 * Observation level 220000 cm, standard (flat) geometry, primary gamma, 1 TeV, zenith 0.3 rad.
 */

#include <stdio.h>
#include <stdlib.h>

#include "mc_tel.h"

static void set_bunch(struct bunch3d *b, float photons, float x, float y, float z, float cx,
                      float cy, float cz, float ctime, float dist, float lambda) {
  b->photons = photons;
  b->x = x;
  b->y = y;
  b->z = z;
  b->cx = cx;
  b->cy = cy;
  b->cz = cz;
  b->ctime = ctime;
  b->dist = dist;
  b->lambda = lambda;
}

int main(int argc, char **argv) {
  FILE *out;
  IO_BUFFER *iobuf;
  IO_ITEM_HEADER array_header;
  float runh[273] = {0.0F};
  float evth[273] = {0.0F};
  float evte[273] = {0.0F};
  float rune[3] = {0.0F};
  double tel_x[2] = {0.0, 1000.0};
  double tel_y[2] = {0.0, -500.0};
  double tel_z[2] = {1800.0, 2000.0};
  double tel_r[2] = {1200.0, 800.0};
  double toff = 0.0;
  double xoff[1] = {0.0};
  double yoff[1] = {0.0};
  struct bunch3d tel0[8];
  struct bunch3d tel1[3];
  double sum0 = 0.0;
  double sum1 = 0.0;
  int i;

  if (argc != 2) {
    fprintf(stderr, "usage: %s OUTPUT.iact\n", argv[0]);
    return 2;
  }
  out = fopen(argv[1], "wb");
  if (out == NULL) {
    perror(argv[1]);
    return 1;
  }
  iobuf = allocate_io_buffer(100000);
  if (iobuf == NULL) {
    return 1;
  }
  iobuf->output_file = out;

  runh[1] = 1.0F;        /* run number */
  runh[3] = 7.41F;       /* version */
  runh[4] = 1.0F;        /* number of observation levels */
  runh[5] = 220000.0F;   /* observation level [cm] */
  runh[15] = -2.0F;      /* spectral slope */
  runh[16] = 10.0F;      /* energy min [GeV] */
  runh[17] = 100000.0F;  /* energy max [GeV] */
  evth[1] = 1.0F;        /* event number */
  evth[2] = 14.0F;       /* primary (gamma) */
  evth[3] = 1000.0F;     /* energy [GeV] */
  evth[6] = 2500000.0F;  /* first interaction [cm] */
  evth[10] = 0.3F;       /* zenith [rad] */
  evth[11] = 1.2F;       /* azimuth [rad] */
  evth[76] = 5136.0F;    /* atmosphere 5, refraction */
  evth[78] = 0.0F;       /* CURVED flag: standard */
  evth[95] = 250.0F;     /* wavelength range [nm] */
  evth[96] = 800.0F;
  evte[1] = 1.0F;
  rune[1] = 1.0F;
  rune[2] = 2.0F;

  /* Telescope 1, sphere radius 1200 cm. Columns: photons, x y z, cx cy cz, ctime, dist, lambda. */
  /* 0: downward, arrival on the sphere top; emitted 42 m above the centre. */
  set_bunch(&tel0[0], 2.0F, 0.0F, 0.0F, 1200.0F, 0.0F, 0.0F, -1.0F, 100.0F, 3000.0F, 0.0F);
  /* 1: downward, emitted inside the sphere (300, 0, 400) cm; arrival inside. */
  set_bunch(&tel0[1], 5.0F, 300.0F, 0.0F, -400.0F, 0.0F, 0.0F, -1.0F, 50.0F, 800.0F, 350.0F);
  /* 2: upward, emitted 30 m below the centre; the line enters the sphere from below. */
  set_bunch(&tel0[2], 3.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 1.0F, 200.0F, 3000.0F, 0.0F);
  /* 3: horizontal, emitted 20 m from the centre along the x axis. */
  set_bunch(&tel0[3], 4.0F, 0.0F, 0.0F, 0.0F, 1.0F, 0.0F, 0.0F, 300.0F, 2000.0F, 0.0F);
  /* 4: arrival (900, 0, -1200) cm outside the sphere; the line from (-900, 0, 1200) enters first. */
  set_bunch(&tel0[4], 6.0F, 900.0F, 0.0F, -1200.0F, 0.6F, 0.0F, -0.8F, 400.0F, 3000.0F, 0.0F);
  /* 5: misses the sphere (passes at y = 1500 cm). */
  set_bunch(&tel0[5], 9.0F, 0.0F, 1500.0F, 0.0F, 1.0F, 0.0F, 0.0F, 500.0F, 3000.0F, 0.0F);
  /* 6: misses the sphere (points away from it). */
  set_bunch(&tel0[6], 7.0F, 0.0F, 0.0F, 3000.0F, 0.0F, 0.0F, 1.0F, 600.0F, 1000.0F, 0.0F);
  /* 7: downward, arrival on the sphere; emitted 30 m above the centre at (0, 0, 3000) cm. */
  set_bunch(&tel0[7], 1.5F, 0.0F, 0.0F, 1200.0F, 0.0F, 0.0F, -1.0F, 700.0F, 1800.0F, 532.0F);

  /* Telescope 2, sphere radius 800 cm. */
  /* 0: downward, arrival on the sphere top, emission 28 m above the telescope centre. */
  set_bunch(&tel1[0], 1.0F, 0.0F, 0.0F, 800.0F, 0.0F, 0.0F, -1.0F, 120.0F, 2000.0F, 0.0F);
  /* 1: upward (cz = 0.8), emitted at (-940, 0, -880) cm, outside; enters at s = 500 cm. */
  set_bunch(&tel1[1], 4.0F, -340.0F, 0.0F, -80.0F, 0.6F, 0.0F, 0.8F, 150.0F, 1000.0F, 532.0F);
  /* 2: downward, emitted inside the sphere at (0, 0, 300) cm. */
  set_bunch(&tel1[2], 2.0F, 0.0F, 0.0F, -200.0F, 0.0F, 0.0F, -1.0F, 60.0F, 500.0F, 0.0F);

  for (i = 0; i < 8; i++) {
    sum0 += tel0[i].photons;
  }
  for (i = 0; i < 3; i++) {
    sum1 += tel1[i].photons;
  }

  if (write_tel_block(iobuf, IO_TYPE_MC_RUNH, 1, runh, 273) != 0 ||
      write_tel_pos(iobuf, 2, tel_x, tel_y, tel_z, tel_r) != 0 ||
      write_tel_block(iobuf, IO_TYPE_MC_EVTH, 1, evth, 273) != 0 ||
      write_tel_offset(iobuf, 1, toff, xoff, yoff) != 0) {
    fprintf(stderr, "write error (header)\n");
    return 1;
  }
  if (begin_write_tel_array(iobuf, &array_header, 0) != 0 ||
      write_tel_photons3d(iobuf, 0, 0, sum0, tel0, 8, 0, NULL) != 0 ||
      write_tel_photons3d(iobuf, 0, 1, sum1, tel1, 3, 0, NULL) != 0 ||
      end_write_tel_array(iobuf, &array_header) != 0) {
    fprintf(stderr, "write error (telescope array)\n");
    return 1;
  }
  if (write_tel_block(iobuf, IO_TYPE_MC_EVTE, 1, evte, 273) != 0 ||
      write_tel_block(iobuf, IO_TYPE_MC_RUNE, 1, rune, 3) != 0) {
    fprintf(stderr, "write error (trailer)\n");
    return 1;
  }
  fclose(out);
  return 0;
}
