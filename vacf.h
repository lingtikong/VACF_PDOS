#ifndef VACF_H
#define VACF_H

#include "complex.h"
#include "fftw3.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "read_data.h"
#include "memory.h"

class VACF {
public:
  VACF(int, char **);
  ~VACF();

  ReadVel *vels;
  Memory  *memory;

private:
  int ncount, nlag, sysdim;
  int dstep; // (steps[nalg] - steps[0])/nlag
  int dlag;
  double **sum, *vv0;
  double **pdos;
  double dt; // MD time step size in unit of ps
  int ismear; // smearing method
  double tau0, dtau; // smearing parameters, in unit of dt

  char  *outacf, *outdos;

  void read_acf(char *);
  void compute_acf(char *);
  double smearing(int);
  void compute_dos();

  int ndos;
  double vmin, vmax, dv, vstep;

  void write_acf();
  void write_dos();

  void help();
  int  count_words(const char *);
  void normal_dos();
  void normal_acf();
};

#endif
