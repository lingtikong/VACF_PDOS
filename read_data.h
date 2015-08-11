#ifndef READ_DATA_H
#define READ_DATA_H

#include <stdio.h>
#include <stdlib.h>
#include "memory.h"

using namespace std;

class ReadVel {
public:
  ReadVel(char *, int);
  ~ReadVel();

  int natom, ntm, nmax, nstep;
  Memory *memory;
  double **velx, **vely, **velz;

  int *steps;
};

#endif
