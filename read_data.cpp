#include "read_data.h"
#include "string.h"

#define MAXLINE 256

/* -----------------------------------------------------------------------------
 * constructor, to read lammps dump file, only vx, vy, vz should be dumped
 * in order
 * -------------------------------------------------------------------------- */
ReadVel::ReadVel(char *infile, int sdim)
{
  velx = vely = velz = NULL;
  steps = NULL;
  memory = NULL;

  FILE *fp = fopen(infile,"r");
  if (fp == NULL) {
    printf("\nError while opening file %s to read!\n", infile);
    exit(1);
  }
  printf("\nNow to read file %s... ", infile); fflush(stdout);

  nstep = 0;
  memory = new Memory();
  // read the first image
  char str[MAXLINE], *ptr;
  fgets(str, MAXLINE, fp); fgets(str, MAXLINE, fp);
  int istep = atoi(strtok(str," \t\n\r\f"));
  fgets(str, MAXLINE, fp); fgets(str, MAXLINE, fp);
  ntm = atoi(strtok(str," \t\n\r\f"));

  natom = ntm; nmax = 200001;
  velx = memory->create(velx, nmax, natom, "ReadVel_velx");
  if (sdim > 1) vely = memory->create(vely, nmax, natom, "ReadVel_vely");
  if (sdim > 2) velz = memory->create(velz, nmax, natom, "ReadVel_velz");
  steps = memory->create(steps, nmax, "ReadVel_steps");

  double vx, vy, vz;
  double *xsum, *ysum, *zsum;
  xsum = memory->create(xsum, natom, "xsum");
  ysum = memory->create(ysum, natom, "xsum");
  zsum = memory->create(zsum, natom, "xsum");
  for (int i = 0; i < natom; ++i) xsum[i] = ysum[i] = zsum[i] = 0.;

  for (int i = 0; i < 5; ++i) fgets(str, MAXLINE, fp);
  for (int i = 0; i < ntm; ++i){
    fgets(str, MAXLINE, fp);

    vx = atof(strtok(str, " \t\n\r\f"));
    vy = atof(strtok(NULL," \t\n\r\f"));
    vz = atof(strtok(NULL," \t\n\r\f"));
    velx[nstep][i] = vx; xsum[i] += vx;
    if (sdim > 1) {vely[nstep][i] = vy; ysum[i] += vy;}
    if (sdim > 2) {velz[nstep][i] = vz; zsum[i] += vz;}

  }
  steps[nstep] = istep;
  ++nstep;

  // read the remaining 
  while ( 1 ){
    fgets(str, MAXLINE, fp);
    if (feof(fp)) break;
    fgets(str, MAXLINE, fp);
    istep = atoi(strtok(str," \t\n\r\f"));
    fgets(str, MAXLINE, fp);
    fgets(str, MAXLINE, fp);
    ntm = atoi(strtok(str," \t\n\r\f"));

    if (nstep >= nmax){
      nmax += 10000;
      velx = memory->grow(velx, nmax, natom, "ReadVel_velx");
      if (sdim > 1) vely = memory->grow(vely, nmax, natom, "ReadVel_vely");
      if (sdim > 2) velz = memory->grow(velz, nmax, natom, "ReadVel_velz");
      steps = memory->grow(steps, nmax, "ReadVel_steps");
    }

    for (int i = 0; i < 5; ++i) fgets(str, MAXLINE, fp);
    for (int i = 0; i < ntm; ++i){
      fgets(str, MAXLINE, fp);

      vx = atof(strtok(str, " \t\n\r\f"));
      vy = atof(strtok(NULL," \t\n\r\f"));
      vz = atof(strtok(NULL," \t\n\r\f"));
      velx[nstep][i] = vx; xsum[i] += vx;
      if (sdim > 1) {vely[nstep][i] = vy; ysum[i] += vy;}
      if (sdim > 2) {velz[nstep][i] = vz; zsum[i] += vz;}

    }
    if (istep < steps[nstep-1]) break;
    steps[nstep] = istep;
    ++nstep;
  }
  fclose(fp);

  // shrink memory if too much wasted
  if (nstep+1000 < nmax){
    velx = memory->grow(velx, nstep, natom, "ReadVel_velx");
    if (sdim > 1) vely = memory->grow(vely, nstep, natom, "ReadVel_vely");
    if (sdim > 2) velz = memory->grow(velz, nstep, natom, "ReadVel_velz");
    steps = memory->grow(steps, nstep, "ReadVel_steps");
  }
  for (int i = 0; i < natom; ++i) xsum[i] /= nstep;
  for (int is = 0; is < nstep; ++is)
  for (int i = 0; i < natom; ++i) velx[is][i] -= xsum[i];

  if (sdim > 1){
    for (int i = 0; i < natom; ++i) ysum[i] /= nstep;
    for (int is = 0; is < nstep; ++is)
    for (int i = 0; i < natom; ++i) vely[is][i] -= ysum[i];
  }

  if (sdim > 2){
    for (int i = 0; i < natom; ++i) zsum[i] /= nstep;
    for (int is = 0; is < nstep; ++is)
    for (int i = 0; i < natom; ++i) velz[is][i] -= zsum[i];
  }

  memory->destroy(xsum);
  memory->destroy(ysum);
  memory->destroy(zsum);

  printf("Done, %d timesteps for %d atoms read.\n\n", nstep, natom);
return;
}

/* -----------------------------------------------------------------------------
 * destructor, to free memory
 * -------------------------------------------------------------------------- */
ReadVel::~ReadVel()
{
  if (velx) memory->destroy(velx);
  if (vely) memory->destroy(vely);
  if (velz) memory->destroy(velz);
  if (steps) memory->destroy(steps);

  delete memory;
}

/* -------------------------------------------------------------------------- */
