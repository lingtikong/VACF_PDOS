#include "vacf.h"
#include "string.h"
#include "timer.h"
//#include <gsl/gsl_spline.h>

#define MAXLINE 256
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/* -----------------------------------------------------------------------------
 * constructor, to analyze command line option and to drive the code
 * -------------------------------------------------------------------------- */
VACF::VACF(int narg, char **arg)
{
  if (narg < 2) help();

  dt = 1.e-3;
  sysdim = 3;
  nlag = ncount = ndos = 0;
  memory = new Memory();
  sum = NULL; pdos = NULL;
  outacf = outdos = NULL;
  dlag = 100;
  npad = 10;

  Timer *timer = new Timer();

  int flag_read_acf = 0;
  vmin = vmax = 0.;

  int iarg = 1;
  while (iarg < narg){
    if (strcmp(arg[iarg],"-oa") == 0) {
      if (++iarg >= narg) help();

      if (outacf) delete []outacf;
      outacf = new char[strlen(arg[iarg])+1];
      strcpy(outacf, arg[iarg]);
    
    } else if (strcmp(arg[iarg],"-od") == 0) {
      if (++iarg >= narg) help();

      if (outdos) delete []outdos;
      outdos = new char[strlen(arg[iarg])+1];
      strcpy(outdos, arg[iarg]);
    
    } else if ( strcmp(arg[iarg], "-dt") == 0) {
      if (++iarg >= narg) help();
      dt = atof(arg[iarg]);

    } else if ( strcmp(arg[iarg], "-l") == 0) {
      if (++iarg >= narg) help();
      dlag = atoi(arg[iarg]);
      if (dlag < 1) dlag = 100;

    } else if ( strcmp(arg[iarg], "-dim") == 0) {
      if (++iarg >= narg) help();
      sysdim = atoi(arg[iarg]);
      if (sysdim<1 || sysdim>3) help();

    } else if ( strcmp(arg[iarg], "-fr") == 0) {
      if (iarg+2 >= narg) help();
      vmin  = atof(arg[++iarg]);
      vmax  = atof(arg[++iarg]);

      if (vmin < 0.)    vmin = 0.;
      if (vmax <= vmin) vmax = 0.;

    } else if ( strcmp(arg[iarg], "-p") == 0) {
      if (++iarg >= narg) help();
      npad = atoi(arg[iarg]);
      if (npad < 1 ) npad = 10;

    } else if ( strcmp(arg[iarg], "-h") == 0) {
      help();

    } else if ( strcmp(arg[iarg], "-r") == 0) {
      flag_read_acf = 1;

    } else {
      timer->start();
      if (flag_read_acf){
        read_acf(arg[iarg]);

      } else {
        vels = new ReadVel(arg[iarg], sysdim);
        compute_acf();

        dstep = (vels->steps[nlag] - vels->steps[0])/nlag;
        delete vels; vels = NULL;
      }
      timer->stop(); timer->print();
    }
    iarg++;
  }
  if (ncount < 1) return;

  if (!flag_read_acf) write_acf();

  timer->start();
  printf("\nNow to compute the pdos..."); fflush(stdout);
  compute_dos();
  write_dos();
  printf("Done!");
  timer->stop(); timer->print();
 
return;
}

/* -----------------------------------------------------------------------------
 * destructor, to free memory
 * -------------------------------------------------------------------------- */
VACF::~VACF()
{
  if (vels) delete vels;
  if (sum)  memory->destroy(sum);
  if (pdos) memory->destroy(pdos);
  if (outacf) delete[]outacf;
  if (outdos) delete[]outdos;

  if (memory) delete memory;
}

/* -----------------------------------------------------------------------------
 * private method to write out the computed acf
 * -------------------------------------------------------------------------- */
void VACF::write_acf()
{

  if (outacf == NULL){
    outacf = new char[8];
    strcpy(outacf,"acf.dat");
  }

  printf("\nNow to write the acf data.\n\n");
  FILE *fp=fopen(outacf, "w");
  if (sysdim == 1){
    fprintf(fp,"# lag  ave-vacf\n");
    fprintf(fp,"# ps    ------ \n");
    for (int i=0; i<nlag; i++) fprintf(fp,"%lg %lg\n", i*dstep*dt, sum[i][0]);

  } else if (sysdim == 2){
    fprintf(fp,"# lag  vacf-ave vacf-x  vacf-y\n");
    fprintf(fp,"# ps    ----     ----    ----  \n");
    for (int i=0; i<nlag; i++) fprintf(fp,"%lg %lg %lg %lg\n", i*dstep*dt, (sum[i][0]+sum[i][1])*0.5, sum[i][0], sum[i][1]);

  } else {
    fprintf(fp,"# lag  vacf-ave vacf-x  vacf-y  vacf-z\n");
    fprintf(fp,"# ps    ----     ----    ----    ---- \n");
    for (int i=0; i<nlag; i++) fprintf(fp,"%lg %lg %lg %lg %lg\n", i*dstep*dt, (sum[i][0]+sum[i][1]+sum[i][2])/3., sum[i][0], sum[i][1], sum[i][2]);

  }

  fclose(fp);

}

/* -----------------------------------------------------------------------------
 * private method to compute the acf, using fftw
 * -------------------------------------------------------------------------- */
void VACF::compute_acf()
{
  int ndat = vels->nstep;
  if (ndat < 1) {
    printf("\nError: No data read!\n\n");
    return;
  }

  int nttt = ndat/dlag; // set nlag as 1/10 of the total MD dumps
  if (nlag == 0){
    nlag = nttt;
    if (nlag < 10){
      printf("\nError: your velocity file contains too few data! At least %d is required!\n", 10*dlag);
      exit(2);
    }
    if (sum)  memory->destroy(sum);
    sum = memory->create(sum, nlag, sysdim, "sum");
    for (int i=0; i<nlag; i++)
    for (int j=0; j<sysdim; j++) sum[i][j] = 0.;
  } else if (ndat < nlag) {
    printf("Error: Too few data points read from this file!\n");
    return;
  }
  int ntotal = ndat + nlag; 

  printf("Now to compute the acf... "); fflush(stdout);

  // prepare fftw
  double *fftw_real;
  fftw_complex *fftw_cmpx;

  fftw_real = memory->create(fftw_real, ntotal, "fftw_real");
  fftw_cmpx = memory->create(fftw_cmpx, ntotal/2+1, "fftw_cmpx");

  fftw_plan r2c = fftw_plan_dft_r2c_1d(ntotal, fftw_real, fftw_cmpx, FFTW_MEASURE);
  fftw_plan c2r = fftw_plan_dft_c2r_1d(ntotal, fftw_cmpx, fftw_real, FFTW_MEASURE);

  double fac = 1./sqrt(double(ndat));

  // acf via fft
  for (int itm=0; itm < vels->natom; itm++){
    for (int i=0; i<ndat; i++) fftw_real[i] = vels->velx[i][itm] * fac;
    for (int i=ndat; i<ntotal; i++) fftw_real[i] = 0.;

    fftw_execute(r2c);

    for (int i=0; i<ntotal/2+1; i++) fftw_cmpx[i] *= conj(fftw_cmpx[i]);

    fftw_execute(c2r);

    for (int i=0; i<nlag; i++) sum[i][0] += fftw_real[i]/double(ndat-i);

    ncount++;
  }

  if (sysdim > 1){
    for (int itm=0; itm < vels->natom; itm++){
      for (int i=0; i<ndat; i++) fftw_real[i] = vels->vely[i][itm] * fac;
      for (int i=ndat; i<ntotal; i++) fftw_real[i] = 0.;
  
      fftw_execute(r2c);
  
      for (int i=0; i<ntotal/2+1; i++) fftw_cmpx[i] *= conj(fftw_cmpx[i]);
  
      fftw_execute(c2r);
  
      for (int i=0; i<nlag; i++) sum[i][1] += fftw_real[i]/double(ndat-i);
  
      ncount++;
    }
  }

  if (sysdim > 2){
    for (int itm=0; itm < vels->natom; itm++){
      for (int i=0; i<ndat; i++) fftw_real[i] = vels->velz[i][itm]*fac;
      for (int i=ndat; i<ntotal; i++) fftw_real[i] = 0.;
  
      fftw_execute(r2c);
  
      for (int i=0; i<ntotal/2+1; i++) fftw_cmpx[i] *= conj(fftw_cmpx[i]);
  
      fftw_execute(c2r);
  
      for (int i=0; i<nlag; i++) sum[i][2] += fftw_real[i]/double(ndat-i);
  
      ncount++;
    }
  }

  memory->destroy(fftw_real);
  memory->destroy(fftw_cmpx);


  // normalize acf
/*
  double ave[3]; ave[0] = ave[1] = ave[2] = 0.;
  for (int i=0; i<nlag; i++)
    for (int j=0; j<sysdim; j++) ave[j] += sum[i][j];
  for (int j=0; j<sysdim; j++) ave[j] /= nlag;

  //double icr0 = 1./(sum[0] - ave);
  for (int i=1; i<nlag; i++)
    for (int j=0; j<sysdim; j++) sum[i][j] = (sum[i][j] - ave[j])/(sum[0][j]-ave[j]);
  for (int j=0; j<sysdim; j++) sum[0][j] = 1.;
*/
  for (int i=1; i<nlag; ++i)
    for (int j=0; j<sysdim; j++) sum[i][j] /= sum[0][j];
  for (int j=0; j<sysdim; ++j) sum[0][j] = 1.;

  printf("Done!\n"); fflush(stdout);
return;
}

/* -----------------------------------------------------------------------------
 * private method to compute the phonon dos, by fft the velocity-velocity acf
 * -------------------------------------------------------------------------- */
void VACF::compute_dos()
{
  double *fftw_real;
  fftw_complex *fftw_cmpx;

  int ntotal = nlag * npad;

  dv = 1./double(dstep*ntotal*dt);
  ndos = ntotal/2;
  if (vmax > vmin) ndos = MIN(int(vmax/dv)+1, ndos); 

  if (pdos) memory->destroy(pdos);
  pdos = memory->create(pdos, ndos, sysdim, "pdos");

  fftw_real = memory->create(fftw_real, ntotal, "fftw_real");
  fftw_cmpx = memory->create(fftw_cmpx, ntotal/2+1, "fftw_cmpx");

  fftw_plan r2c = fftw_plan_dft_r2c_1d(ntotal, fftw_real, fftw_cmpx, FFTW_ESTIMATE);

  for (int idim=0; idim<sysdim; idim++){

    for (int i=0; i<nlag; i++) fftw_real[i] = sum[i][idim];
    for (int i=nlag; i<ntotal; i++) fftw_real[i] = 0.;

    // FFT
    fftw_execute(r2c);

    for (int i=0; i<ndos; i++) pdos[i][idim] = creal(fftw_cmpx[i]) - creal(fftw_cmpx[0]);
  }

return;
}
  
/* -----------------------------------------------------------------------------
 * private method to write the computed pdos to file
 * -------------------------------------------------------------------------- */
void VACF::write_dos()
{
  // to output the computed dos
  if (outdos == NULL){
    outdos = new char[13];
    strcpy(outdos,"dos_vacf.dat");
  }
  FILE *fp = fopen(outdos, "w");
  int istr, iend;
  if (vmin <= 0.) istr = 0;
  else istr = int(vmin/dv);

  if (vmax <= vmin) iend = ndos;
  else iend = int(vmax/dv);
  iend = MIN(iend, ndos);

  if (sysdim == 1){
    fprintf(fp,"#freq   dos\n");
    fprintf(fp,"#THz    ---\n");
    for (int i=istr; i<iend; i++) fprintf(fp,"%lg %lg\n", i*dv, pdos[i][0]);

  } else if (sysdim == 2){
    fprintf(fp,"#freq dos-total dos-x dos-y\n");
    fprintf(fp,"#THz   -------   ---   ---\n");
    for (int i=istr; i<iend; i++) fprintf(fp,"%lg %lg %lg %lg\n", i*dv, (pdos[i][0] + pdos[i][1])*0.5, pdos[i][0], pdos[i][1]);

  } else {
    fprintf(fp,"#freq dos-total dos-x dos-y dos-z\n");
    fprintf(fp,"#THz   -------   ---   ---   ---\n");
    for (int i=istr; i<iend; i++)
      fprintf(fp,"%lg %lg %lg %lg %lg\n", i*dv, (pdos[i][0] + pdos[i][1] + pdos[i][2])/3., pdos[i][0], pdos[i][1], pdos[i][2]);
  }

  fclose(fp);

return;
}

/* -----------------------------------------------------------------------------
 * private method to read acf info from file
 * The acf should be in the same format as the one written by this code, and
 * its first column must be time lag in unit of fs, begin from 0.
 * -------------------------------------------------------------------------- */
void VACF::read_acf(char *acf_file)
{
  printf("\nNow to read acf info from file: %s ...", acf_file);
  FILE *fp = fopen(acf_file, "r");
  int nmax = 10000;
  nlag = 0;
  sum = memory->grow(sum, nmax, "sum");

  double t0;
  char oneline[MAXLINE];
  fgets(oneline, MAXLINE, fp);
  fgets(oneline, MAXLINE, fp);
  fgets(oneline, MAXLINE, fp);
  sysdim = count_words(oneline);
  if (sysdim < 2) {
    printf("Error: not enough columns in file %s!\n", acf_file);
    exit(1);
  } else if (sysdim > 2) sysdim -= 2;
  else sysdim--;
  sum = memory->grow(sum, nmax, sysdim, "sum");

  while (! feof(fp) ){
    if (nlag >= nmax){
      nmax += 10000;
      sum = memory->grow(sum, nmax, sysdim, "sum");
    }

    t0 = atof(strtok(oneline," \r\t\n\f"));
    if (sysdim > 1) strtok(NULL," \r\t\n\f");
    for (int i=0; i<sysdim; i++) sum[nlag][i] = atof(strtok(NULL," \r\t\n\f"));
    nlag++;

    fgets(oneline, MAXLINE, fp);
  }
  fclose(fp);
  dstep = 1;
  dt = t0/double(nlag-1);

  printf("Done! %d data read.\n", nlag);
  ncount = 1;

return;
}

/* -----------------------------------------------------------------------------
 * private method to print help info
 * -------------------------------------------------------------------------- */
void VACF::help()
{
  printf("\nvacf\n\nCode to compute the velocity-velocity autocorrelation based on\n");
  printf("velocities dumped from LAMMPS. The dump command should look like:\n");
  printf("  dump         2 all custom N file vx vy vz\n");
  printf("  dump_modify  2 sort id\n\n");
  printf("Usage:\n  vacf [options] lmp-dump-file\n\n");
  printf("Available options:\n");
  printf("  -oa vacf-file      : to define the output acf file name; default: acf.dat\n");
  printf("  -od dos-file       : to define the output dos file name; default: dos_vacf.dat\n");
  printf("  -dt timestep       : to define the MD time step in ps; default: 1e-3\n");
  printf("  -l  lag-fraction   : to define the lag for acf as nlag=nsteps/lag-fraction; default: 100\n");
  printf("  -p  padding-acf    : to define the padding of acf as nlag * padding-acf, for dos; default: 10\n");
  printf("  -fr vmin vmax      : to define the frequency range to output the PDOS; default: 0 to max. (THz)\n");
  printf("  -r                 : to indicate the file to read is vacf instead of velocities.\n");
  printf("  -dim sysdim        : to define the system dimension, [1, 3]; default: 3\n");
  printf("  -h                 : to print this help info.\n\n\n");

  exit(0);
}
/* -------------------------------------------------------------------------- */

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
int VACF::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy = memory->create(copy,n,"copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->destroy(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->destroy(copy);
  return n;
}
