#include "vacf.h"
#include "string.h"
#include "timer.h"

#define ZERO 1.e-8
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
  ismear = 0;
  nlag = ncount = ndos = 0;
  memory = new Memory();
  vv0 = NULL;
  sum = pdos = NULL;
  outacf = outdos = NULL;
  dlag = 10;

  Timer *timer = new Timer();

  int flag_read_acf = 0;
  vmin = vmax = 0.;
  vstep = 0.01; // Default frequency step size, in unit of THz.

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
      if (dlag < 1) dlag = 10;

    } else if ( strcmp(arg[iarg], "-s") == 0) {
      if (iarg+3 >= narg) help();
      ismear = atoi(arg[++iarg]);
      tau0 = atof(arg[++iarg]);
      dtau = atof(arg[++iarg]);

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

    } else if ( strcmp(arg[iarg], "-df") == 0) {
      if (++iarg >= narg) help();
      vstep = atof(arg[iarg]);
      if (vstep <= 0.) help();

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
        compute_acf(arg[iarg]);

        dstep = (vels->steps[nlag] - vels->steps[0])/nlag;
        delete vels; vels = NULL;
      }
      timer->stop(); timer->print();
    }
    ++iarg;
  }
  if (ncount < 1) return;

  normal_acf();
  if (!flag_read_acf) write_acf();

  // compute phDOS
  timer->start();
  printf("\nNow to compute the pdos..."); fflush(stdout);
  compute_dos();
  normal_dos();
  write_dos();

  printf("Done! ");
  timer->stop(); timer->print();
  delete timer;
 
return;
}

/* -----------------------------------------------------------------------------
 * destructor, to free memory
 * -------------------------------------------------------------------------- */
VACF::~VACF()
{
  memory->destroy(sum);
  memory->destroy(pdos);
  if (outacf) delete []outacf;
  if (outdos) delete []outdos;
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
    fprintf(fp,"# <v*v>(0): %lg\n", vv0[0]);
    fprintf(fp,"# lag  vacf-x\n");
    fprintf(fp,"# ps    ------ \n");
    for (int i = 0; i < nlag; ++i) fprintf(fp,"%lg %lg %lg\n", i*dstep*dt, sum[i][0], sum[i][0]*smearing(i));

  } else if (sysdim == 2){
    fprintf(fp,"# <v*v>(0): %lg %lg\n", vv0[0], vv0[1]);
    fprintf(fp,"# lag  vacf-x  vacf-y vacf-ave\n");
    fprintf(fp,"# ps    ----     ----    ----  \n");
    for (int i = 0; i < nlag; ++i) fprintf(fp,"%lg %lg %lg %lg\n", i*dstep*dt, sum[i][0], sum[i][1], (sum[i][0]+sum[i][1])*0.5*smearing(i));

  } else {
    fprintf(fp,"# <v*v>(0): %lg %lg %lg\n", vv0[0], vv0[1], vv0[2]);
    fprintf(fp,"# lag  vacf-x  vacf-y  vacf-z vacf-ave\n");
    fprintf(fp,"# ps    ----     ----    ----    ---- \n");
    for (int i = 0; i < nlag; ++i) fprintf(fp,"%lg %lg %lg %lg %lg\n", i*dstep*dt, sum[i][0], sum[i][1], sum[i][2], (sum[i][0]+sum[i][1]+sum[i][2])/3.*smearing(i));

  }

  fclose(fp);
  return;
}

/* -----------------------------------------------------------------------------
 * private method to compute the acf, using fftw
 * -------------------------------------------------------------------------- */
void VACF::compute_acf(char *fname)
{
  int ndat = vels->nstep;
  if (ndat < 1) {
    printf("\nError: No data read!\n\n");
    return;
  }

  int nttt = ndat / dlag; // set nlag as 1/10 of the total MD dumps
  if (nlag == 0){
    nlag = nttt;
    if (nlag < 10){
      printf("\nError: your velocity file contains too few data! At least %d is required!\n", 10*dlag);
      exit(2);
    }
    if (sum)  memory->destroy(sum);
    sum = memory->create(sum, nlag, sysdim, "sum");
    for (int i = 0; i < nlag; ++i)
    for (int j = 0; j < sysdim; ++j) sum[i][j] = 0.;

  } else if (ndat < nlag) {
    printf("Error: Too few data points read from this file!\n");
    return;
  }
  int ntotal = ndat + nlag; 

  printf("Now to compute acf for data read from %s... ", fname); fflush(stdout);

  // prepare fftw
  double *fftw_real;
  fftw_complex *fftw_cmpx;

  fftw_real = memory->create(fftw_real, ntotal, "fftw_real");
  fftw_cmpx = memory->create(fftw_cmpx, ntotal/2+2, "fftw_cmpx");

  fftw_plan r2c = fftw_plan_dft_r2c_1d(ntotal, fftw_real, fftw_cmpx, FFTW_MEASURE);
  fftw_plan c2r = fftw_plan_dft_c2r_1d(ntotal, fftw_cmpx, fftw_real, FFTW_MEASURE);

  double fac = 1./sqrt(double(ndat));

  // acf via fft
  for (int itm = 0; itm < vels->natom; ++itm){
    for (int i = 0; i < ndat; ++i) fftw_real[i] = vels->velx[i][itm] * fac;
    for (int i = ndat; i < ntotal; ++i) fftw_real[i] = 0.;

    fftw_execute(r2c);

    for (int i = 0; i < ntotal/2+1; ++i) fftw_cmpx[i] *= conj(fftw_cmpx[i]);

    fftw_execute(c2r);

    for (int i = 0; i < nlag; ++i) sum[i][0] += fftw_real[i] / double(ndat-i);

    ++ncount;
  }

  if (sysdim > 1){
    for (int itm = 0; itm < vels->natom; ++itm){
      for (int i = 0; i < ndat; ++i) fftw_real[i] = vels->vely[i][itm] * fac;
      for (int i = ndat; i < ntotal; ++i) fftw_real[i] = 0.;
  
      fftw_execute(r2c);
  
      for (int i = 0; i < ntotal/2+1; ++i) fftw_cmpx[i] *= conj(fftw_cmpx[i]);
  
      fftw_execute(c2r);
  
      for (int i = 0; i < nlag; ++i) sum[i][1] += fftw_real[i] / double(ndat-i);
  
      ++ncount;
    }
  }

  if (sysdim > 2){
    for (int itm = 0; itm < vels->natom; ++itm){
      for (int i = 0; i < ndat; ++i) fftw_real[i] = vels->velz[i][itm]*fac;
      for (int i = ndat; i < ntotal; ++i) fftw_real[i] = 0.;
  
      fftw_execute(r2c);
  
      for (int i = 0; i < ntotal/2+1; ++i) fftw_cmpx[i] *= conj(fftw_cmpx[i]);
  
      fftw_execute(c2r);
  
      for (int i = 0; i < nlag; ++i) sum[i][2] += fftw_real[i] / double(ndat-i);
  
      ++ncount;
    }
  }

  memory->destroy(fftw_real);
  memory->destroy(fftw_cmpx);
  fftw_destroy_plan(r2c);
  fftw_destroy_plan(c2r);

  printf("Done!\n"); fflush(stdout);
return;
}

/* -----------------------------------------------------------------------------
 * private method to compute the phonon dos, by fft the velocity-velocity acf
 * -------------------------------------------------------------------------- */
void VACF::compute_dos()
{
  double *fftw_in, *fftw_out;

  int ntotal = int(1./(2.*dstep*dt*vstep)) + 1;
  if (ntotal < nlag) ntotal = nlag;

  dv = 1./double(dstep*2*(ntotal-1)*dt);
  if (vmax <= ZERO){
    ndos = ntotal / 3;

  } else {
    ndos = MIN(int(vmax/dv)+1, ntotal);

  }

  if (pdos) memory->destroy(pdos);
  pdos = memory->create(pdos, ndos, sysdim, "pdos");

  fftw_in  = memory->create(fftw_in,  ntotal, "fftw_in");
  fftw_out = memory->create(fftw_out, ntotal, "fftw_out");

  fftw_plan r2r = fftw_plan_r2r_1d(ntotal, fftw_in, fftw_out, FFTW_REDFT00, FFTW_ESTIMATE);

  for (int idim = 0; idim < sysdim; ++idim){

    for (int i = 0; i < nlag; ++i) fftw_in[i] = sum[i][idim];
    for (int i = nlag; i < ntotal; ++i) fftw_in[i] = 0.;

    // FFT
    fftw_execute(r2r);

    for (int i = 0; i < ndos; ++i) pdos[i][idim] = fftw_out[i];
  }

  memory->destroy(fftw_in);
  memory->destroy(fftw_out);
  fftw_destroy_plan(r2r);

  // calculate diffusion coefficients; if metal unit, D will be in cm^2/s.
  // The expression here might be wrong.
  for (int j = 0; j < sysdim; ++j) vv0[j] = pdos[0][j]*vv0[j]*dstep*dt*1.e-4 / double(ntotal);

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
  else iend = int(vmax/dv) + 1;
  iend = MIN(iend, ndos);

  if (sysdim == 1){
    // if (vv0) fprintf(fp,"#Diffusion coefficients: %lg cm^2/s\n", vv0[0]);
    fprintf(fp,"#freq   dos\n");
    fprintf(fp,"#THz    ---\n");
    for (int i = istr; i < iend; ++i) fprintf(fp,"%lg %lg\n", i*dv, pdos[i][0]);

  } else if (sysdim == 2){
    // if (vv0) fprintf(fp,"#Diffusion coefficients: %lg, %lg, %lg (ave) cm^2/s\n", vv0[0], vv0[1], (vv0[0]+vv0[1])*0.5);
    fprintf(fp,"#freq dos-x dos-y dos-total\n");
    fprintf(fp,"#THz  ----- ----- ---------\n");
    for (int i = istr; i < iend; ++i) fprintf(fp,"%lg %lg %lg %lg\n", i*dv, pdos[i][0], pdos[i][1], (pdos[i][0] + pdos[i][1])*0.5);

  } else {
    // if (vv0) fprintf(fp,"#Diffusion coefficients: %lg, %lg, %lg, %lg (ave) cm^2/s\n", vv0[0], vv0[1], vv0[2], (vv0[0]+vv0[1]+vv0[2])/3.);
    fprintf(fp,"#freq dos-x dos-y dos-z dos-total\n");
    fprintf(fp,"#THz  ----- ----- ----- ---------\n");
    for (int i = istr; i < iend; ++i)
      fprintf(fp,"%lg %lg %lg %lg %lg\n", i*dv, pdos[i][0], pdos[i][1], pdos[i][2], (pdos[i][0] + pdos[i][1] + pdos[i][2])/3.);
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

  double t0;
  char oneline[MAXLINE];
  fgets(oneline, MAXLINE, fp);
  if (vv0 == NULL){
    vv0 = memory->create(vv0, 3, "vv0");
    sprintf(oneline, "%s 0. 0.\n", oneline);
    char *ptr = strtok(oneline," \r\t\n\f");
    ptr = strtok(NULL, " \r\t\n\f");
    for (int j = 0; j < 3; ++j) vv0[j] = atof(strtok(NULL, " \r\t\n\f"));
  }
  fgets(oneline, MAXLINE, fp);
  fgets(oneline, MAXLINE, fp);

  fgets(oneline, MAXLINE, fp);
  sysdim = count_words(oneline);
  if (sysdim < 2) {
    printf("Error: not enough columns in file %s!\n", acf_file);
    exit(1);
  } else if (sysdim > 2) sysdim -= 2;
  else --sysdim;

  sum = memory->create(sum, nmax, sysdim, "sum");

  while (! feof(fp) ){
    if (nlag >= nmax){
      nmax += 10000;
      sum = memory->grow(sum, nmax, sysdim, "sum");
    }

    t0 = atof(strtok(oneline," \r\t\n\f"));
    for (int i = 0; i < sysdim; ++i) sum[nlag][i] = atof(strtok(NULL," \r\t\n\f"));
    ++nlag;

    fgets(oneline, MAXLINE, fp);
  }
  fclose(fp);

  dstep = 1;
  dt = t0 / double(nlag-1);

  printf("Done! %d data read.\n", nlag);
  ncount = 1;

return;
}

/* -----------------------------------------------------------------------------
 * private method to bring the end part of the acf to zero.
 * 1: Fermi-Dirac smearing;
 * 2: Gaussian smearing.
 * -------------------------------------------------------------------------- */
double VACF::smearing(int istep)
{
  if (dtau <= ZERO) dtau = 0.2;
  if (tau0 <= ZERO) tau0 = 2.0;

  double x = istep * dstep * dt;
  double dx = (x - tau0) / dtau;

  if (ismear == 1){
     return 1./(1. + exp(dx));

  } else if (ismear == 2){
     return 0.5*(1.-erf(dx));

  } else if (ismear == 3){
     if (dx < 0.) {
        double x2 = dx * dx;
        double x4 = x2 * x2;
        return x4 / (1. + x4);
     } else
        return 0.;

  } else if (ismear == 4){
     return 1./(1. + exp(dx*dx*dx));

  } else if (ismear == 5){
     return 0.5*(1.-erf(dx*dx*dx));

  }
  return 1.;
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
  printf("  -l  lag-fraction   : to define the lag for acf as nlag = nsteps/lag-fraction; default: 100\n");
  printf("  -s i tau0 dtau     : to define the smearing method and the parameters; default: not set\n");
  printf("                       For x = (t - tau0)/dtau, the available smearing functions are:\n");
  printf("                        i = 1: S(x) = 1 / ( 1 + exp(x) )    (Fermi-Dirac); \n");
  printf("                        i = 2: S(x) = (1 - erf(x))/2        (Gaussian); \n");
  printf("                        i = 3: S(x) = x**4 / (1 + x**4)     (Power 4); \n");
  printf("                        i = 4: S(x) = 1 / ( 1 + exp(x**3) ) (Cubic Fermi-Dirac); \n");
  printf("                        i = 5: S(x) = (1 - erf(x**3))/2     (Cubic Gaussian); \n");
  printf("                       tau0 should be on order of 2-5 times of T, dtau roughly one tenth of T,\n");
  printf("                       with T the relaxation time.\n");
  printf("  -fr vmin vmax      : to define the frequency range to output the PDOS; default: 0 to max. (THz)\n");
  printf("  -df df             : to define the frequency stepsize to calculate/output the PDOS; default: 0.01 THz.\n");
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

/* ----------------------------------------------------------------------------
 * Private method to normalize the acf; smearing method will be used to bring
 * long time term to zero, if defined.
 * ---------------------------------------------------------------------------- */
void VACF::normal_acf()
{
  if (vv0 == NULL){
    vv0 = memory->create(vv0, sysdim, "vv");
    for (int j = 0; j < sysdim; ++j) vv0[j] = sum[0][j];
  }

  for (int j = 0; j < sysdim; ++j){
    for (int i = 1; i < nlag;   ++i) sum[i][j] /= sum[0][j];
    sum[0][j] = 1.;
  }

  // apply smearing
  if (ismear > 0){
    for (int i = 0; i < nlag; ++i){
      double scale = smearing(i);
      for (int j = 0; j < sysdim; ++j) sum[i][j] *= scale;
    }
  }

  return;
}

/* ----------------------------------------------------------------------------
 * Private method to normalize the DOS; Simpson's rule is used for the integration.
 * ---------------------------------------------------------------------------- */
void VACF::normal_dos()
{
  double odd, even, psum;

  for (int idim = 0; idim < sysdim; ++idim){
    odd = even = psum = 0.;
    for (int i = 1; i < ndos-1; i +=2) odd  += pdos[i][idim];
    for (int i = 2; i < ndos-1; i +=2) even += pdos[i][idim];
    psum = pdos[0][idim] + pdos[ndos-1][idim];
    psum += 4.*odd + 2.*even;
    psum = 3./(psum*dv);
    for (int i = 0; i < ndos; ++i) pdos[i][idim] *= psum;
  }
 
  return;
}
