# vacf

###  Introduction

`vacf` is a home written code to evaluate the velocity-velocity autocorrelation
function and then to compute the phonon density of states for an atomistic system
based on the velocities of atoms. The velocity information will be read from a file
in the LAMMPS dump format, which however can be generated from any molecular
dynamics simulation code, especially LAMMPS.

The code is written in C++ by Lingti Kong at Shanghai Jiao Tong University.

### Theoretical background

The autocorrelation function (ACF) reveals how the value of a quantity at a given
time/position correlates with its value at another time/position. Especially, for
the time autocorrelation *C(t)* of a quantity *v(t)*:

$$C(t) = \frac{\sum_{t_0} v(t + t_0) v(t_0)}{\sum_{t_0} v(t_0) v(t_0)}.$$

Where the summation goes over all time origins. In the case of velocity-velocity
autocorrelation, we generally average over all atoms (of the same element) as well:

$$C(t) = \frac{\sum_{i}\sum_{t_0} v_i(t + t_0) v_i(t_0)}{\sum_i\sum_{t_0} v_i(t_0) v_i(t_0)}.$$

The phonon density of states (DOS) can be given from the velocity autocorrelation function
by:

$$g(\nu) = \int_{-\infty}^{\infty}e^{2\pi i\nu} C(t)\mathrm{d}t.$$

With `vacf`, the $g(\nu)$ will be normalized to one:
$$\int g(\nu)\mathrm{d}\nu = 1.$$

### Installation

Preferrably, this code runs under Linux. A C++ compiler and the fftw library are
required to compile this code.

The [fftw](https://www.fftw.org) library can be accessed from [https://www.fftw.org](https://www.fftw.org).
Once installed, one can revise the `Makefile` to provid the locations of both the
header files (`include`) and library files (`lib`).

The generaly steps needed is:
- Install fftw; if it is already available, one can skip this step;
- Untar the `vacf.tar.gz` by: `tar -xzvf vacf.tar.gz`;
- Get into `vacf` by `cd vacf`;
- Revise the `Makefile` to define the compiler and the locations of [fftw](https://www.fftw.org);
- Compile the code by:

`make`

If everything goes succesfully, you should have an executable `vacf` now. You can copy it to
your home `bin` directory:

`cp vacf ~/bin`

### Usage
#### Obtaining the velocity dump file.
The velocity dump file can be obtained from molecular dynamics simulations by, e.g., LAMMPS.
Any other MD code also works, but the dump file must be converted into the LAMMPS dump format.

Generally, the velocity file will be huge, you'd better split your system into multiple groups,
and dump out the velocity for each group into a separate file. The dumped velocities for the atoms
must be in the same order for different frames. Especially, if your system has more than one kind
of elements, it is suggested to dump the velocities for atoms of different elements separately.
With LAMMPS, this can be achieved by:

```
group       C type 1

dump        1 C custom 5 dump.vc vx vy vz
dump_modify 1 sort id
```

which will dump the velocities of atoms in the group `C` (defined for atoms of `type 1`)
every `5` steps to the file `dump.vc`. `dump_modify 1 sort id` makes sure that the velocities
have been sorted according to the atom `id`s.

The frequency to dump the velocities (`5` in the above example) and the MD time step determines
the resolution of your autocorrelation function. Generally, it is recommended to be on the order
of `nevery * timestep = 5 fs.`

#### Calculating the velocity-velocity autocorrelation function and phonon DOS
The usage of the `vacf` code is:

`vacf [options] dump-files`

For example:
```
  vacf -dt 2e-3 -fr 0 20 -oa acf.dat -od dos.dat -df 0.01 -s 1 1 0.1 dump.vc
```

For explanation of the options, please refer to the section below. The above example
with try to compute the velocity-velocity autocorrelation function based on the velocity
information read from `dump.vc` and write to `acf.dat` (`-oa acf.dat`). The MD time
step is `2e-3 ps`, the DOS will be calculated from 0 to 20 THz (`-fr 0 20`) at a
stepsize of 0.01 (`-df 0.01`) and written to `dos.dat` (`-od dos.dat`).
Before calculating the phonon DOS, the ACF is smoothed by method `1` with a threshold
frequency of `1` with a width of `0.1` (`-s 1 1 0.1`).

The unit of the velocities does not matter here, as it will be cancelled out; while the
MD time step is expected to be in unit of `ps`.

#### Calculating the phonon DOS based on previously calculated velocity autocorrelation function
Sometimes one would like to calculate the phonon DOS based on a previously calculated
velocity autocorrelation function. In this case, one can use a command like the following to
achieve the goal:
```
   vacf -fr 0 20 -od dos.dat -df 0.01 -s 1 1 0.1 -r acf.dat
```

#### Help on the options
If one invoke `vacf -h`, the following help information will be provided:
```
  -oa vacf-file      : to define the output acf file name; default: acf.dat
  -od dos-file       : to define the output dos file name; default: dos_vacf.dat
  -dt timestep       : to define the MD time step in ps; default: 1e-3
  -l  lag-fraction   : to define the lag for acf as nlag=nsteps/lag-fraction; default: 100
  -s i tau0 dtau     : to define the smearing method and the parameters; default: not set
                       Suggested: tau0 = 5.*T, dtau = 0.5*T, with T the relaxation time.
  -fr vmin vmax      : to define the frequency range to output the PDOS; default: 0 to max. (THz)
  -df df             : to define the frequency stepsize to calculate/output the PDOS; default: 0.01 THz.
  -r                 : to indicate the file to read is vacf instead of velocities.
  -dim sysdim        : to define the system dimension, [1, 3]; default: 3
  -h                 : to print this help info.
```
Here some explanations on some of the options are provided:

#### `-dt`
This option defines the time step used in your molecular dynamics simulations in unit of `ps`.
If this option is not provided, a default time step of `1e-3 ps` is assumed.

#### `-l`
This option defines the maximum time lag used to calculate the autocorrelation function.
If defined, the time lag will be `nsteps / lag-fraction`, where `nstep` is the total number
of time intervals for the velocity data. If not defined, a default value of 100 will be
assumed.

Generally, the final time lag in terms of time (nsteps/lat-fraction * timestep) should be
an order greater than the decay time of the autocorrelation function, so as to guarantee
the reliability of the resultant phonon DOS.

#### `-s`
This option defines the method to smear the autocorrelation function before calculating the
phonon DOS. This is sometimes useful in that the ACF is accurate for small `t`, while inaccurate
for large `t`, as the number of data employed to calculate the ACF decreases with the increasing
of `t`. In turn, large fluctuations in ACF at large `t` might be observed. In order to reduce
such unexpected fluctuations, one can employ a smearing technique to smooth the ACF:
$$ C(t) = C0(t) * S(t).$$

The first parameter defines the method used, the second parameter defines the location where the
smearing factor $S(t)$ is 0.5, and the third parameter defines the width of the smearing.
The available smearing methods are ($ x = (t - \tau_0) / \Delta\tau$):

 - 1, Fermi-Dirac smearing:
$$ S(x) = \frac{1}{1 + \exp(x)}. $$

 - 2, Gaussian smearing:
$$ S(x) = \frac{1-\erf(x)}{2}. $$

 - 3, Power-4 smearing:
$$ S(x) = x^4/(1 + x^4)~~(x < 0). $$
$$ S(x) = 0.           ~~(x >= 0). $$

 - 4, Cubic-Fermi-Dirac smearing:
$$ S(x) = \frac{1}{1 + \exp(x^3)}. $$

 - 5, Cubic Gaussian smearing:
$$ S(x) = \frac{1-\erf(x^3)}{2}. $$

Here $\tau0$ should be in unit of time (`ps`). Generally, a value of several times the decay
time of the autocorrelation function is suggested. While the width parameter is recommned to
by one tenth of the decay time.

#### `-dim`
This defines the dimension of your system. For example, `-dim 1` suggests that your system
is one dimensionnal and the atoms move along the x direction only. The dump file is expected
to contain the `vx` data for all atoms only.

### Citation
For publications using this tool, the following paper is suggested to be cited:

Dynamical stability of iron under high-temperature and high-pressure conditions
L.T. Kong, J.F. Li, Q.W. Shi, H.J. Huang and K. Zhao
EPL 97 (2012) 56004 doi: 10.1209/0295-5075/97/56004

### Copyright
All copyrights are reserved to Prof. Lingti Kong ([:love_letter:](mailto:konglt@sjtu.edu.cn)).

Nov 2022
