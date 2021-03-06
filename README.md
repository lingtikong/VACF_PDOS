# vacf

###  Introduction

`vacf` is a home written code to evaluate the velocity-velocity autocorrelation
function and then to compute the phonon density of states for an atomistic system
based on the velocities of atoms. The velocity information will be read from a file
in the LAMMPS dump format, which however can be generated from any molecular
dynamics simulation code, especially LAMMPS.

### Theoretical background

The autocorrelation function (ACF) reveals how the value of a quantity at a given
time/position correlates with its value at another time/position. Especially, for
the time autocorrelation *C(t)* of a quantity *v(t)*:
![formula](https://render.githubusercontent.com/render/math?math=
C(t) = \frac{\sum_{t_0} v(t + t_0) v(t_0)}{\sum_{t_0} v_(t_0) v(t_0)}.)
Where the summation goes over all time origins.

### Installation

### Usage

vacf vel-dump-file-name(s)

The computed velocity-velocity autocorrelation will be written to
acf.dat, unless one specify a file name via the "-oa acf-file-name"
command line option. While the computed phonon density of states
will be written to file "dos_vacf.dat", unless one specify a different
file name via "-od dos-file-name".

The velocities should be dumped via:
dump  2  all custom n  dump-file vx vy vz
dump_modify 2 sort id

You can decompose your system into several groups and dump each
group speparately, so as to avoid super-huge files.

Note, the "dump_modify sort id" is essential to guarantee the
correctness of the acf and dos; it is also important to provide
the correct MD time step via command line option "-dt" to make
sure the units of time (in acf file) and frequency (in dos file)
to be correct. The desired unit for timestep is ps, in which
case the corresponding frequency unit will be THz. For unit
other than ps, the unit for frequency will by inverse of your
time unit.
### Copyright
