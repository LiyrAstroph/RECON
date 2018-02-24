# RECON
A method for measuring spectral power and reconstructing time series in active galactic nuclei

reference: Yan-Rong Li & Jian-Min Wang, MNRAS Letter in press, [arXiv: 1802.07958](https://arxiv.org/abs/1802.07958)

# Compiling
To compile ``RECON``, the following third-party packages are required:
- **MPICH**, an MPI implementation (version 1.0 or higher), available at http://www-unix.mcs.anl.gov/mpi/mpich/
- **FFTW**, a fast Fourier transform library (version 3.0 or higher), available at http://www.fftw.org/
- **GSL**, the GNU Scientific Library (version 2.2.1 or higher), available at http://www.gnu.org/software/gsl
- **DNest**, the diffusive nested sampling library, developed by Yan-Rong Li, which is a C version of the DNest code by Brendon Brewer (https://github.com/eggplantbren/DNest3), available at https://github.com/LiyrAstroph/DNest_C

After installing the above packages, edit the corresponding library paths in ``Makefile``, and then use ``make`` to compile the ``RECON``.

# Running

```bash
mpiexec -n np ./recon param
```

where ``np`` is the number of cores and ``param`` is the parameter file.  

# Command-line Options
``RECON`` provides several command-line options:
- **-g**, generate a mock light curve using the PSD specified in ``param`` file, e.g.,
```bash
mpiexec -n 4 ./recon param -g
```
- **-s**, set a seed for random number generator, e.g., set a seed with a value of 100,
```bash
mpiexec -n 4 ./recon param -s100
```
- **-l**, implement level-dependent sampling, which usually improves sampling efficiency, e.g.,
```bash
mpiexec -n 4 ./recon param -l
```
- **-p**, only do posterior processing when MCMC sampling is available, e.g.,
```bash
mpiexec -n 4 ./recon param -p
```
- **-t**, set a temperature for posterior processing, e.g.,
```bash
mpiexec -n 4 ./recon param -t2
```
- **-d**, calculate the periodogram (PSD) of input data, e.g.,
```bash
mpiexec -n 4 ./recon param -d
```

Once can also combine the above options, e.g.,
```bash
mpiexec -n 4 ./recon param -ls100
```
or
```bash
mpiexec -n 4 ./recon param -pt2
```

# Parameter File
A typical parameter file looks like,
```
# parameter file
# lines beginning with '#' are regarded as comments and are neglected
#

PSDType    0                # psd type

#===============================================
# data

FileDir         .                # file direcotry

FileName        data/ngc5548.txt     # file name

V               2.0              # factor for extension of time interval in term of data
W               2.0              # factor for time resolution in terms of data  

FreqLimit       5.0e-4           # frequency limit below which PSD flattens off
FlagEndMatch    0                # end matching to reduce leakage
#===============================================
# only for simulation (-g option turned on)

FileSim         data/sim.txt       #
ND_Sim          200                # number of points
DT              1.0                # time resolution
FBad            0.33                # fraction of bad points
FErr            0.1                # measurement noise

# PSD argument
PSDArg          1.0e-3:2.5:0.0                # type 0, single power-law:  A:alpha:noise
#PSDArg          1.0e3:1.0e-2:0.0            # type 1, DRW: A:fknee:noise                         
#PSDArg           1.0e1:2.5:1.5:1.5e-2:0.01   # type 2, bending power-law: A:a1:a2:freq:noise
#===============================================
```

# Outputs 
``RECON`` generates posterior samples ``posterior_sample.txt`` for PSD parameters and frequency series for reconstructing input time series. The reconstructions based on posterior samples are stored in ``recon.txt``.

``RECON`` also calculates evidence at the end of running, printed out on terminal screen.

# An Exemplary Case
Application to the 5100A light curve of NGC 5548, see arXiv: 1802.07958.
![Application to the 5100A light curve of NGC 5548](https://github.com/liyropt/MyGithubPic/blob/master/ngc5548_lc_github.jpg)

**A more detailed usage guideline is coming.**
