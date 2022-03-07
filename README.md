# nmr
# NMR from Molecular Dynamics!

This code is part of the codes wrote by Sylvia to get NMR properties from LAMMPS dump files.

# Files

### main files
- main.f90;
- routines.f90;
- t2calc.f90;

### example file
- corr-acid.in : input example file;
- NMR-T2 : detailed pdf about how code works;

## Used Equations

Selfcorrelation function for spin I = $\frac{1}{2}$, $\Omega(\theta,\phi)$ (with $\theta$ being the polar angle and $\phi$ azimutal one), m=0, $Y_{2,0} = \sqrt{\frac{5}{16\pi}}(3cos^2\theta-1)$ : 

$$G_{R,T}^{m} \implies G_{R,T}^{0}(t) = \frac{3\mu_0^2\hbar^2\gamma^4}{256\pi^2N}\sum_{i\neq j}^{N}\left(\frac{9(cos^2\theta_{ij}(t+\tau)-1)(cos^2\theta_{ij}(\tau)-1)}{r^3_{ij}(t+\tau)r_{ij}^3(\tau)}\right)$$

The translational and rotations proton time $T_{2,T,R}$ can be determined as

$$\frac{1}{T_{2,R,T}} = \frac{3}{2}J_{R,T}(0)+\frac{5}{2}J_{R,T}(\omega_0)+J_{R,T}(2\omega_0)$$

where $\omega_0$ is the Larmor frequency and $J_{R,T} = \mathcal{F}(G_{R,T}(t))$ (the fourier transform need to be used with foward propagation from [0,2$\pi$))

$$J_{R,T}(\omega) = 2\int_{R^+}G_{R,T}cos(\omega t)dt $$ 

since we have a fast regime, $\omega \lt \lt 1$ the $T_2$ time will be give by

$$\frac{1}{T_2} =\frac{10(\Delta \omega^2_R\tau_R + \Delta \omega^2_T\tau_T)}{9}$$


## How to proceed ?

In summary, the following are the steps to calculate T 2 relaxation time from MD
simulations for bulk liquids:

- i.
Save the position trajectories (x,y,z) of the hydrogen atoms from a molecular
dynamics simulation e.g lammpstrj file from a MD simulation using
LAMMPS.
- ii.
Compute the HH autocorrelation function G R ,T (t ) to monitor the evolution of
r ij (t) and θ ij (t) for both intermolecular and intramolecular interactions.
- iii.
Apply equations (3) or (8) to obtain the respective T 2,R,T relaxation rates.
- iv.
T 2 relaxation time is finally obtained from equations (5) or (9)

All your files and folders are presented as a tree in the file explorer. You can switch from one to another by clicking a file in the tree.

## Input and Output

The input file:
tsim 25000 ! Total number of snapshots
dt 10 ! number of windows
twin 2500! size of window
dumpfile dump-water.lammpstrj ! the lammps trajectory file
outfile water-out.dat ! the output file

In the main.f90 file, modify the following depending on the system
MOL1 ! number of atoms per molecule
HNO ! number of hydrogen atoms per molecule
NATOMS ! Total number of atoms
NMOL ! Total number of molecules


The output file contains the following:

Line 1: G R (0), G T (0)
Column 1: time (t)
Column 2: G R (t)/G R (0)
Column 3: G T (t)/G T (0)
Column 4: G R (t)
Column 5 : G T (t)
To run the code : ./compile
: ./corr.x < input
