# CBSsiesta
This code takes the .HSX file of the code Siesta (https://launchpad.net/siesta) in order to calculate the Complex Band Stucture 
(CBS) of a material in one of the three directions perpendicular to the plane formed by two lattice vectors.

## Installation
This code is meant to work as a post-processing tool of Siesta. Therefore it requires a Siesta installation. I tested the code 
only using Siesta-4.0.1. However I do not see any reason why Siesta >= 3.2 distributions
should not work. Previous versions of Siesta do not produce the HSX file, so they are not supported.
The code uses the arch.make of Siesta for the installation. It is usually placed in Obj folder of Siesta. If this is the case,
it is sufficient to download the present code and place it in the folder Utils of the Siesta distribution.
Then make clean, make cbs should do the rest.
In case the arch.make is somewhere else, some modifications needs to be applied (see Makefile).

## Usage
The code requires the .HSX, .ORB_INDX and the .XV file obtained from the Siesta calculation. This means that the keyword
SaveHS must be set to true in the original Siesta calculation.
The input file for the cbs code must be called cbs.in and it should containd the following keywords:
* SystemLabel. It is the prefix of the .HSX (and .XV, .ORB_INDX) files.
* CBS.dir. Accepts 1,2,3. It defines the direction along which to calculate the CBS. 1 means that the CBS is calculated along
the direction perpendicular to the plane formed by the lattice vectors a_2 and a_3. 2 gives the CBS along the direction 
perpendicular to the the plane formed by the lattice vectors a_3 and a_1. 3 direction perpendicular to
the plane formed by the lattice vectors a_1 and a_2.
* CBS.Emin. Requires a float and units (for instance -35 eV). It is the minimum energy analysed by the code. 
* CBS.Emax. A flot and units as above. It is the maximum energy analysed by the code.
* CBS.bins. Number of energies in the range [CBS.Emin,CBS.Emax] analysed by the code.
* CBS.2DBZk1. Float. Defines, together with the next keyword, the k parallel at which the CBS analysis is performed. 
* CBS.2DBZk2. Float. 
