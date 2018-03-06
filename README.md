Abstract and general info

# Content 
* Schrodinger equation 
    * Integration method **(TBI)**
    * Diagonalization method 
* Generation of coupling scheme 
    * Uniform
    * Dimerised
        * SSH(a) and SSH(b) - (see https://arxiv.org/abs/1609.07516)
        * ABC - (see https://arxiv.org/abs/1612.05097)
    * PST - (see https://arxiv.org/abs/1603.01881)        
* Generation of XY Hamiltonian
    * Linear chain
    * Star networks
    * Lattices **(TBI)**
* Adding errors
    * Diagonal random disorder 
    * Off diagonal random disorder 
    * Injection time delays **(TBI)**
* Diagonalization of the Hamiltonian
    * Eigenvectors and Eigenenergies
    * Site occupation probabilities
* Dynamics
    * Fidelity
    * Site occupation probabilities
* Entanglement 
    * EOF **(TBI)**
    * Entropy **(TBI)**

***

# Compile and run

**Requirements:**
This is prepared be compiled by: 
* GNU-GCC 6.3.9 gfortran (https://gcc.gnu.org/wiki/GFortran)
* Python 2.7. (https://www.python.org/download/releases/2.7/)
* BLAS libraries (http://www.netlib.org/blas/) 
* LAPACK libraries (http://www.netlib.org/lapack/)
* MATPLOTLIB (https://matplotlib.org/)

NOTE: Still not tested in other versions/environments.

**Download and run:**
1. git clone https://github.com/estaremp/spinchain.git
2. set initial conditions in module PARAMETERS.f90
3. ./run.sh

***

# Output
* General output _spinchain.out_
* Data files _*.data_
* Plots _*.png_
