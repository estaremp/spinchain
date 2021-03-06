This code solves the Schrodinger equation of a set of interacting qubits in a system defined under the XY Heisenberg Hamiltonian. The dynamics of the system is obtained by evolving the eigenstates (obtained from direct diagonalization of the Hamiltonian) and other parameters such as Fidelity or Entanglement of Formation are computed.

It also allows for the addition of static perturbations and time delays between injections. It computes several disorder realisations and outputs the average of the resultant parameters. 

Plots created through python (matplotlib) scripts are also generated upon request. 

See https://github.com/estaremp/spinchain/wiki for more info.

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
* Python 3.7. (https://www.python.org/download/releases/python-370/)
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
