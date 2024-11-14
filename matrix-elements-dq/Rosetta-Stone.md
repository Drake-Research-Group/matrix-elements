# Rosetta Stone #

This document includes all of the changes made from Drake's original Fortran77 code to the modernized Fortran2018 present within this GitHub repo.
It will detail all of the variable name changes, common practices, and structure changes to the code.

## Variable Name Changes ##

- PA, PB --> powers_R1_and_R2_A, powers_R1_and_R2_B.
  - Contains the powers of R1 and R2 for wavefunction A and wavefunction B respectively.

- SA, SB --> powers_R12_A, powers_R12_B
  - Contains the powers of R12 for wavefunction A and wavefunciton B respectively.

- IZA, IZB --> nuclearCharge_A, nuclearCharge_B
  - Variables for the nuclear charge of wavefunctions A and B respectively.

- LRGLA, LRGLB --> angularMomentum_A, angularMomentum_B
  - Variables for the angular momentum of wavefunctions A and B.

- NSPNA, NSPNB --> totalSpin_A, totalSpin_B
  - Varaibles for the total spin of wavefunction A and B.

- NEIGA, NEIGB --> eigenstate_A, eigenstate_B
  - Represents the eigenstates for wavefunctions A and B.

- WA, WB --> waveFnEnergies_A, waveFnEnergies_B
  - An array containing the energy of the wavefunction as well as the screened hydrogenic energy for wavfunctions A and B respectively.

- DA, DB --> waveFnCoefficients_A, waveFnCoefficients_B
  - An array containing all of the coefficients for wavefunctions A and B

- NBLKA, NBLKB --> waveFnSectorStart_A, waveFnSectorStart_B
  - An array containing the term number where the next sector of the wavefunction starts.

- NBXA, NBXB --> listCounter_A, listCounter_B
  - Variable containing the total number of non-linear parameters in the wavefunction.

- NBLXA, NBLXB --> sectorList_A, sectorList_B
  - An array containing a list of the sectors of the wavefunctions for wavefunctions A and B. Ex. (5 sectors -> [1, 2, 3, 4, 5]) Seems be to changed slightly so that the array is actually [0, 1, 3, 4, 5] and skips 2, but I'm not sure why yet. (Looking into it currently)

- NB1, NB2, NB3, NB4 --> sectorStart_A, sectorStop_A, sectorStart_B, sectorStop_B
  - Variables containing the term number where the sector of the wavefunction starts and stops. Used in *breit* commonly in the do loops to loop over the sectors of the wavefunctions.

- ERMAC --> machineError
  - Variable that determines the precision error on the CPU you are working on. Updated the name for further clarity.

- NW --> numTerms
  - Variable containing the total number of terms in the basis set.

- NAMA, NAMB --> waveFnName_A, waveFnName_B
  - Variables containing the name of the wavefunction being processed.

- AMM, AMMB --> reducedMassRatio_A, reducedMassRatio_B
  - Variable containing the reduced mass ratio. Renamed for clarity.

- NEIG0 --> eigenstateZero
  - The first eigenstate. For large basis set sizes this is equivalent to choosing the principle quantum number.

- YA, YB --> nonLinParams_A, nonLinParams_B
  - Arrays containing the $\alpha$'s and $\beta$'s for the respective wavefunction.

- MP1, MP2 --> lowest_pow_r1, highest_pow_r1
  - variables containing the lowest and highest power of R1

- MQ1, MQ2 --> lowest_pow_r2, highest_pow_r2
  - variables conatining the lowest and highest power of R2

- MS1, MS2 --> lowest_pow_R12, highest_pow_R12
  - Variables containing the lowest and highest power of R12

- L1P, L2P --> l1_prime, l2_prime
  - Renamed the $l^\prime_1$ and $l^\prime_2$ variables for more clarity

- LTOT, LTOTP --> l_total, l_total_prime
  - Renamed for further clarity, these variables contain the total $L$ and total $L^\prime$ for the system.

- M1P, M2P --> m1_prime, m2_prime
  - Renamed the $m_1^\prime$ and $m_2^\prime$ varaibles for more clarity.

- MTOT, MTOTP --> m_total, m_total_prime
  - Renamed for further clarity, these variables contain the total $M$ and total $M^\prime$ for the system.

- ISUM --> sumOfPowers
  - Variable which contains the sum of powers of R1, R2, and R12. This is for use of the expression $i + j + k \le \Omega$

- FAC --> factorial
  - Array containing the factorial numbers, after calculated once, the values are stored here for faster calculation. renamed for clarity.

- R --> harmonicSeries
  - Array containing the numbers of the harmonic Series. After calculated once, the values are stored here for faster calulation. Renamed for clarity.

- BIN --> binomialCoefficient
  - Array containing the numbers for binomial coefficients. After calculated once, the values are stored here for faster calculation. Renamed for clarity.

- F21X --> hyperGeometricFn
  - Array containing the numbers for the hypergeometric function. After calculated once, the values are stored here for faster calculation. Renamed for clarity.

- LA, LB --> orbitalAngMomentum_A, orbitalAngMomentum_B
  - Array containing the orbital angular momenta for wavefunctions A and B respectively. It is a 2D array of size (2, :), the first row (1, :) is for electron 1, the second row (2, :) is for electron 2.

- LINCS --> increasedOrbitalAngMomentum
  - Variable pertaining to the increased Orbital Angular Momentum. renamed for clarity.

- HTERM, HTSCR --> matrix_element, screened_hydrogrenic
  - Variables containing the matrix element and screened hydrogenic results. Printed at the end of the program, these variables contain the programs calculation results.

- LTOTA, LTOTB --> totalAngMomentum_A, totalAngMomentum_B
  - Variables containing the total angular momentum for state A and state B.

- P, Q, S --> R1, R2, R12
  - Variables containing the R1, R2, and R12 values. Renamed for Clarity. Only changed in matrix-elements as the letters mean something else in the other library files.

## Code Structure Changes ##

- All instances of COMMON blocks have been removed. Global variables are instead initialized in a global-variables-dq.f90 module, which is imported into the program. This allows us to only write the variable definitions in one place instead of 
re-typing the variables in every function we want to use. It also makes it so that a size or type change in an array or variable only needs to be changed in one place, and will act everywhere at once, instead of having to change every single COMMON block.
- The process present in **SPIN.f** have been moved into seperate files. I am trying to avoid having files with thousands of lines in them as it makes the programs hard to edit, understand, and track down errors. Instead the different processes inside **SPIN.f**
have been moved to seperate files such as **recursive-relations-dq.f90**, **generate-integrals-dq.f90**, **Hypergeometric.f90**, etc.
- Similarly to **SPIN.f**, **CROSS.f** has been modified and split into separate files which are only concerned with one process at a time. In this case **CROSS.f** has been split into **wigner-coefficients-dq.f90** and **angular-coefficients0dq.f90**.
- All *GO TO* statements have been removed as they make the code progress non-linearly. Instead, if statements were modified to encapsulate code instead of jumping past blocks of code. This gives the code a more modern structure. 
- Include statements have been removed and replaced by thier respective modules.

## Common Practices ##

- All if statements are structured as *if* *end if* 
- All do loops are structed as *do* *end do*. No line label termination on do loops should be present. 
