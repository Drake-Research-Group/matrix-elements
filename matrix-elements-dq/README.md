# Matrix Elements #

## Usage ##

Ensure that the wavefunction files you wish to run are in the main directory of this folder. This is a feature that will eventually be removed and pathed to the normal wave directory, but for now the files must be present in the matrix elements folder.

## Dependencies ##

## Matl88 ##

## Cross188 ##

## DQfun ##

## DQModule ##


## GenInt ##

$ia, ib, ic$ are the powers of $r_1.r_2,$ and $r_{12}$.

- t(ia, ib, ic) $\rightarrow$ The General Integration Function, calculates all positive powers for $r_1, r_2,$ and $r_{12}$
- tab(ia, ib, ic) $\rightarrow$ if $ia$ or $ib$ are negative powers.
- tc(ia, ib, ic) $\rightarrow$ if $ic$ is a negative power.
- tabx(ia, ib, ic) $\rightarrow$ Same functionality as tab but omits cancelled out terms. Used for integrals that appear in pairs.
- tlog(ia, ib, ic) $\rightarrow$ same functionality as t(ia, ib, ic) but used in the case that an additional $\log(r_{12})$ appears in the integrand.
- tc3(ia, ib, ic) $\rightarrow$ calculates $r_{12}^{-3}$ integrals.
- tc32(ia, ib, ic) $\rightarrow$ second method of calculating $r_{12}^{-3}$ integrals.
- tc33(ia, ib, ic) $\rightarrow$ version of tc3 for $ia = -3$ or $ib = -3$.
- tc3x(ia, ib, ic) same as tc3(ia, ib, ic) but omits cancelled out terms.
- tc3new(ia, ib, ic) $\rightarrow$ new version of tc3.
- tc3aux(ia, ib, ic) $\rightarrow$ special calculation of tc3
- CalcBinomialCoefficients $\rightarrow$ calculates the binomial coefficients. Useful for run speed.
- calc $\rightarrow$ Calculates frequently used common integrals.
- calcHypergeometricFn $\rightarrow$ $2F_1$ Calculates hypergeometric functions.
- csl(i, j, k) $\rightarrow$ Integrals containing $\cos (r_{12}) \log(r_{12}) functions.
- atlf $\rightarrow$ Integrals stored in the $at[ ]$ array. 
- atl $\rightarrow$ stores all of the integrals with a log term within them.
- calcRecursiveRelations $\rightarrow$ Calculates recursion relations.
- spl $\rightarrow$ Calls calcRecursiveRelations if a few extra conditions are met.
- splg $\rightarrow$ calls log function. Behaves the same as spl but for logs

## Spin88 ##
## Wavext ##
