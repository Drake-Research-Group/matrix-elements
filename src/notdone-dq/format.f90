!? This file is used to store all of the format statements used in various programs.
!? It is designed to help make the code appear clearer and more readable for the user.
!* Developed by Evan Petrimoulx June 21 2023

module format
    implicit none

    character(len = 25),  parameter :: fmt1 = '(" TERM",i2," =",2d20.12)'
    character(len = 22),  parameter :: fmt2 = '(1x,1p,d18.11,2d22.15)'
    character(len = 6),   parameter :: fmt3 = '(20i3)'
    character(len = 28),  parameter :: fmt4 = '(i4,1p,4d18.11/(4x,4d18.11))'
    character(len = 38),  parameter :: fmt5 = '(" GENINT",i1,6i3,2f10.6,4i2,3x,2i2,$)'
    character(len = 18),  parameter :: fmt6 = '(2i4,2x,3i4,d15.7)'
    character(len = 40),  parameter :: fmt7 = '(" NO CONTRIBUTION DUE TO SPIN ORTHOG.")'
    character(len = 45),  parameter :: fmt8 = '("RETURN CALLED DUE TO VANISHING 3-J SYMBOL")'
    character(len = 23),  parameter :: fmt9 = '(" TRAP 1", 9i4,2d12.5)'
    character(len = 10),  parameter :: fmt10 = '(a51)'
    character(len = 30),  parameter :: fmt11 = '(3i3,i5,i3,1x,a12,a52)'
    character(len = 8) ,  parameter :: fmt12 = '(16i5)'
    character(len = 20),  parameter :: fmt13 = '(3i5,i2,2f20.12)'
    character(len = 20),  parameter :: fmt14 = '(10(i3, 2i2))'
    character(len = 132), parameter :: fmt15 = '(1h ,3hz =,i3,6h   l =,i3,6h   s =,i3,6h   n =,i4,"   NEIGA =",i3/1x,a12,1x,a52)'
    character(len = 30),  parameter :: fmt16 = '(1x, 3d26.19)'
    character(len = 40),  parameter :: fmt17 = '(" no contribution from",2x,4i3,3x,4i3)'
    character(len = 56),  parameter :: fmt18 = '(28h error- negative j in threej )'
    character(len = 66),  parameter :: fmt19 = '(33h error- m larger than j in threej )'
    character(len = 40),  parameter :: fmt20 = '(" ERROR IN CROSS, lt = ",i4)'
    character(len = 40),  parameter :: fmt21 = '(" exit called, itot =",i6)'
    character(len = 40),  parameter :: fmt22 = '(" NO CONVERGENCE IN TSUM",2d15.6)'
    character(len = 40),  parameter :: fmt23 = '(12h  tabb error)'
    character(len = 30),  parameter :: fmt24 = '(" ", A10, " tab error")'
    character(len = 30),  parameter :: fmt25 = '(" TC3AUX ERROR",3i4)'
    character(len = 30),  parameter :: fmt26 = '(" TC3X ERROR",3i4)'
    character(len = 30),  parameter :: fmt27 = '(" TC33 ERROR",3i4)'
    character(len = 30),  parameter :: fmt28 = '(" TC32 ERROR",3i4)'
    character(len = 30),  parameter :: fmt29 = '(" TC3 ERROR",3i4)'
    character(len = 40),  parameter :: fmt30 = '(2i5," OUT OF RANGE IN F21")'
    character(len = 30),  parameter :: fmt31 = '("(",i2,",",i6,")")'
    character(len = 40),  parameter :: fmt32 = '(" J1+J2+J3 NOT EVEN", 3i5)'
    character(len = 85),  parameter :: fmt33 = '(43h ***** negative angular momentum in del2vw.)'
    character(len = 73),  parameter :: fmt34 = '(36h ***** fault in call of wbar01,isig=,i2)'
    character(len = 40),  parameter :: fmt35 = '(" error in cross, lmax1, lmax2 = ",2i4)'
    character(len = 14),  parameter :: fmt36 = '(9i5)'
    character(len = 40),  parameter :: fmt37 = '(" j1+j2+j3 not even", 3i5)'
    character(len = 60),  parameter :: fmt38 = '(1p,4d22.14)'
    character(len = 60),  parameter :: fmt39 = '(1x,1pd18.11,2d22.15,d10.3)'
    character(len = 40),  parameter :: fmt40 = '(d22.16)'
    character(len = 40),  parameter :: fmt41 = '(5i4,3d12.5)'
    character(len = 40),  parameter :: fmt42 = '(1x,1p,4d18.11)'
    
    contains
end module format