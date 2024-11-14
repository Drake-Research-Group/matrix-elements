! /**
! * @author Evan Petrimoulx
! * @date July 10th 2024
! *
! * @brief This code calculates various matrix elements for 2 electron wavefunctions including:
! *   1. 1/R1^2
! *   2. 1/R1
! *   3. R1
! *   4. R1^2
! *   5. 1/R12^2
! *   6. 1/R12
! *   7. R12 
! *   8. R12^2 
! *   9. R1.R2 
! *   10. 1/(R1R12)
! *   11. 1/(R1R2)
! * 
! * @cite This code is based off of G.W.F. Drake'R12 original Fortran77 program "spin68".
! * @note The mathematical formulas can be found in the "Atomic Physics Handbook" - Chapter 11 - High Precision Calculations for Helium
! * 
! * @note Note the naming scheme of the program - variables with the suffix _A correspond to the first wavefunction, 
! * variables with the suffix _B correspond to the second wavefunction.
! */


program matrix_elements
  use dqmodule
  use wavExt
  use a1a_block
  use waveFileData
  use f1f_block
  use maxpow_block
  use radial_integrals
  use generate_integrals
  use binomial_coefficients

  implicit none

  ! Variable Initializations
  character(len = 4) :: dumpfile_type
  character(len = 12) :: waveFnName_A, waveFnName_B
  character(len = 17) :: waveFn_A, waveFn_B, matout, dumpfile_path
  character(len = 51), dimension(100) :: line

  logical :: screen, dumpfile_exists

  integer :: numTerms
  integer :: eigenstateZero
  integer :: i, j, k, l, m, n
  integer :: sectorStart, sectorStop ! Start and stop position corresponding to the number of terms in each spatial sector of the waveFn
  integer :: increasedOrbitalAngMomentum
  integer :: sumOfPowers_A, sumOfPowers_B, ih, lincs
  integer :: end_flag, end_of_file, iostat, ios
  integer, dimension(50) :: maxPartialWaveFnSize_A

  type(dq_real) :: da0, db0, XX
  type(dq_real) :: dqsign
  type(dq_real) :: reducedMassRatio_A, reducedMassRatio_B
  type(dq_real), dimension(40) :: matrix_element, screened_hydrogenic

  ! Determine OS in order to grab path for calling correct C++ program for Unicode Symbols
  character(len = 10) :: char_integer
  character(len = 20) :: OS_NAME
  character(len = 50) :: build_path
  character(len = 11), dimension(11) :: element_symbols

  ! Conversion values for reading in Real(16) and converting to type(dq_real)
  real(16) :: a16, b16, c16
  real(16), allocatable, dimension(:) :: da16, db16
  real(16), dimension(40) :: matrix_element_16, screened_hydrogenic_16
  
  ! Variable Initializations
  matrix_element = dqreal(0)
  screened_hydrogenic = dqreal(0)
  iostat = 0
  ih = 2
  da0 = dqreal(1)
  db0 = dqreal(1)
  screen = .false.
  reducedMassRatio_A = dqreal(0)

  element_symbols = ['<1/r_1^2>  ', '<1/r_1>    ', '<r_1>      ', '<r_1^2>    ', '<1/r_12^2> ', '<1/r_12>   ', '<r_12>     ', '<r_12^2>   ', '<r_1 . r_2>', '<1/r_1r_12>', '<1/r_1r_12>']

  ! Open input.dat for reading and output.dat for writing
  call routedOpen(5, FILE = 'input.dat', STATUS = 'UNKNOWN')
  call routedOpen(4, FILE = 'output.dat', STATUS = 'UNKNOWN')
  
  read(5, '(A51)') line(1)
  read(5,*) matout

  write(*, '(A51)') line(1)
  write(4, '(A51)') line(1)

  do while(.true.)
    read(5, *, iostat = end_of_file) waveFn_A, waveFn_B

    if(end_of_file /= 0) then 
      stop 
    end if 

    if(waveFn_A == 'EXIT') then 
      stop
    end if 

    open(1,file = waveFn_A, status = 'OLD')
    read(1,'(3i3,i5,i3,1x,a12)') nuclearCharge_A, angularMomentum_A, totalSpin_A, numTerms, eigenstate_A, waveFnName_A

    if(nuclearCharge_A == 1) then 
      eigenstate_A = 1
      ih = 1
    else
      ih = 2
    end if 

    eigenstateZero = eigenstate_A - 1
    numTerms_A = numTerms + eigenstateZero

    read(1, '(1x,3d30.23)')a16,b16,c16 
    read(1, '(I5)', advance = 'no') numNonLinParams_A
    
    ! Convert read in values to type(dq_real)
    waveFnEnergies_A(1) = dqreal(a16)
    waveFnEnergies_A(2) = dqreal(b16)
    reducedMassRatio_A = dqreal(c16)


    ! Now that the number of terms is known, we can allocate the proper sizes to the arrays pertaining to waveFn_A
    allocate(waveFnSectorStart_A(numNonLinParams_A + 1))
    allocate(orbitalAngMomentum_A(2, numTerms_A))
    allocate(powers_R1_and_R2_A(2, numTerms_A))
    allocate(waveFnCoefficients_A(numTerms_A))
    allocate(sectorList_A(numNonLinParams_A))
    allocate(nonLinParams_A(2, numTerms_A))
    allocate(powers_R12_A(numTerms_A))
    allocate(da16(numTerms_A))

    nonLinParams_A(1, 1) = dqreal(0)
    nonLinParams_A(2, 1) = dqreal(0)
    waveFnSectorStart_A(1) = 0

    read(1, '(15I5)')  (waveFnSectorStart_A(i + 1), i = 1, numNonLinParams_A), (maxPartialWaveFnSize_A(i + 1), i = 1, numNonLinParams_A - ih + 1)


    ! Print important wavefunction data to screen
    write(*, *)
    write(*, '(A)') waveFnName_A
    write(*, '("Z = ", i3.3, 3x, "L = ", i3.3, 3x, "R12 = ", i3.3, 3x, "N = ", i4.4)') nuclearCharge_A, angularMomentum_A, totalSpin_A, numTerms
    write(*, *)
    write(*, '(A, 5x)', advance = "no") "Screened Energy:"
    write(*, '(1p2e40.32)') waveFnEnergies_A(1)
    write(*, '(A)', advance = "no") "Wavefunction Energy: " 
    write(*, '(1p2e40.32)') waveFnEnergies_A(2) 
    write(*, '(A)', advance = "no") "Reduced Mass Ratio:  "
    write(*, '(1p2e40.32)') reducedMassRatio_A

    ! Output file records
    write(4, *)
    write(4, '(A)') waveFnName_A
    write(4, '("Z = ", i3.3, 3x, "L = ", i3.3, 3x, "R12 = ", i3.3, 3x, "N = ", i4.4)') nuclearCharge_A, angularMomentum_A, totalSpin_A, numTerms
    write(4, *)
    write(4, '(A, 5x)', advance = "no") "Screened Energy:"
    write(4, '(1p2e40.32)') waveFnEnergies_A(1)
    write(4, '(A)', advance = "no") "Wavefunction Energy: " 
    write(4, '(1p2e40.32)') waveFnEnergies_A(2) 
    write(4, '(A)', advance = "no") "Reduced Mass Ratio:  "
    write(4, '(1p2e40.32)') reducedMassRatio_A
    
    if(reducedMassRatio_A /= dqreal(0)) then
      dumpfile_type = '.DML' ! Finite Mass Case
    else 
      dumpfile_type = '.DMP' ! Infinite Mass Case
    end if 

    ! Check if DumpFile exists:
    dumpfile_path = trim('dump/' // waveFnName_A(:len(waveFnName_A)-4) // '.DMP')
    inquire(file = dumpfile_path, exist = dumpfile_exists)
    

    ! If DumpFile doesn't exist, use the coefficients stored in the wavefunction instead
    if(.not.dumpfile_exists) then 
      write(*, '(A)') "Dumpfile not found. Using coefficients from wavefunction file instead."
      dumpfile_type = ''
    end if

    listCounter_A = 1

    do i=1,numNonLinParams_A
      waveFnSectorStart_A(i+1) = waveFnSectorStart_A(i + 1) + eigenstateZero

      read(1, *) sectorStart, sectorStop, orbitalAngMomentum_A(1, sectorStart + eigenstateZero), orbitalAngMomentum_A(2, sectorStart + eigenstateZero), a16, b16
      nonLinParams_A(1, sectorStart + eigenstateZero) = dqreal(a16)
      nonLinParams_A(2, sectorStart + eigenstateZero) = dqreal(b16)
      
      write(*, *)
      write(*, '(A, I3)') "Basis Set Sector:     ", i
      write(*, '(A)', advance = "no") "Corresponding Alpha: "
      write(*, '(1p2e40.32)') nonLinParams_A(1, sectorStart + eigenstateZero)
      write(*, '(A)', advance = "no") "Corresponding Beta:  "
      write(*, '(1p2e40.32)') nonLinParams_A(2, sectorStart + eigenstateZero)
      write(4, *)
      write(4, '(A, I3)') "Basis Set Sector:     ", i
      write(4, '(A)', advance = "no") "Corresponding Alpha: "
      write(4, '(1p2e40.32)') nonLinParams_A(1, sectorStart + eigenstateZero)
      write(4, '(A)', advance = "no") "Corresponding Beta:  "
      write(4, '(1p2e40.32)') nonLinParams_A(2, sectorStart + eigenstateZero)

      if(sectorStart == 1.and.sectorStop == 1) then 
        screen = .true.
      end if 

      sectorStart = sectorStart + eigenstateZero
      sectorStop = sectorStop + eigenstateZero

      do j=sectorStart,sectorStop
        orbitalAngMomentum_A(:,j) = orbitalAngMomentum_A(:, sectorStart)
        nonLinParams_A(:,j) = nonLinParams_A(:, sectorStart)
      end do

      if(i == 1) then 
        cycle
      end if 

      if(nonLinParams_A(1, sectorStart) == nonLinParams_A(1, sectorStart - 1) .and. nonLinParams_A(2, sectorStart) == nonLinParams_A(2, sectorStart - 1)) then 
        cycle
      end if 

      listCounter_A = listCounter_A + 1
      sectorList_A(listCounter_A) = i - 1
    end do


    ! Original Code was wonky here and hard to follow. I believe this implementation is correct but if there is ever a problem with the sectorList array in the future, check here.
    sectorList_A(listCounter_A + 1) = numNonLinParams_A

    read(1,'(10(i3,2i2))') (powers_R1_and_R2_A(1,k), powers_R1_and_R2_A(2,k), powers_R12_A(k),k=eigenstate_A, numTerms_A)

    if(dumpfile_type == ".DMP") then 
      open(unit = 1000, file = 'dump/' // wavefnName_A(:len(trim(waveFnName_A))-4) // '.DMP', form='unformatted')
      read(1000) (waveFnCoefficients_A(i), i = eigenstate_A, numTerms_A)
      close(1000)
      print*, waveFnCoefficients_A
      stop
    else if(dumpfile_type == ".DML") then
      open(unit = 1000, file = 'dump/' // wavefnName_A(:len(trim(waveFnName_A))-4) // '.DML',  form='unformatted')
      read(1000) (waveFnCoefficients_A(i), i = eigenstate_A, numTerms_A)
      close(1000)
    else 
      read(1,'(1x,3d30.23)') (da16(i), i= eigenstate_A, numTerms_A)
      do i = eigenstate_A, numTerms_A
        waveFnCoefficients_A(i) = da16(i)
      end do 
    end if 
    
    call CalcBinomialCoefficients

    ! Calculation of screened hydrogenic wavefunctions
    if(screen) then
      increasedOrbitalAngMomentum = orbitalAngMomentum_A(1, eigenstate_A) + orbitalAngMomentum_A(2, eigenstate_A) - angularMomentum_A
      nonLinParams_A(1, 1) = nonLinParams_A(1, eigenstate_A)
      lincs = 0

      if(abs(nonLinParams_A(1,1) - dqreal(1) / dqreal(increasedOrbitalAngMomentum + 1)) > (dqreal(1) / dqreal(1000000))) then
        lincs = 1
      end if 

      nonLinParams_A(2,1) = dqreal(nuclearCharge_A-1)/dqreal(((angularMomentum_A + eigenstate_A) - lincs)*nuclearCharge_A)
      da0 = waveFnCoefficients_A(eigenstate_A)

      XX = dqreal(1)
      do m = 1, eigenstate_A
        orbitalAngMomentum_A(:, m) = orbitalAngMomentum_A(:, eigenstate_A)
        nonLinParams_A(:, m) = nonLinParams_A(:, 1)
        powers_R1_and_R2_A(1, m) = increasedOrbitalAngMomentum
        powers_R1_and_R2_A(2, m) = angularMomentum_A + m - 1
        powers_R12_A(m) = 0
        waveFnCoefficients_A(m) = waveFnCoefficients_A(eigenstate_A) * XX * dqreal(2)**(increasedOrbitalAngMomentum + angularMomentum_A + 2) * (dqreal(nuclearCharge_A) - dqreal(1))**(angularMomentum_A + 1) * sqrt((dqreal(nuclearCharge_A - 1)) &
                & * factorial( angularMomentum_A + eigenstate_A + angularMomentum_A + 1) / (dqreal(nuclearCharge_A) * factorial(angularMomentum_A + eigenstate_A - angularMomentum_A) &
                & * factorial(2 * increasedOrbitalAngMomentum + 2))) / (factorial(2 * angularMomentum_A + 2) * dqreal( angularMomentum_A + eigenstate_A) * (dqreal(increasedOrbitalAngMomentum) + dqreal(1))**(increasedOrbitalAngMomentum + 2) * (dqreal( angularMomentum_A + eigenstate_A) * dqreal(nuclearCharge_A))**(angularMomentum_A + 1))
        XX = -XX*dqreal(2*(nuclearCharge_A-1)*(eigenstate_A-M))/dqreal((2*angularMomentum_A+M+1)*M*(angularMomentum_A + eigenstate_A)*nuclearCharge_A)
      end do
    end if 

    ia = powers_R1_and_R2_A(1,1) + powers_R1_and_R2_A(2,1) + powers_R12_A(1)
    dqsign = dqreal(1)

    if(waveFnCoefficients_A(1) < dqreal(0)) then 
      dqsign = dqreal(-1)
    end if 

    da0 = da0 * dqsign

    do k=1,numTerms_A
      waveFnCoefficients_A(k) = waveFnCoefficients_A(k) * dqsign
      sumOfPowers_A = powers_R1_and_R2_A(1,k) + powers_R1_and_R2_A(2,k) + powers_R12_A(k)

      if(sumOfPowers_A > ia) then 
        ia = sumOfPowers_A
      end if 
    end do

    !   LEFT STATE READING ROUTINE - IF LEFT STATE IS THE SAME AS THE RIGHT STATE,
    !   USE THE FILE NAME 'SAME'.
    nuclearCharge_B = 0

    if(waveFn_B  /= 'SAME') then
      open(1, file = waveFn_B(:), status = 'OLD')

      read(1,'(3i3,i5,i3,1x,a12)') nuclearCharge_B, angularMomentum_B, totalSpin_B, numTerms, eigenstate_B, waveFnName_B

      if(nuclearCharge_B == 1) then 
        ih = 1
      else 
        ih = 2
      end if 

      eigenstateZero = eigenstate_B - 1
      numTerms_B = numTerms + eigenstateZero

      call dqread(1, waveFnEnergies_B(1), waveFnEnergies_B(2), reducedMassRatio_B)

      if(reducedMassRatio_A /= reducedMassRatio_B) then
        write(*,'(2D20.10,A)') reducedMassRatio_A,reducedMassRatio_B,' MASSES NOT EQUAL.'
        write(4,'(2D20.10,A)') reducedMassRatio_A,reducedMassRatio_B,' MASSES NOT EQUAL.'
        stop
      endif

      read(1,'(I5)', advance = 'no') numNonLinParams_B

      ! NumTerms_B is now known, so now we can allocate all array sizes for waveFn B
      allocate(waveFnSectorStart_B(numNonLinParams_B + 1))
      allocate(orbitalAngMomentum_B(2, numTerms_B))
      allocate(powers_R1_and_R2_B(2, numTerms_B))
      allocate(waveFnCoefficients_B(numTerms_B))
      allocate(sectorList_B(numNonLinParams_B))
      allocate(nonLinParams_B(2, numTerms_B))
      allocate(Powers_R12_B(numTerms_B))
      allocate(db16(numTerms_B))

      read(1, '(15I5)') (waveFnSectorStart_B(i+1),i=1,numNonLinParams_B),(maxPartialWaveFnSize_A(i+1),i=1,numNonLinParams_B-ih+1)
      
      screen = .false.

      sectorList_B(1) = 0
      listCounter_B = 1

      do i=1,numNonLinParams_B
        waveFnSectorStart_B(i+1) = waveFnSectorStart_B(i+1) + eigenstateZero
        
        read(1,2) sectorStart, sectorStop, orbitalAngMomentum_B(1,sectorStart+eigenstateZero), orbitalAngMomentum_B(2,sectorStart+eigenstateZero)
        call dqread(1, nonLinParams_B(1,sectorStart+eigenstateZero), nonLinParams_B(2,sectorStart+eigenstateZero))
        
        if(sectorStart == 1.and.sectorStop == 1) then
          screen = .true.
        end if

        sectorStart = sectorStart + eigenstateZero
        sectorStop = sectorStop + eigenstateZero

        do k=1,2
          do j=sectorStart,sectorStop
            orbitalAngMomentum_B(k,j) = orbitalAngMomentum_B(k,sectorStart)
            nonLinParams_B(k,j) = nonLinParams_B(k,sectorStart)
          end do
        end do 

        if(i == 1) then
          cycle
        end if 

        if(nonLinParams_B(1,sectorStart) == nonLinParams_B(1,sectorStart-1).and.nonLinParams_B(2,sectorStart) == nonLinParams_B(2,sectorStart-1)) then
          cycle 
        end if 

        listCounter_B = listCounter_B + 1
        sectorList_B(listCounter_B) = i - 1
      end do

      sectorList_B(listCounter_B + 1) = numNonLinParams_B

      read(1,'(10(i3,2i2))')(powers_R1_and_R2_B(1,k),powers_R1_and_R2_B(2,k),powers_R12_B(k),k=eigenstate_B,numTerms_B)

      do i = eigenstate_B, numTerms_B
        call dqread(1, waveFnCoefficients_B(i))
      end do 

      write(4,'("Z = ",i3.3, 3x,"L = ",i3.3, 3x,"R12 = ",i3.3, 3x,"n = ",i4.4)') nuclearCharge_B, angularMomentum_B, totalSpin_B, numTerms
      write(*,'("Z = ",i3.3, 3x,"L = ",i3.3, 3x,"R12 = ",i3.3, 3x,"n = ",i4.4)') nuclearCharge_B, angularMomentum_B, totalSpin_B, numTerms

      !   CALCULATION OF SCREENED HYDROGENIC WAVEFUNCTIONS.
      if((nuclearCharge_B /= 0) .and. (screen)) then
        increasedOrbitalAngMomentum = orbitalAngMomentum_B(1, eigenstate_B) + orbitalAngMomentum_B(2, eigenstate_B) - angularMomentum_B
        nonLinParams_B(1,1) = nonLinParams_B(1,eigenstate_B)
        lincs = 0

        if(abs(nonLinParams_B(1, 1) - dqreal(1) / dqreal(increasedOrbitalAngMomentum + 1)) > (dqreal(1) / dqreal(1000000))) then
          lincs = 1
        end if 

        nonLinParams_B(2,1) = dqreal((nuclearCharge_B - 1) / (( angularMomentum_A + eigenstate_A - lincs) * nuclearCharge_B))
        db0 = waveFnCoefficients_B(eigenstate_B)

        XX = dqreal(1)
        do m=1,eigenstate_B
          do k=1,2
            orbitalAngMomentum_B(k,m) = orbitalAngMomentum_B(k,eigenstate_B)
            nonLinParams_B(k,m) = nonLinParams_B(k,1)
          end do 

          powers_R1_and_R2_B(1,m) = increasedOrbitalAngMomentum
          powers_R1_and_R2_B(2,m) = angularMomentum_B + m - 1
          powers_R12_B(m) = 0
          waveFnCoefficients_B(m) = waveFnCoefficients_B(eigenstate_B) * dqreal(1) * dqreal(2) ** (increasedOrbitalAngMomentum + angularMomentum_B + 2) * dqreal(nuclearCharge_B - 1) ** (angularMomentum_B + 1) * sqrt(dqreal(nuclearCharge_B - 1) * factorial( angularMomentum_A + eigenstate_A + angularMomentum_B + 1) & 
                & / (dqreal(nuclearCharge_B) *factorial( angularMomentum_A + eigenstate_A - angularMomentum_B) * factorial(2 * increasedOrbitalAngMomentum + 2))) / (factorial(2 * angularMomentum_B + 2) * dqreal( angularMomentum_A + eigenstate_A) * dqreal(increasedOrbitalAngMomentum + 1) ** (increasedOrbitalAngMomentum + 2) * (dqreal( angularMomentum_A + eigenstate_A) * dqreal(nuclearCharge_A)) ** (angularMomentum_B + 1))
          XX = -XX*dqreal(2*(nuclearCharge_A-1)*(eigenstate_A-M))/dqreal((2*angularMomentum_A+M+1)*M*(angularMomentum_A + eigenstate_A)*nuclearCharge_A)
        end do 
      end if 

    else
      angularMomentum_B = angularMomentum_A
      totalSpin_B = totalSpin_A
      eigenstate_B = eigenstate_A
      waveFnName_B = 'SAME'
      numTerms_B = numTerms_A
      numNonLinParams_B = numNonLinParams_A
      listCounter_B = listCounter_A
      

      ! Allocate Array sizes using numTerms_B
      allocate(waveFnSectorStart_B(numNonLinParams_B + 1))
      allocate(orbitalAngMomentum_B(2, numTerms_B))
      allocate(powers_R1_and_R2_B(2, numTerms_B))
      allocate(waveFnCoefficients_B(numTerms_B))
      allocate(sectorList_B(numNonLinParams_B))
      allocate(nonLinParams_B(2, numTerms_B))
      allocate(Powers_R12_B(numTerms_B))
      allocate(db16(numTerms_B))

      do i=1,numNonLinParams_B
        waveFnSectorStart_B(i+1) = waveFnSectorStart_A(i+1)
        sectorList_B(i) = sectorList_A(i)
      end do

      do k=1,numTerms_B
        do i=1,2
          orbitalAngMomentum_B(i,k) = orbitalAngMomentum_A(i,k)
          nonLinParams_B(i,k) = nonLinParams_A(i,k)
          powers_R1_and_R2_B(i,k) = powers_R1_and_R2_A(i,k)
        end do 
        powers_R12_B(k) = powers_R12_A(k)
        waveFnCoefficients_B(k) = waveFnCoefficients_A(k)
      end do

      db0 = da0
    end if 

    ib = powers_R1_and_R2_B(1,1) + powers_R1_and_R2_B(2,1) + powers_R12_B(1)
    dqsign = dqreal(1)

    if(waveFnCoefficients_B(1) < dqreal(0)) then 
      dqsign = dqreal(-1)
    end if 

    db0 = db0 * dqsign

    do k = 1, numTerms_B
      waveFnCoefficients_B(k) = waveFnCoefficients_B(k)*dqsign
      sumOfPowers_B = powers_R1_and_R2_B(1,k) + powers_R1_and_R2_B(2,k) + powers_R12_B(k)

      if(sumOfPowers_B > ib) then 
        ib = sumOfPowers_B
      end if 
    end do
    
    do n=1,numNonLinParams_A
      sectorStart = waveFnSectorStart_A(n) + 1
      sectorStop = waveFnSectorStart_A(n+1)
      iia(n) = 0
      jja(n) = 0
      kka(n) = 0
      ijka(n) = 0

      do i = sectorStart, sectorStop
        if(powers_R1_and_R2_A(1, i) > iia(n)) then
          iia(n) = powers_R1_and_R2_A(1,i)
        end if 

        if(powers_R1_and_R2_A(2,i) > jja(n)) then
          jja(n) = powers_R1_and_R2_A(2,i)
        end if 

        if(powers_R12_A(i) > kka(n)) then 
          kka(n) = powers_R12_A(i)
        end if 

        j = powers_R1_and_R2_A(1,i) + powers_R1_and_R2_A(2,i) + powers_R12_A(i)

        if(j > ijka(n)) then 
          ijka(n) = j
        end if 
      end do
    end do


    waveFnSectorStart_B(1) = 0

    do n=1,numNonLinParams_B
      sectorStart = waveFnSectorStart_B(n) + 1
      sectorStop = waveFnSectorStart_B(n+1)
      iib(n) = 0
      jjb(n) = 0
      kkb(n) = 0
      ijkb(n) = 0

      do i=sectorStart, sectorStop
        if(powers_R1_and_R2_B(1,i) > iib(n)) then 
          iib(n) = powers_R1_and_R2_B(1,i)
        end if

        if(powers_R1_and_R2_B(2,i) > jjb(n)) then
          jjb(n) = powers_R1_and_R2_B(2,i)
        end if 

        if(powers_R12_B(i) > kkb(n)) then 
          kkb(n) = powers_R12_B(i)
        end if 

        j = powers_R1_and_R2_B(1,i) + powers_R1_and_R2_B(2,i) + powers_R12_B(i)

        if(j > ijkb(n)) then 
          ijkb(n) = j
        end if 
      end do
    end do

    close(1,status='KEEP')

    call breit(matrix_element, screened_hydrogenic)
    
    do i=1,40
      matrix_element(i) = matrix_element(i)
      screened_hydrogenic(i) = screened_hydrogenic(i)
      screened_hydrogenic(i) = screened_hydrogenic(i) / (dqreal(2) * da0 * db0)
    end do


    ! Print results to screen
    write(*, *) ""
    write(*, '(A)') "Matrix Elements:"
    write(*, '(A)') "-------------------"
    write(4, *) ""
    write(4, '(A)') "Matrix Elements:"
    write(4, '(A)') "-------------------"

    do i = 1, 31, 3
      matrix_element_16(i) = matrix_element(i)
      screened_hydrogenic_16(i) = screened_hydrogenic((i + 2) / 3)

      write(char_integer, '(I0)') (i + 2)/3
      write(*, '(A11)') element_symbols((i + 2) / 3)
      write(*, '(E23.15E2)') screened_hydrogenic_16(i)
      write(*, '(E23.15E2)') matrix_element_16(i)
      write(*, *) " "
      write(4, '(A11)') element_symbols((i + 2) / 3)
      write(4, '(E23.15E2)') screened_hydrogenic_16(i)
      write(4, '(E23.15E2)') matrix_element_16(i)
      write(4, *) " "
    end do
    
    ! Calculation is finished, deallocate arrays to prevent memory leaks
    deallocate(da16)
    deallocate(db16)
    deallocate(powers_R12_A)
    deallocate(powers_R12_B)
    deallocate(sectorList_A)
    deallocate(sectorList_B)
    deallocate(nonLinParams_A)
    deallocate(nonLinParams_B)
    deallocate(powers_R1_and_R2_A)
    deallocate(powers_R1_and_R2_B)
    deallocate(waveFnSectorStart_A)
    deallocate(waveFnSectorStart_B)
    deallocate(orbitalAngMomentum_A)
    deallocate(orbitalAngMomentum_B)
    deallocate(waveFnCoefficients_A)
    deallocate(waveFnCoefficients_B)

    open(1, file = 'data/other/' // matout, status = 'OLD', access = 'append')
    backspace 1
    write(1,'(A14,1X,A14)') waveFn_A, waveFn_B

    ! Format in real(16)
    do i = 1, 31, 3
      write(1, '(A, I2, A, 2D21.9)') 'TERM', (i+2)/3, ' =', matrix_element_16(i), screened_hydrogenic_16(i)
    end do 

    write(*,'(A51)') line(1)
    write(*,'(A51)')
    write(4,'(A51)') line(1)
    write(4,'(A51)')
    close(1,status='KEEP')

    ! RESET INPUT FILE FOR NEXT CALCULATION.
    cycle
    rewind 5

    do i=1,100
      read(5,'(A)', iostat = end_flag) line(i)
      if(end_flag /= 0) then 
        rewind 5
        exit 
      end if 
    end do

    write(5,'(A)') line(1)
    write(5,'(A16)') line(2)

    do l = 4, i - 1
      write(5,'(A32)') line(l)
    end do

    write(5,'(A32)') line(3)

    rewind 5

    read(5,'(/)')

  end do 

  2  format(3i5,i2,2f20.12)

end program matrix_elements


subroutine breit(matrix_element,screened_hydrogenic)
  use dqmodule
  use a1a_block
  use waveFileData
  use c1c_block
  use d1d_block
  use f1f_block
  use g1g_block
  use maxpow_block
  use generate_integrals
  use recursive_relations
  use angular_coefficients
  implicit none

  logical :: lgo, exchange
  
  character(len = 80) :: buffer
  
  type(dq_real), dimension(40) :: matrix_element, screened_hydrogenic
  type(dq_real) :: psi2, dqsign, sum, test

  integer :: R1, R2, R12, ii, jj
  integer :: idx, jdx_, i, j, m, n ! Loop iterators
  integer :: id, jd, idep, jdep
  integer :: sectorStart_A, sectorStop_A
  integer :: totalAngMomentum_A, totalAngMomentum_B
  integer :: lmax, temp, sectorStart_B, sectorStop_B
  integer :: loop_counter,l1x,l2x,l3x,l4x

  integer, save :: sum_d, highest_pow_r1d, highest_pow_r2d, highest_pow_r12d

 data l1x,l2x,l3x,l4x/4*100/
  loop_counter = 0
  ! FOR D STATES
  sum_d = 2
  highest_pow_r1d = 3
  highest_pow_r2d = 3 
  highest_pow_r12d = 1 
 
  ! FOR R1 STATES
  ! sum_d = 2
  ! highest_pow_r1d = 3
  ! highest_pow_r2d = 3
  ! highest_pow_r12d = 1

  machineError = dqreal(1)
  test = dqreal(1)

  ! Calculate Machine Precision
  open(2, file = 'data/machineError.dat')
  do while(.true.)
    machineError = machineError / dqreal(2)
    test = test + machineError

    rewind(2)
    write(2, *) test
    rewind(2)
    read(2,'(A)') buffer

    if (buffer(4:4) == '2') then 
      exit
    end if 
  end do 
  close(2, status = 'DELETE')

  machineError = dqreal(2) * machineError

  write(*,'(A)', advance = "no") "Machine Error:       "
  write(*, *) dqreal(machineError)
  write(*, *)
  write(4,'(A)', advance = "no") "Machine Error:       "
  write(4, *) dqreal(machineError)
  write(4, *)

  kross = 9
  klog = -2

  ! Set matrix_element and screened_hydrogenic to 0
  matrix_element(:) = dqreal(0)
  screened_hydrogenic(1:11) = dqreal(0)

  ! Summation over the right hand state direct and exchange terms (idx == integer, direct-exchange)
  do idx = 1, 2
    do i = 1, numNonLinParams_A
      sectorStart_A = waveFnSectorStart_A(i) + 1
      sectorStop_A = waveFnSectorStart_A(i + 1)

      do j = 1, numNonLinParams_B
        sectorStart_B = waveFnSectorStart_B(j) + 1
        sectorStop_B = waveFnSectorStart_B(j + 1)

        lgo = .false.

        if(i == 1 .and. j == 1) then
          lgo = .true.
        end if 

        if(nuclearCharge_B == 0 .and. j < i) then 
          cycle 
        end if 

        !   SUMMATION OVER LEFT HAND STATE DIRECT AND EXCHANGE TERMS.
        !   THE ORDERING IS D-D, E-E, D-E AND E-D.

        do jdx_ = 1, 2
          ide =  jdx_
          jde =  jdx_

          if(idx == 2) then 
            ide = 3 - jde
          end if 

          id = ide - 1
          jd = jde - 1
          idep = 3 - ide
          jdep = 3 - jde
          dqsign = dqreal(1)
          
          if(totalSpin_A == 1) then 
            dqsign = dqreal(-1)**(ide-1)
          end if 

          if(totalSpin_B == 1) then 
            dqsign = dqsign * dqreal(-1)**(jde - 1)
          end if 

          l1 = orbitalAngMomentum_A(ide,sectorStart_A)
          l2 = orbitalAngMomentum_A(idep,sectorStart_A)
          l1_prime = orbitalAngMomentum_B(jde,sectorStart_B)
          l2_prime = orbitalAngMomentum_B(jdep,sectorStart_B)

          if(l1x.ne.l1.or.l2x.ne.l2.or.l3x.ne.l1_prime.or.l4x.ne.l2_prime) then
            call cross(l1_prime, l2_prime, l1, l2, angularMomentum_B, angularMomentum_B, angularMomentum_A, angularMomentum_A, ide, jde)
            l1x = l1
            l2x = l2
            l3x = l1_prime
            l4x = l2_prime
          endif

          do jj = 1, 8
            mlt(jj) = lt(jj)
            do ii = lt(jj), 1, -2
              if(cplt(ii, jj) /= dqreal(0.q0)) then
                mlt(jj) = ii
              end if
            end do
          end do

          if(l1_prime == l1) then
            !   for cos(theta1)cos(theta2)-direct terms sp-sp or ps-ps
            mlt(9) = 2
            lt(9) = 2
            cplt(1,9) = dqreal(0)
            cplt(2,9) = dqreal(1)/dqreal(10)
            cplt(3,9) = dqreal(0)

            if(l1 == 0) then
              !   for cos^2(theta1)-direct terms sp-sp
              mlt(10) = 1
              lt(10) = 3
              cplt(1,10) = dqreal(1) / dqreal(6)
              cplt(2,10) = dqreal(0)
              cplt(3,10) = dqreal(-1) / dqreal(15)
              !   for cos^2(theta2)-direct terms sp-sp
              mlt(11) = 1
              lt(11) = 1
              cplt(1,11) = dqreal(1) / dqreal(10)
              cplt(2,11) = dqreal(0)
              cplt(3,11) = dqreal(0)
            else
              !   for cos^2(theta1)-direct terms ps-ps
              mlt(10) = 1
              lt(10) = 1
              cplt(1,10) = dqreal(1) / dqreal(10)
              cplt(2,10) = dqreal(0)
              cplt(3,10) = dqreal(0)
              !   for cos^2(theta2)-direct terms ps-ps
              mlt(11) = 1
              lt(11) = 3
              cplt(1,11) = dqreal(1) / dqreal(6)
              cplt(2,11) = dqreal(0)
              cplt(3,11) = dqreal(-1) / dqreal(15)
            end if
          else
            !    for cos(theta1)cos(theta2)-exchange terms sp-ps or ps-sp
            mlt(9) = 1
            lt(9) = 3
            cplt(1,9) = dqreal(1) / dqreal(6)
            cplt(2,9) = dqreal(0)
            cplt(3,9) = dqreal(-1) / dqreal(15)
            !    for cos^2(theta1)  - exchange terms sp - ps or ps - sp
            mlt(10) = 2
            lt(10) = 2
            cplt(1,10) = dqreal(0)
            cplt(2,10) = dqreal(1) / dqreal(10)
            cplt(3,10) = dqreal(0)
            !    for cos^2(theta2)  - exchange terms sp - ps or ps - sp
            mlt(11) = 2
            lt(11) = 2
            cplt(1,11) = dqreal(0)
            cplt(2,11) = dqreal(1) / dqreal(10)
            cplt(3,11) = dqreal(0)
          end if

          if(jdx_ /= 2) then 
            if(i <= i.and.j <= j) then
              y11 = nonLinParams_A(ide, sectorStart_A) + nonLinParams_B(jde, sectorStart_B)
              y22 = nonLinParams_A(idep, sectorStart_A) + nonLinParams_B(jdep, sectorStart_B)
              exchange = .false.


              if(y22 > y11) then
                y22 = nonLinParams_A(ide,sectorStart_A) + nonLinParams_B(jde,sectorStart_B)
                y11 = nonLinParams_A(idep,sectorStart_A) + nonLinParams_B(jdep,sectorStart_B)
                exchange = .true.
              end if

              lmax = lt(2)
              lowest_pow_r1 = 3
              lowest_pow_r12 = 3
              msl = lt(1)
              imin = -2 + angularMomentum_A + angularMomentum_B
              totalAngMomentum_A = angularMomentum_A + eigenstate_A - 1
              totalAngMomentum_B = angularMomentum_B + eigenstate_B - 1

              if(nuclearCharge_A <= 1.or.(sectorStart_A /= 1.and.sectorStart_B /= 1)) then 
                sumOfPowers = ijkb(j) + ijka(i) + sum_d
                highest_pow_r1 = iib(j) + (1-id)*iia(i) + id*jja(i) + lowest_pow_r1 + highest_pow_r1d+ 1
                highest_pow_r2 = jjb(j) + (1-id)*jja(i) + id*iia(i)+ lowest_pow_r1+highest_pow_r2d + id
                highest_pow_r12 = lowest_pow_r12 + kkb(j) + kka(i) + 2*lt(1) + highest_pow_r12d

              else if(sectorStart_B /= 1.or.sectorStart_A /= 1) then
                if(sectorStart_B /= 1) then 
                  !   LIMITS FOR HYDROGENIC '1,N' TERMS.
                  sumOfPowers = totalAngMomentum_A + ijkb(j) + sum_d
                  highest_pow_r1 = iib(j) + id*totalAngMomentum_A + lowest_pow_r1 + highest_pow_r1d + 1 - id
                  highest_pow_r2 = jjb(j) + (1-id)*totalAngMomentum_A + lowest_pow_r1 + highest_pow_r2d
                  highest_pow_r12 = lowest_pow_r12 + kkb(j) + 2*lt(1) + highest_pow_r12d + id

                else
                  !   LIMITS FOR HYDROGENIC 'N,1' TERMS.
                  sumOfPowers = totalAngMomentum_B + ijka(i) + sum_d
                  highest_pow_r1 = (1-id)*iia(i) + id*(jja(i)+totalAngMomentum_B) + lowest_pow_r1 + highest_pow_r1d + 1- id
                  highest_pow_r2 = (1-id)*(jja(i)+totalAngMomentum_B) + id*iia(i) + lowest_pow_r1 + highest_pow_r2d
                  highest_pow_r12 = lowest_pow_r12 + kka(i) + 2*lt(1) + highest_pow_r12d + id
                end if
              else
                sumOfPowers = totalAngMomentum_B + totalAngMomentum_A + sum_d
                highest_pow_r1 = (1-id)*totalAngMomentum_A + id*(totalAngMomentum_B+totalAngMomentum_A) + lowest_pow_r1 + highest_pow_r1d + 1 - id
                highest_pow_r2 = (1-id)*(totalAngMomentum_B+totalAngMomentum_A) + id*totalAngMomentum_A + lowest_pow_r1 + highest_pow_r2d
                highest_pow_r12 = lowest_pow_r12 + 2*lt(1) + highest_pow_r12d + id
              end if 

              if(id == 0.and.(j > j.or.i > i)) then 
                highest_pow_r12 = highest_pow_r12 +angularMomentum_A
              end if 

              if(id == 0.and.(j > j.and.i > i)) then 
                highest_pow_r12 = highest_pow_r12 +angularMomentum_A
              end if 

              if(id == 1.and.i > i.and.l1_prime > 0) then 
                highest_pow_r12 = highest_pow_r12 +angularMomentum_A
              end if 

              if(id == 0.and.(j > j.or.i > i)) then 
                msl = msl + 1
              end if 

              if((j > j.or.i > i)) then 
                lmax = min(angularMomentum_A +2,lmax+1)
              end if 

              if((j > j.and.i > i)) then 
                lmax = min(angularMomentum_A +2,lmax+1)
              end if 

              if(msl > 12) then 
                stop 180
              end if 

              if(exchange) then
                temp = highest_pow_r1
                highest_pow_r1 = highest_pow_r2
                highest_pow_r2 = temp
              endif

              if(j == i) then 
                write(*, *)
                write(*, '(A)') "Generating Integrals...  "
                write(*, '(A)') "--------------------------------"
                write(*, '(A, I3)') "Lowest Power of R1:   ", lowest_pow_r1
                write(*, '(A, I3)') "Lowest Power of R12:  ", lowest_pow_r12
                write(*, '(A, I3)') "Highest Power of R1:  ", highest_pow_r1
                write(*, '(A, I3)') "Highest Power of R2:  ", highest_pow_r2
                write(*, '(A, I3)') "Highest Power of R12: ", highest_pow_r12
                write(*, *)
                write(*, '(A, I2.2, 3x, A, I2.2, 3x, A, I2.2, 3x, A, I2.2)') "L1 = ", l1, "L2 = ", l2, "L1' = ", l1_prime, "L2' = ", l2_prime
                write(*, *)
                write(*, '(A, A)') "Exchange: ", trim(adjustl(merge('True ', 'False', exchange)))
                write(*, '(A, I3)') "Sum of Powers:  ", sumOfPowers ! i + j + k 
                write(*, '(A, 10X)', advance = "no") "Y11:"
                write(*, '(1p2e40.32)') y11
                write(*, '(A, 10X)', advance = "no") "Y22:"
                write(*, '(1p2e40.32)') y22

                write(4, *)
                write(4, '(A)') "Generating Integrals...  "
                write(4, '(A)') "--------------------------------"
                write(4, '(A, I3)') "Lowest Power of R1:   ", lowest_pow_r1
                write(4, '(A, I3)') "Lowest Power of R12:  ", lowest_pow_r12
                write(4, '(A, I3)') "Highest Power of R1:  ", highest_pow_r1
                write(4, '(A, I3)') "Highest Power of R2:  ", highest_pow_r2
                write(4, '(A, I3)') "Highest Power of R12: ", highest_pow_r12
                write(4, *)
                write(4, '(A, I2.2, 3x, A, I2.2, 3x, A, I2.2, 3x, A, I2.2)') "L1 = ", l1, "L2 = ", l2, "L1' = ", l1_prime, "L2' = ", l2_prime
                write(4, *)
                write(4, '(A, A)') "Exchange: ", trim(adjustl(merge('True ', 'False', exchange)))
                write(4, '(A, I3)') "Sum of Powers:  ", sumOfPowers ! i + j + k 
                write(4, '(A, 10X)', advance = "no") "Y11:"
                write(4, '(1p2e40.32)') y11
                write(4, '(A, 10X)', advance = "no") "Y22:"
                write(4, '(1p2e40.32)') y22
              end if 

              lmax = lmax + 1
              highest_pow_r12 = highest_pow_r12 + 2
              call genint(j, i, lmax, factorial, dqreal(0))
            end if
          end if 

          jdx = jdx_

          if(exchange) then 
            jdx = 3 - jdx_
          end if 

          do m = sectorStart_A, sectorStop_A !   SUMMATION OVER RIGHT HAND STATE BASIS FUNCTIONS.
            do n = sectorStart_B, sectorStop_B !   SUMMATION OVER LEFT HAND STATE BASIS FUNCTIONS.
              if(nuclearCharge_B == 0 .and. n < m) then 
                cycle
              end if 

              psi2 = waveFnCoefficients_A(m) * waveFnCoefficients_B(n) * dqsign * dqreal(2)

              if(nuclearCharge_B == 0.and.n /= m) then
                psi2 = psi2 * dqreal(2)
              end if 

              R1 = powers_R1_and_R2_A(ide, m) + powers_R1_and_R2_B(jde, n)
              R2 = powers_R1_and_R2_A(idep, m) + powers_R1_and_R2_B(jdep, n)
              R12 = powers_R12_A(m) + powers_R12_B(n)

              !   1/R1**2
              sum = psi2 * (spl(R1-2,R2,R12,1) + spl(R1,R2-2,R12,1)) / dqreal(2)
              matrix_element(1) = matrix_element(1) + sum
              
              if(lgo) then 
                screened_hydrogenic(1) = screened_hydrogenic(1) + sum
              end if 

              !   1/R1
              sum = psi2 * (spl(R1-1,R2,R12,1) + spl(R1,R2-1,R12,1)) / dqreal(2)
              matrix_element(4) = matrix_element(4) + sum

              if(lgo) screened_hydrogenic(2) = screened_hydrogenic(2) + sum

              !   R1
              sum = psi2 * (spl(R1+1,R2,R12,1) + spl(R1,R2+1,R12,1)) / dqreal(2)
              matrix_element(7) = matrix_element(7) + sum

              if(lgo) then 
                screened_hydrogenic(3) = screened_hydrogenic(3) + sum
              end if 

              !   R1**2
              sum = psi2 * (spl(R1+2,R2,R12,1) + spl(R1,R2+2,R12,1)) / dqreal(2)
              matrix_element(10) = matrix_element(10) + sum

              if(lgo) then 
                screened_hydrogenic(4) = screened_hydrogenic(4) + sum
              end if 

              !   1/R12**2
              sum = psi2 * spl(R1,R2,R12-2,1)
              matrix_element(13) = matrix_element(13) + sum

              if(lgo) then
                screened_hydrogenic(5) = screened_hydrogenic(5) + sum
              end if 

              !   1/R12
              sum = psi2 * spl(R1,R2,R12-1,1)
              matrix_element(16) = matrix_element(16) + sum

              if(lgo) then 
                screened_hydrogenic(6) = screened_hydrogenic(6) + sum
              end if 

              !   R12
              sum = psi2 * spl(R1,R2,R12+1,1)
              matrix_element(19) = matrix_element(19) + sum

              if(lgo) then 
                screened_hydrogenic(7) = screened_hydrogenic(7) + sum
              end if 

              !   R12**2
              sum = psi2 * spl(R1,R2,R12+2,1)
              matrix_element(22) = matrix_element(22) + sum

              if(lgo) then 
                screened_hydrogenic(8) = screened_hydrogenic(8) + sum
              end if 

              !   R1.R2
              sum = psi2 * spl(R1+1,R2+1,R12,2)
              matrix_element(25) = matrix_element(25) + sum

              if(lgo) then 
                screened_hydrogenic(9) = screened_hydrogenic(9) + sum
              end if 

              !   1/(R1R12)
              sum = psi2 * (spl(R1-1,R2,R12-1,1)+spl(R1,R2-1,R12-1,1)) / dqreal(2)
              matrix_element(28) = matrix_element(28) + sum

              if(lgo) then 
                screened_hydrogenic(10) = screened_hydrogenic(10) + sum
              end if 

              !   1/(R1R2)
              sum = psi2 * spl(R1-1,R2-1,R12,1)
              matrix_element(31) = matrix_element(31) + sum

              if(lgo) then 
                screened_hydrogenic(11) = screened_hydrogenic(11) + sum
              end if 

            end do
          end do
        end do
      end do
    end do
  end do
  return
end subroutine 
