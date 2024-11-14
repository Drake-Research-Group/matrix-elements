program h1  
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
  
    type(dq_real) :: da0, db0
    type(dq_real) :: dqsign
    type(dq_real) :: reducedMassRatio_A, reducedMassRatio_B
    type(dq_real), dimension(40) :: matrix_element, screened_hydrogenic
  
    ! Determine OS in order to grab path for calling correct C++ program for Unicode Symbols
    character(len = 10) :: char_integer
    character(len = 20) :: OS_NAME
    character(len = 50) :: build_path
    character(len = 50) :: sys_call
  
    ! Variable Initializations
    matrix_element = dqreal(0)
    screened_hydrogenic = dqreal(0)
    iostat = 0
    ih = 2
    da0 = dqreal(1)
    db0 = dqreal(1)
    screen = .false.
    reducedMassRatio_A = dqreal(0)
  
    ! Grab system's Operating System and store it in a txt
    call execute_command_line('python3 ./scripts/getOS.py')
  
    ! Determine the pathing to the right build folder depending on detected OS
    open(1, file='data/operating-system.txt', iostat = ios)
    read(1, '(A)') OS_NAME
    close(1)
  
    if(trim(OS_NAME) == "macOS") then 
      build_path = "build/macOS/"
    else if(trim(OS_NAME) == "Linux") then
      build_path = "build/linux/"
    else if(trim(OS_NAME) == "Windows") then
      build_path = ".\build\windows\"
    else 
      write(*, '(A)') "Unknown Operating System, closing the program."
      stop
    end if 
  
    ! Open matl88.dat for reading and matl88.out for writing
    call routedOpen(5,FILE='matl88.dat',STATUS='UNKNOWN')
    call routedOpen(4,FILE='matl88.out',STATUS='UNKNOWN')
    
    read(5,5) line(1)
    read(5,*) matout
  
    write(*,5) line(1)
    write(4,5) line(1)
  
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
  
      call dqread(1, waveFnEnergies_A(1), waveFnEnergies_A(2), reducedMassRatio_A)
      read(1, '(I5)', advance = 'no') numNonLinParams_A
  
  
      ! Now that the number of terms is known, we can allocate the proper sizes to the arrays pertaining to waveFn_A
      allocate(waveFnSectorStart_A(numNonLinParams_A + 1))
      allocate(orbitalAngMomentum_A(2, numTerms_A))
      allocate(powers_R1_and_R2_A(2, numTerms_A))
      allocate(waveFnCoefficients_A(numTerms_A))
      allocate(nonLinParams_A(2, numTerms_A))
      allocate(powers_R12_A(numTerms_A))
  
      nonLinParams_A(1, 1) = dqreal(0)
      nonLinParams_A(2, 1) = dqreal(0)
      waveFnSectorStart_A(1) = 0
  
      read(1, '(15I5)')  (waveFnSectorStart_A(i + 1), i = 1, numNonLinParams_A), (maxPartialWaveFnSize_A(i + 1), i = 1, numNonLinParams_A - ih + 1)
  
  
      ! Print important wavefunction data to screen
      write(*, *)
      write(*, '(A)') waveFnName_A
      write(*, '("Z = ", i3.3, 3x, "L = ", i3.3, 3x, "S = ", i3.3, 3x, "N = ", i4.4)') nuclearCharge_A, angularMomentum_A, totalSpin_A, numTerms
      write(*, *)
      write(*, '(A, 5x)', advance = "no") "Screened Energy:"
      call dqwrite(6, waveFnEnergies_A(1))
      write(*, '(A)', advance = "no") "Wavefunction Energy: " 
      call dqwrite(6, waveFnEnergies_A(2)) 
      write(*, '(A)', advance = "no") "Reduced Mass Ratio:  "
      call dqwrite(6, reducedMassRatio_A)
      
      if(reducedMassRatio_A /= dqreal(0)) then
        dumpfile_type = '.DML' ! Finite Mass Case
      else 
        dumpfile_type = '.DMP' ! Infinite Mass Case
      end if 
  
      ! Check if DumpFile exists:
      dumpfile_path = trim('dump/' // waveFnName_A // '.DMP')
      inquire(file = dumpfile_path, exist = dumpfile_exists)
      
      ! If DumpFile doesn't exist, use the coefficients stored in the wavefunction instead
      if(.not.dumpfile_exists) then 
        write(*, '(A)') "Dumpfile not found. Using coefficients from wavefunction file instead."
        dumpfile_type = ''
      end if
  
      do i=1,numNonLinParams_A
        waveFnSectorStart_A(i+1) = waveFnSectorStart_A(i + 1) + eigenstateZero
  
        read(1, *) sectorStart, sectorStop, orbitalAngMomentum_A(1, sectorStart + eigenstateZero), orbitalAngMomentum_A(2, sectorStart + eigenstateZero)
        call dqread(1, nonLinParams_A(1, sectorStart + eigenstateZero))
        call dqread(1, nonLinParams_A(2, sectorStart + eigenstateZero))
        write(*, *)
        write(*, '(A, I3)') "Basis Set Sector:     ", i
        write(*, '(A)', advance = "no") "Corresponding Alpha: "
        call dqwrite(6, nonLinParams_A(1, sectorStart + eigenstateZero))
        write(*, '(A)', advance = "no") "Corresponding Beta:  "
        call dqwrite(6, nonLinParams_A(2, sectorStart + eigenstateZero))
  
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
      end do
  
      read(1,'(10(i3,2i2))') (powers_R1_and_R2_A(1,k), powers_R1_and_R2_A(2,k), powers_R12_A(k),k=eigenstate_A, numTerms_A)
  
      if(dumpfile_type == ".DMP") then 
        open(unit = 1000, file = 'dump/' // wavefnName_A(:len(trim(waveFnName_A))-4) // '.DMP', form='unformatted')
        read(1000) (waveFnCoefficients_A(i), i = eigenstate_A, numTerms_A)
        close(1000)
      else if(dumpfile_type == ".DML") then
        open(unit = 1000, file = 'dump/' // wavefnName_A(:len(trim(waveFnName_A))-4) // '.DML',  form='unformatted')
        read(1000) (waveFnCoefficients_A(i), i = eigenstate_A, numTerms_A)
        close(1000)
      else 
        do i = eigenstate_A, numTerms_A
          call dqread(1, waveFnCoefficients_A(i))
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
  
        do m = 1, eigenstate_A
          orbitalAngMomentum_A(:, m) = orbitalAngMomentum_A(:, eigenstate_A)
          nonLinParams_A(:, m) = nonLinParams_A(:, 1)
          powers_R1_and_R2_A(1, m) = increasedOrbitalAngMomentum
          powers_R1_and_R2_A(2, m) = angularMomentum_A + m - 1
          powers_R12_A(m) = 0
          waveFnCoefficients_A(m) = waveFnCoefficients_A(eigenstate_A) * dqreal(1) * dqreal(2)**(increasedOrbitalAngMomentum + angularMomentum_A + 2) * (dqreal(nuclearCharge_A) - dqreal(1))**(angularMomentum_A + 1) * sqrt((dqreal(nuclearCharge_A - 1)) &
                  & * factorial( angularMomentum_A + eigenstate_A + angularMomentum_A + 1) / (dqreal(nuclearCharge_A) * factorial(angularMomentum_A + eigenstate_A - angularMomentum_A) &
                  & * factorial(2 * increasedOrbitalAngMomentum + 2))) / (factorial(2 * angularMomentum_A + 2) * dqreal( angularMomentum_A + eigenstate_A) * (dqreal(increasedOrbitalAngMomentum) + dqreal(1))**(increasedOrbitalAngMomentum + 2) * (dqreal( angularMomentum_A + eigenstate_A) * dqreal(nuclearCharge_A))**(angularMomentum_A + 1))
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
          stop
        endif
  
        read(1,'(I5)', advance = 'no') numNonLinParams_B
  
        ! NumTerms_B is now known, so now we can allocate all array sizes for waveFn B
        allocate(waveFnSectorStart_B(numNonLinParams_B + 1))
        allocate(orbitalAngMomentum_B(2, numTerms_B))
        allocate(powers_R1_and_R2_B(2, numTerms_B))
        allocate(waveFnCoefficients_B(numTerms_B))
        allocate(nonLinParams_B(2, numTerms_B))
        allocate(Powers_R12_B(numTerms_B))
  
        read(1, '(15I5)') (waveFnSectorStart_B(i+1),i=1,numNonLinParams_B),(maxPartialWaveFnSize_A(i+1),i=1,numNonLinParams_B-ih+1)
        
        screen = .false.
  
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
        end do
  
        read(1,'(10(i3,2i2))')(powers_R1_and_R2_B(1,k),powers_R1_and_R2_B(2,k),powers_R12_B(k),k=eigenstate_B,numTerms_B)
  
        do i = eigenstate_B, numTerms_B
          call dqread(1, waveFnCoefficients_B(i))
        end do 
  
        write(4,'("Z = ",i3.3, 3x,"L = ",i3.3, 3x,"S = ",i3.3, 3x,"n = ",i4.4)') nuclearCharge_B, angularMomentum_B, totalSpin_B, numTerms
        write(*,'("Z = ",i3.3, 3x,"L = ",i3.3, 3x,"S = ",i3.3, 3x,"n = ",i4.4)') nuclearCharge_B, angularMomentum_B, totalSpin_B, numTerms
  
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
          end do 
        end if 
  
      else
        angularMomentum_B = angularMomentum_A
        totalSpin_B = totalSpin_A
        eigenstate_B = eigenstate_A
        waveFnName_B = 'SAME'
        numTerms_B = numTerms_A
        numNonLinParams_B = numNonLinParams_A
        
  
        ! Allocate Array sizes using numTerms_B
        allocate(waveFnSectorStart_B(numNonLinParams_B + 1))
        allocate(orbitalAngMomentum_B(2, numTerms_B))
        allocate(powers_R1_and_R2_B(2, numTerms_B))
        allocate(waveFnCoefficients_B(numTerms_B))
        allocate(nonLinParams_B(2, numTerms_B))
        allocate(Powers_R12_B(numTerms_B))
  
        do i=1,numNonLinParams_B
          waveFnSectorStart_B(i+1) = waveFnSectorStart_A(i+1)
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
      open(1000, file='mat-output.txt')
  
      write(*, '(A)') "Matrix Elements:"
      call execute_command_line("python3 ./scripts/prettify-output.py")
  
      do i = 1, 31, 3
        write(char_integer, '(I0)') (i + 2)/3
        sys_call = trim(build_path)//'print-equations '//char_integer
        call execute_command_line(sys_call)
        call dqwrite(6, screened_hydrogenic((i + 2)/3))
        call dqwrite(6, matrix_element(i))
        print*, " "
        call dqwrite(1000, screened_hydrogenic((i + 2)/3))
        call dqwrite(1000, matrix_element(i))
        write(1000, *) " "
      end do
      
      ! Calculation is finished, deallocate arrays to prevent memory leaks
      deallocate(powers_R12_A)
      deallocate(powers_R12_B)
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
  
  
      call routedOpen(1, file = matout, status = 'OLD')
  
      backspace 1
  
      write(1,'(A14,1X,A14)') waveFn_A,waveFn_B
  
      do i = 1, 31, 3
        call dqwrite(1, matrix_element(i))
      end do 
  
      write(*,5) line(1)
      write(*,5)
  
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
    5  format(a51)  
end program h1


subroutine breit(matrix_element,screened_hydrogenic)
   !? Note: x is an overwritten variable from the common blocks in a1a. Because of this there were array rank mismatch arguments.
   !? To keep it simple the 'x' real was renamed to 'x_'. This prevents the program from trying to write and compare the 'x' array instead of the 'x' real value.
   !? This problem only exists within this function or any other function with this comment. All other subroutines or functions are fine.
   !? This problem exists from the renaming of common block variables from an older technique. This is the quick fix. A better fix should be implemented later
   !* Evan Petrimoulx June 13 2023

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
    use format 
    use wavExt
    implicit none

    logical :: lgo, exchange
    
    character(len = 80) :: buffer
    
    type(dq_real), dimension(40) :: matrix_element, screened_hydrogenic
    type(dq_real) :: psi2, dqsign, sum, test
    type(dq_real) :: a1, a4, b1, b4, c1, c4, d2d1, de, del1, del2, delt12, e, e0, fac, fp
    type(dq_real) :: ov, ovmax, s1, s2, s3, sgn1, sgn2, t1, t10, t23, tt, x_, xx, y1, y2, y3, y4, y_, yaa, yaa2, ybb, ybb2, yp_, z, zz1, zz2

    integer :: p, q, s
    integer :: idx, jdx_, i, j, m, n ! Loop iterators
    integer :: id, jd, idep, jdep
    integer :: sectorStart_A, sectorStop_A
    integer :: totalAngMomentum_A, totalAngMomentum_B
    integer :: lmax, temp, sectorStart_B, sectorStop_B
    integer :: iaa, ibb, ic, ic2, ic3, imax, ip, ipow, iq, ir1, is, isum, isumd, itl, itr, jdxx, kang, kex, kkl, kkr, l13, l3, l4, ll, ll1, ll3
    integer :: ir2, jmax, ll4, lp2, lp3, lss, lsuma, lsumb, md, mla, mlb, mlp2, mlp3, mp2d, mq2d, ms2d, msld, mtemp, nbl1, nbl3, nbx1, nbx3, nbx4 
    integer :: nbxa, nbxb, nhi, nhj, nn, nbx2

    integer, save :: sum_d, highest_pow_r1d, highest_pow_r2d, highest_pow_r12d

    !? FOR S STATES
    !? DATA ISUMD,MP2D,MQ2D,MS2D,MSLD/0,1,2,-2,1/
    !? FOR P STATES
    !? DATA ISUMD,MP2D,MQ2D,MS2D,MSLD/0,2,2,0,2/
    !? FOR D, F, G STATES
    data isumd,mp2d,mq2d,ms2d,msld/0,2,2,0,2/

    !? POWERS OF R1 AND R12.  POWER OF R2 = 0.
    data ipow/-2,0, -1,0, 0,0, 0,-2, 1,-2 ,0,-1/
    type(dq_real), save :: xc 
    
    xc = dqreal(0)

    !? IST CHOOSES THE COORDINATES OF ELECTRON 1 OR 2 AND 'I' SPECIFIES WHICH
    !? SET OF COEFFICIENTS OF PL(COSTHETA12) IS TO BE USED.
    !? SET KROSS = 5 OR 9 FOR USE IN 'CROSS'.

    machineError = dqreal(1.1)

    do while(test .gt. 1.1q0)
        machineError = machineError/dqreal(2)
        test = dqreal(1.1) + machineError
    end do 

    kross = 5
    klog = 0
    nuclearCharge_B = 0
    kang = 1
    kang = 0
    mla = angularMomentum_A
    mlb = angularMomentum_B
    md = - mla + mlb

    do i=1,10
        matrix_element(i) = dqreal(0)
        screened_hydrogenic(i) = dqreal(0)
    end do 

    zz1 = (z-1.q0)/z
    zz2 = zz1*zz1
    e = wa(2)/z**2
    de = 0.5q0*wa(1)
    nn = angularMomentum_A + eigenstate_A
    e0 = - 0.5q0 - 0.5q0*zz2/nn**2
    t10 = (2.q0/(nn**3*(angularMomentum_A+0.5q0)) - 1.5q0/nn**4)*zz1**4
    ar(6) = dqreal(0)
    al(6) = dqreal(0)
    xc(3,3) = dqreal(2)*e*e - 2.5q0
    xc(2,3) = dqreal(4)*e
    xc(3,2) = dqreal(4)*e
    xc(3,6) = dqreal(-4)*e/z
    xc(1,3) = dqreal(2)
    xc(2,2) = dqreal(4)
    xc(3,1) = dqreal(2)
    xc(2,6) = dqreal(-4)/z
    xc(6,2) = dqreal(-4)/z
    xc(6,6) = dqreal(2)/z**2

    !? SUMMATION OVER RIGHT HAND STATE DIRECT AND EXCHANGE TERMS.
    do idx=1,2

        !? SUMMATION OVER BLOCKS OF TERMS WITH THE SAME NON-LINEAR PAREMETERS.
        do nhi=1,nbxa
            nbx1 = numNonLinParams_A(nhi) + 1
            nbx2 = numNonLinParams_A(nhi+1)
            do nhj=1,nbxb
                nbx3 = numNonLinParams_B(nhj) + 1
                nbx4 = numNonLinParams_B(nhj+1)

                !   THE OUTER LOOP CORRESPONDS TO THE RIGHT HAND STATE AND IS DENOTED BY
                !   'A' OR'1' OR '2' IN VARIOUS PLACES.
                do nbl1=nbx1,nbx2
                    sectorStart_A = waveFnSectorStart_A(nbl1) + 1
                    sectorStop_A = waveFnSectorStart_A(nbl1+1)

                    !   THE INNER LOOP CORRESPONDS TO THE LEFT HAND STATE AND IS DENOTED BY
                    !   'B' OR '3' OR '4' IN VARIOUS PLACES (3 MEANS ELECTRON 1 AND 4 MEANS ELECTRON 2).
                    do nbl3=nbx3,nbx4
                        sectorStart_B = waveFnSectorStart_B(nbl3) + 1
                        sectorStop_B = waveFnSectorStart_B(nbl3 + 1)
                        lgo = .false.

                        if(nbl1.eq.1.and.nbl3.eq.1.and.z.ne.1.) then 
                            lgo = .true.
                        end if 

                        if(nuclearCharge_B.eq.0.and.nbl3.lt.nbl1) then 
                            cycle
                        end if 

                        !   SUMMATION OVER LEFT HAND STATE DIRECT AND EXCHANGE TERMS.
                        !   THE ORDERING IS D-D, E-E, D-E AND E-D.
                        ovmax = dqreal(0)
                        jdx_ = 1
                        ide =  jdx_
                        jde =  jdx_

                        if(idx.eq.2) then 
                            ide = 3 - jde
                        end if 

                        id = ide - 1
                        jd = jde - 1
                        idep = 3 - ide
                        jdep = 3 - jde
                        sgn1 = dqreal(1)

                        if(totalSpin_A.eq.1) then 
                            sgn1 = (-1.q0)**(ide-1)
                        end if 

                        sgn2 = sgn1

                        if(totalSpin_B.eq.1) then 
                            sgn2 = sgn1*(-1.)**(jde-1)
                        end if 

                        l1 = la(ide,sectorStart_A)
                        l2 = la(idep,sectorStart_A)
                        l3 = lb(jde,sectorStart_B)
                        l4 = lb(jdep,sectorStart_B)
                        ll1 = l1*(l1+1)
                        ll3 = l3*(l3+1)
                        ll4 = l4*(l4+1)
                        l13 = l1 + l3

                        call cross(l3,l4,l1,l2,angularMomentum_B,mlb,angularMomentum_A,mla,ide,jde)

                        mlp2 = mlt(2)
                        mlp3 = mlt(3)
                        lp2 = lt(2)
                        lp3 = lt(3)
                        mlt(6) = mlp2
                        mlt(9) = mlp2
                        mlt(11) = mlp3
                        mlt(12) = mlp3
                        mlt(13) = mlp3

                        if(kang.ne.0) then 
                            write(4,fmt3) lss,l3,l4,l1,l2,angularMomentum_B,mlb,angularMomentum_A,mla,ide,jde,sectorStart_B,sectorStop_B,sectorStart_A,sectorStop_A,nuclearCharge_B

                            write(4,fmt4)

                            do ic=1,16
                                if(lt(ic).gt.0) then 
                                    ic2 = lt(ic)

                                    write(4,fmt4) ic,(cplt(ic3,ic), ic3=  1,ic2)
                                end if 
                            end do 

                        end if 

                        if((nbl1.gt.nbx1.or.nbl3.gt.nbx3) .or. jdx_ .eq. 2) then 
                            go to 77
                        end if 

                        y11 = nonLinParams_A(ide,sectorStart_A) + nonLinParams_B(jde,sectorStart_B)
                        y22 = nonLinParams_A(idep,sectorStart_A) + nonLinParams_B(jdep,sectorStart_B)
                        kex = 0

                        if(y22.gt.y11) then
                            y22 = nonLinParams_A(ide,sectorStart_A) + nonLinParams_B(jde,sectorStart_B)
                            y11 = nonLinParams_A(idep,sectorStart_A) + nonLinParams_B(jdep,sectorStart_B)
                            kex = 1
                        end if

                        lmax = lt(3)
                        lowest_pow_r1 = 3
                        lowest_pow_r12 = 3
                        msl = lt(1) + msld
                        imin = -2 + angularMomentum_A + angularMomentum_B
                        lsuma = 2*id
                        lsumb = l2*(iabs(l1-l2)/2)
                        totalAngMomentum_A = angularMomentum_A + eigenstate_A - 1
                        totalAngMomentum_B = angularMomentum_B + eigenstate_B - 1

                        if((sectorStart_A.ne.1.and.sectorStart_B.ne.1).or.z.le.1.) then 
                            isum = ijkb(nbl3) + ijka(nbl1) + isumd
                            highest_pow_r1 = iib(nbl3) + (1-id)*iia(nbl1) + id*jja(nbl1) + lowest_pow_r1 + mp2d + 3*(1 - id) + id
                            highest_pow_r2 = jjb(nbl3) + (1-id)*jja(nbl1) + id*iia(nbl1)+ lowest_pow_r1+mq2d - 1 + id + id
                            highest_pow_r12 = lowest_pow_r12 + kkb(nbl3) + kka(nbl1) + 2*lt(1) + ms2d

                            go to 80
                        end if 

                        if(sectorStart_B.ne.1.or.sectorStart_A.ne.1) then 
                            if(sectorStart_B.ne.1) then 
                                isum = totalAngMomentum_A + ijkb(nbl3) + isumd
                                highest_pow_r1 = iib(nbl3) + id*totalAngMomentum_A + lowest_pow_r1 + mp2d
                                highest_pow_r2 = jjb(nbl3) + (1-id)*totalAngMomentum_A + lowest_pow_r1 + mq2d - 1
                                highest_pow_r12 = lowest_pow_r12 + kkb(nbl3) + 2*lt(1) + ms2d

                                if(id.eq.1.and.nbx4.gt.nbx3) then 
                                    highest_pow_r1 = highest_pow_r1 +1
                                end if 

                                go to 80
                            end if 

                            isum = totalAngMomentum_B + ijka(nbl1) + isumd
                            highest_pow_r1 = (1-id)*iia(nbl1) + id*(jja(nbl1)+totalAngMomentum_B) + lowest_pow_r1 + mp2d
                            highest_pow_r2 = (1-id)*(jja(nbl1)+totalAngMomentum_B) + id*iia(nbl1) + lowest_pow_r1 + mq2d - 1
                            highest_pow_r12 = lowest_pow_r12 + kka(nbl1) + 2*lt(1) + ms2d

                            go to 80
                        end if 

                        isum = totalAngMomentum_B + totalAngMomentum_A + isumd
                        highest_pow_r1 = (1-id)*totalAngMomentum_A + id*(totalAngMomentum_B+totalAngMomentum_A) + lowest_pow_r1 + mp2d
                        highest_pow_r2 = (1-id)*(totalAngMomentum_B+totalAngMomentum_A) + id*totalAngMomentum_A + lowest_pow_r1 + mq2d
                        highest_pow_r12 = lowest_pow_r12 + 2*lt(1) + ms2d

                        80 continue

                        if(id.eq.0.and.(nbx4.gt.nbx3.or.nbx2.gt.nbx1)) then
                            highest_pow_r12 = highest_pow_r12 + 2
                        end if 

                        if(id.eq.0.and.(nbx4.gt.nbx3.and.nbx2.gt.nbx1)) then
                            highest_pow_r12 = highest_pow_r12 + 2
                        end if 

                        if(id.eq.1.and.nbx2.gt.nbx1.and.l3.gt.0) then
                            highest_pow_r12 = highest_pow_r12 + 2
                        end if 

                        if(id.eq.0.and.(nbx4.gt.nbx3.or.nbx2.gt.nbx1)) then
                            msl = msl+1
                        end if 

                        if((nbx4.gt.nbx3.or.nbx2.gt.nbx1)) then
                            lmax =min(angularMomentum_A +3,lmax+1)
                        end if 

                        if((nbx4.gt.nbx3.and.nbx2.gt.nbx1)) then
                            lmax =min(angularMomentum_A +3,lmax+1)
                        end if 

                        if(msl.gt.12) then 
                            stop 180
                        end if 

                        if(kex.eq.1) then
                            mtemp = highest_pow_r1
                            highest_pow_r1 = highest_pow_r2
                            highest_pow_r2 = mtemp
                        end if

                        if(nbl3.eq.nbl1) then 
                            write(*, fmt5) kex,lowest_pow_r1,highest_pow_r1,highest_pow_r2,lowest_pow_r12,highest_pow_r12,isum,y11,y22,l3,l4,l1,l2,nbl1,nbl3
                        end if 

                        call genint(id,nbl3,nbl1,lmax,fac,0.q0)

                    77 delt12 = dqreal(0)
                        jdxx = jdx_

                        if(kex.eq.1) then 
                            jdxx = 3 - jdx_
                        end if 

                        ! SUMMATION OVER RIGHT HAND STATE BASIS FUNCTIONS.
                        do i=sectorStart_A,sectorStop_A
                            a1 = powers_R1_and_R2_A(ide,i)
                            b1 = powers_R1_and_R2_A(idep,i)
                            c1 = powers_R12_A(i)         
                            y1 = nonLinParams_A(ide,i)
                            y2 = nonLinParams_A(idep,i)
                            ar(1)= a1*(a1+1) - l1*(l1+1)
                            ar(2)= -2.*y1*(a1+1)
                            ar(3) = y1*y1
                            ar(4) = c1*(c1+1) + 2*a1*c1
                            ar(5) = -2.*y1*c1

                            if(c1.ne.0) then 
                                do ll=mlp2,lp2,2
                                    cplt(ll,9) = a1*cplt(ll,2) + cplt(ll,10)
                                end do 

                                lt(9) = lp2

                                if(abs(cplt(lp2,9)).lt.1.d-06) then 
                                    lt(9) = lp2 - 2
                                end if 

                                do ll=mlp3,lp3,2
                                    cplt(ll,11) = a1*cplt(ll,3) + cplt(ll,8)
                                end do 

                                lt(11) = lp3

                                if(abs(cplt(lp3,11)).lt.1.d-06) then 
                                    lt(11) = lp3 - 2
                                end if 
                            end if 

                            ! SUMMATION OVER LEFT HAND STATE BASIS FUNCTIONS.
                            do j=sectorStart_B,sectorStop_B

                                if(nuclearCharge_B.eq.0.and.j.lt.i) then 
                                    cycle
                                end if 

                                psi2 = 2.*da(i)*db(j)*sgn2

                                if(nuclearCharge_B.eq.0.and.j.ne.i) then 
                                    psi2 = psi2* dqreal(2)
                                end if 

                                a4 = powers_R1_and_R2_B(jde,j)
                                b4 = powers_R1_and_R2_B(jdep,j)
                                c4 = powers_R12_B(j)
                                y3 = nonLinParams_B(jde,j)
                                y4 = nonLinParams_B(jdep,j)
                                p = powers_R1_and_R2_A(ide,i) + powers_R1_and_R2_B(jde,j)
                                q = powers_R1_and_R2_A(idep,i) + powers_R1_and_R2_B(jdep,j)
                                s = powers_R12_A(i) + powers_R12_B(j)
                                y_ = y1 + y3
                                yp_ = y2 + y4

                                if(c4.ne.0) then 
                                    do ll=mlp2,lp2,2
                                      cplt(ll,6) = b4*cplt(ll,2) + cplt(ll,5)
                                    end do 

                                    lt(6) = lp2

                                    if(abs(cplt(lp2,6)) .lt. 1.q-06) then 
                                      lt(6) = lp2 - 2
                                    end if 

                                    do ll=mlp3,lp3,2
                                      cplt(ll,12) = b4*cplt(ll,3) + cplt(ll,7)
                                      cplt(ll,13) = a1*b4*cplt(ll,3) + a1*cplt(ll,7) + b4*cplt(ll,8) + cplt(ll,4)
                                    end do 

                                    lt(12) = lp3
                                    lt(13) = lp3

                                    if(abs(cplt(lp3,12)).lt.1.q-06) then 
                                      lt(12) = lp3 - 2
                                    end if 

                                    if(abs(cplt(lp3,13)).lt.1.q-06) then 
                                      lt(13) = lp3 - 2
                                    end if 
                                end if 

                                al(1)= b4*(b4+1) - l4*(l4+1)
                                al(2)= dqreal(-2)*y4*(b4+1)
                                al(3) = y4*y4
                                al(4) = c4*(c4+1) + 2*b4*c4
                                al(5) = dqreal(-2)*y4*c4

                                if(l13+s.ne.0) then 
                                        
                                    sum = dqreal(0)
                                    ov = spl(p,q,s,1)
                                    kkr = -1

                                    do itr=1,6
                                        kkr = kkr + 2
                                        x_ = ar(itr)
                                        ir1 = p + ipow(kkr)
                                        ir2 = q
                                        ir3 = s + ipow(kkr+1)
                                        kkl = -1

                                        do itl=1,6
                                            kkl = kkl + 2
                                            xx = x_*al(itl) - xc(itr,itl)

                                            if(abs(xx).ge.1.q-10) then 
                                                ip = ir1
                                                iq = ir2 + ipow(kkl)
                                                is = ir3 + ipow(kkl+1)
                                                sum = sum + xx*spl(ip,iq,is,1)
                                            end if 
                                        end do 
                                    end do 

                                    if(c4.ne.0) then 
                                        kkr = -1

                                        do itr=1,5
                                            kkr = kkr + 2
                                            x_ = dqreal(-2)*c4*ar(itr)

                                            if(x_.ne. dqreal(0)) then 
                                                ip = p + ipow(kkr) + 1
                                                iq = q - 1
                                                is = s + ipow(kkr+1) - 2
                                                sum = sum + x_*(spl(ip,iq,is,6)- y4*spl(ip,iq+1,is,2))
                                            end if 
                                        end do 
                                    end if 

                                    if(c1.ne.0) then 
                                      kkl = -1

                                      do itl=1,5
                                        kkl = kkl + 2
                                        x_ = dqreal(-2)*c1*al(itl)

                                        if(x_.ne.dqreal(0)) then 
                                          ip = p - 1
                                          iq = q + ipow(kkl) + 1
                                          is = s + ipow(kkl+1) - 2
                                          sum = sum + x_*(spl(ip,iq,is,9) - y1*spl(ip+1,iq,is,2))
                                        end if 
                                      end do 
                                    end if 

                                    x_ = dqreal(4)*c1*c4
                                    is = s - 4

                                    if(x_.ne.dqreal(0)) then 
                                        sum = sum + x_*(spl(p,q,is,13) - y1*spl(p+1,q,is,12) - y4*spl(p,q+1,is,11) + y1*y4*spl(p+1,q+1,is,3))
                                    end if 

                                    tt = sum*psi2
                                    matrix_element(1) = matrix_element(1) - tt

                                    if(lgo) then 
                                        screened_hydrogenic(1) = screened_hydrogenic(1) - tt
                                    end if 

                                    matrix_element(4) = matrix_element(4) + ov*psi2
                                    x_ = abs(ov*psi2)

                                    if(x_.gt.ovmax) then 
                                        ovmax = x_
                                        imax = i
                                        jmax = j
                                    end if 

                                    if(lgo) then 
                                        screened_hydrogenic(4) = screened_hydrogenic(4) + ov*psi2
                                    end if 

                                    cycle
                                end if 

                                ! SPECIAL CODING FOR UNCORRELATED TERMS.
                                iaa = a1 + a4 + 2
                                ibb = b1 + b4 + 2

                                if(jdxx.ne.1) then 
                                    ip = p
                                    p = q
                                    q = ip
                                end if 

                                yaa = y_/dqreal(iaa)
                                ybb = yp_/dqreal(ibb)
                                yaa2 = y_*yaa/dqreal(iaa-1)
                                ybb2 = yp_*ybb/dqreal(ibb-1)
                                fp = at(jat(p+lowest_pow_r1,q+lowest_pow_r1,1)+s+lowest_pow_r12)/dqreal(2)

                                if(lgo.or.angularMomentum_A.eq.0) then 
                                    del1 = ar(1)*yaa2 + ar(2)*yaa + ar(3)
                                    del2 = al(1)*ybb2 + al(2)*ybb + al(3)
                                    d2d1 = del1*del2
                                    t1 = dqreal(2)*((e0 + yaa + zz1*ybb)**2 + yaa2/dqreal(iaa) + zz2*ybb2/dqreal(ibb) - 1.25q0) - d2d1

                                    go to 42
                                end if 

                                t1 = t10
                             42 t1 = t1*fp
                                t23 = dqreal(2)*fp*de*(dqreal(2)*(e0 + yaa + zz1*ybb) + de)
                                s1 = dqreal(4)*cplt(1,1)*(e*(tx(p+1-kex,q+kex,s+1,2-kex) -  tx(p+1,q+1,s,2-kex)) + (tx(p,q,s+1,2-kex) -  tx(p+kex,q+1-kex,s,2-kex)))/z
                                s2 = dqreal(2)*(tx(p+1-2*kex,q-1+2*kex,s+1,2-kex) - tx(p+1-kex,q+kex,s,2-kex))/z

                                ! S3 = 2*<1/R**2 - 1/R2**2>
                                if(kex.eq.0) then 
                                    s3 = at2(p+lowest_pow_r1,q+lowest_pow_r1)/(z*z)
                                end if 

                                if(kex.eq.1) then 
                                    s3 = (at(jat(p+lowest_pow_r1,q+lowest_pow_r1,1)+s+lowest_pow_r12-2) -at(jat(p+lowest_pow_r1-2,q+lowest_pow_r1,1)+s+lowest_pow_r12))/(z*z)
                                end if 

                                tt = (s1 + s2 + s3 + t1 + t23)*psi2
                                matrix_element(1) = matrix_element(1) + tt

                                if(lgo) then 
                                    screened_hydrogenic(1) = screened_hydrogenic(1) + tt
                                end if 

                                matrix_element(4) = matrix_element(4) + fp*psi2
                                x_ = abs(fp*psi2)

                                if(x_.gt.ovmax) then 
                                    ovmax = x_
                                    imax = i
                                    jmax = j
                                end if 

                                if(lgo) then 
                                    screened_hydrogenic(4) = screened_hydrogenic(4) + fp*psi2
                                end if 
                            end do 
                        end do 
                    end do 
                end do 
            end do 
        end do 
    end do 
    
    return
end subroutine breit
