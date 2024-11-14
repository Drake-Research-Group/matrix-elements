include 'cross1.f90'
include 'spin68.f90'

program del1

  !? Import CommonBlocks
  use a1a_block
  use b1b_block
  use f1f_block
  use h1h_block
  use maxpow_block

  !? Import Formatting
  use format

  !? Import File Management for POW and POL files.
  use wavExt

  implicit real(16) (a-h,o-z)

  real(16), dimension(10) :: hterm, htscr
  integer, dimension(50) :: nn

  character(len = 57) :: lines(100)
  character(len = 16) :: matout
  character(len = 52) :: title
  character(len = 12) :: nama, namb, fmt3a(3)
  character(len = 51) :: line
  character(len = 5)  :: eigen 
  character(len = 4)  :: dotdmp

  integer :: endOfFile = 0
  integer :: errorChecker = 0
  integer :: errorChecker2 = 0

  logical :: lscreen
  logical :: escapeCode = .false.

  !? Open matl.dat for reading and open matl.out for writing
  call routedOpen(5, file = 'matl.dat', status = 'UNKNOWN')
  call routedOpen(4, file = 'matl.out', status = 'UNKNOWN')

  !? Variable Initializations
  na = n1 + 1
  nc = n2 + 1
  n1 = 42
  n2 = 42
  n3 = n1 - 2
  nfa = 1
  itt = 0
  ih = 2
  da0 = 1
  db0 = 1
  amm = 0.0
  dotdmp = '.DMP'

  read(5, fmt10) line

  write(*, fmt10) line
  write(4, fmt10) line

  read(5,*) matout, kontrol

  do while(.true.)
    read(5,*, iostat = endOfFile) waveFn_A, waveFn_B
    
    if(endOfFile /= 0) then
      print*, "Error Occured. Program Stopped."
      stop    
    end if

    if(waveFn_A.eq.'EXIT') then 
      stop
    end if

    !? right state reading routine. 
    if(waveFn_A(1:2).eq.'A:') then
      nfa = 3
    end if

    !? Calls the subroutine openWave to open the waveFunction Files.
    call routedOpen(1, file = waveFn_A(nfa:), status = 'OLD')

    if(waveFn_A(1:4).eq.'pow/') then
      waveFn_A = 'A:'//waveFn_A(5:)
    end if
    
    read(1, fmt11) iza, lrgla, nspna, nw, neiga, nama, title
    
    !? Checks if Title is equal to Q or Z, if it isn't, perform the loop to fix
    do while(title(1:1) .ne. 'Q' .and. title(1:1) .ne. 'Z')
      title = title(2:52)
      itt = itt + 1

      if(itt.gt.10) then 
        exit
      end if 
    end do

    if(iza.eq.1) then 
      neiga = 1
    end if

    if(iza.eq.1) then 
      ih = 1
    end if

    neig0 = neiga - 1
    nwa = nw + neig0

    !? find the right format
    fmt3a(1) = '(1X,3D26.19)'
    fmt3a(2) = '(1X,3D30.23)'
    fmt3a(3) = '(1X,3D38.31)'

    do iff = 1, 3
      if = iff

      if(if.gt.1) then 
        backspace(1)
      end if

      read(1, fmt3a(if), iostat = errorChecker2) wa(1), wa(2), amm
      print*, wa(1), wa(2), amm

      if (iostat /= 0) then 
        write(*, *) "Format match not found."
        print*, iostat
        stop
      end if 
    end do

    read(1, fmt12) nbla, (nblka(i + 1), i = 1, nbla), (nn(i + 1), i = 1, nbla - ih + 1), mar12, mar1, kono

    nblka(1) = 0
    nblxa(1) = 0

    !? Use .dml to denote the finite mass case.
    if(amm.ne.0) then 
      dotdmp = '.DML'
    end if 

    nbxa = 1
    lscreen = .false.

    do i = 1, nbla
      nblka(i + 1) = nblka(i+1) + neig0

      read(1, fmt13) nb1, nb2, la(1, nb1 + neig0), la(2, nb1 + neig0), ya(1, nb1 + neig0), ya(2, nb1 + neig0)
      
      if(nb1.eq.1.and.nb2.eq.1) then 
        lscreen = .true.
      end if

      nb1 = nb1 + neig0
      nb2 = nb2 + neig0

      do k = 1, 2
        do j=nb1,nb2
          la(k, j) = la(k, nb1)
          ya(k, j) = ya(k, nb1)
        end do
      end do

      if(i.eq.1) then 
        cycle
      end if 

      if(ya(1, nb1).eq.ya(1, nb1 - 1).and.ya(2, nb1).eq.ya(2, nb1 - 1))then
        exit
      end if

      nbxa = nbxa + 1
      nblxa(nbxa) = i - 1
    end do

    nblxa(nbxa + 1) = nbla

    read(1, fmt14) (pa(1,k), pa(2,k), sa(k), k = neiga, nwa)
    read(1,fmt3a(if)) (da(i), i = neiga, nwa)

    leiga = .false.

    !? Open .dmp file with square array of wavefunction coefficients.
    errorChecker = 0
    read(1,'(1X,A5)', iostat = errorChecker) eigen

    if(errorChecker .gt. 0)then 
      print*, "Error Occurred, Program Stopped."
      stop
    end if

    if(errorChecker .lt. 0) then
      close(1, status = 'KEEP')
    end if

    if(eigen .eq. 'EIGEN') then 
      read(1,fmt3a(if)) (eiga(i), i = 1, nwa - neiga + 1)

      leiga = .true.

      if(leiga) then
        write(*,'(2A)') 'OPENING ', waveFn_A(1:7)//dotdmp
        call routedOpen(7, file = waveFn_A(1:7)//dotdmp, form = 'UNFORMATTED', status = 'OLD')
      endif
    end if

    close(1, status = 'KEEP')

    z = iza
    
    write(4, fmt15) iza, lrgla, nspna, nw, neiga, nama, title
    write(*, fmt15) iza, lrgla, nspna, nw, neiga, nama, title
    
    call formf

    linca = la(1, neiga) + la(2, neiga) - lrgla

    !? Calculation of screened hydrogenic wavefunctions.
    if(lscreen) then
      n = lrgla + neiga
      linc = la(1, neiga) + la(2, neiga) - lrgla
      l = linc
      ya(1, 1) = ya(1, neiga)
      lincs = 0

      if(abs(ya(1, 1) - 1.0 / (linc + 1)) .gt. 1.0e-06)then 
        lincs = 1
      end if 

      ya(2, 1) = (iza - 1.0) / ((n - lincs) * iza)
      da0 = da(neiga)
      xx= 1
      
      do m = 1, neiga
        la(:, m) = la(:,neiga)
        ya(:, m) = ya(:,1)
        pa(1, m) = linc
        pa(2, m) = lrgla + m - 1
        sa(m) = 0
        da(m) = da(neiga) * xx * 2 ** (l + lrgla + 2) * (iza - 1.0) ** (lrgla + 1) * sqrt((iza - 1) * fac(n + lrgla + 1) / (iza * fac(n - lrgla) * fac(2 * l + 2))) / (fac(2 * lrgla + 2) * n * (l + 1.q0) ** (l + 2) * (n * z) ** (lrgla + 1))
        xx = -xx * 2 * (iza - 1) * (neiga - m) / ((2 * lrgla + m + 1) * m * n * iza) 
      end do
    end if

    maxca = sa(1)
    ia = pa(1, 1) + pa(2, 1) + sa(1)
    sign = 1.0

    if(da(1) .lt. 0.) then 
      sign = -1.0
    end if

    da0 = da0 * sign

    do k = 1, nwa
      da(k) = da(k) * sign
      ia1 = pa(1, k) + pa(2, k) + sa(k)

      if(ia1 .gt. ia) then 
        ia = ia1
      end if 

      if(sa(k) .gt. maxca)then 
        maxca = sa(k)
      end if
    end do

    !? left state reading routine - if left state is the same as the right state,
    !? use the file name 'SAME'.

    izb = 0

    if(waveFn_B .ne. 'SAME') then 
      nfb = 1

      if(waveFn_B(1:2).eq.'A:') then 
        nfb = 3
      end if

      !? Calls the subroutine openWave to open the waveFunction Files.
      call routedOpen(1, file = waveFn_B(nfb:), status = 'OLD')

      if(waveFn_B(1:4).eq.'pow/') then 
        waveFn_B = 'A:'//waveFn_B(5:)
      end if

      read(1, fmt11) izb, lrglb, nspnb, nw, neigb, namb, title

      do while(title(1:1) .ne. 'Q' .or. title(1:1) .ne. 'Z')
        if(title(1:1) .eq. 'Q' .or. title(1:1) .eq. 'Z')then 
          exit
        end if

        title = title(2:52)
        itt = itt + 1

        if(itt .gt. 10)then
          exit
        end if
      end do

      ih = 2

      if(izb .eq. 1) then 
        ih = 1
      end if

      neig0 = neigb - 1
      nwb = nw + neig0

      !? find the right format
      fmt3a(1) = '(1X,3D26.19)'
      fmt3a(2) = '(1X,3D30.23)'
      fmt3a(3) = '(1X,3D38.31)'

      do iff=1,3
        if = iff
        if(if.gt.1) backspace(1)
        read(1,fmt3a(if), iostat = errorChecker) wb(1), wb(2), ammb

        !! MAY NEED TO TEST THIS !!
        if(errorChecker .ne. 0 .and. iff .eq. 3)then 
          write(*,*) 'B Format match not found'
          stop
        end if

        if(errorChecker .ne. 0) then
          cycle 
        end if
      end do

      if(amm .ne. ammb) then
        write(*,'(2D20.10,A)') amm,ammb,' Masses not equal.'
        stop
      endif

      read(1, fmt12) nblb, (nblkb(i + 1), i = 1, nblb), (nn(i + 1), i = 1, nblb - ih + 1), mar12, mar1, kono

      nblkb(1) = 0
      nblxb(1) = 0
      nbxb = 1
      lscreen = .false.

      do i = 1, nblb
        nblkb(i + 1) = nblkb(i + 1) + neig0

        read(1, fmt13) nb1, nb2, lb(1, nb1 + neig0), lb(2, nb1 + neig0), yb(1, nb1 + neig0), yb(2, nb1 + neig0)

        if(nb1.eq.1.and.nb2.eq.1)then
          lscreen = .true.
        end if

        nb1 = nb1 + neig0
        nb2 = nb2 + neig0

        do k=1,2
          do j=nb1,nb2
            lb(k,j) = lb(k,nb1)
            yb(k,j) = yb(k,nb1)
          end do
        end do

        if(i.eq.1) then 
          exit
        end if

        if(yb(1, nb1) .ne. yb(1, nb1 - 1) .or. yb(2, nb1) .ne. yb(2, nb1 - 1))then 
          nbxb = nbxb + 1
          nblxb(nbxb) = i-1
        end if
      end do

      nblxb(nbxb + 1) = nblb

      read(1, fmt14)(pb(1,k),pb(2,k),sb(k),k=neigb,nwb)
      read(1,fmt3a(if)) (db(i), i = neigb, nwb)
      write(4, fmt15) izb, lrglb, nspnb, nw, neigb, namb, title
      write(*, fmt15) izb, lrglb, nspnb, nw, neigb, namb, title

            
      !?   calculation of screened hydrogenic wavefunctions.
      if(izb .ne. 0) then 
        if(lscreen) then
          linc = lb(1, neigb) + lb(2, neigb) - lrglb
          l = linc
          n = lrglb + neigb
          yb(1, 1) = yb(1, neigb)
          lincs = 0
          if(abs(yb(1, 1) -1.q0 / (linc + 1)) .gt. 1.d-06) lincs = 1
          yb(2, 1) = (izb - 1.q0) / ((n - lincs) * izb)
          db0 = db(neigb)
          xx = 1
          
          do m=1,neigb
            lb(:, m) = lb(:, neigb)
            yb(:, m) = yb(:, 1)
            pb(1,m) = linc
            pb(2,m) = lrglb + m - 1
            sb(m) = 0
            db(m) = db(neigb) * xx * 2 ** (l + lrglb + 2) * (izb - 1) ** (lrglb + 1) * sqrt((izb - 1) * fac(n + lrglb + 1) / (izb * fac(n - lrglb) * fac(2 * l + 2))) / (fac(2 * lrglb + 2) * n * (l + 1.q0) ** (l + 2) * (n * z) ** (lrglb + 1))
            xx = -xx * 2 * (izb - 1) * (neigb - m) / ((2 * lrglb + m + 1) * m * n * izb)
          end do
        end if
      end if
    end if

    lrglb = lrgla
    nspnb = nspna
    neigb = neiga
    namb = 'SAME'
    nwb = nwa
    nblb = nbla
    nblkb(1) = 0

    do i = 1, nblb
      nblkb(i + 1) = nblka(i + 1)
    end do

    nbxb = nbxa

    do i = 1, nbxb
      nblxb(i) = nblxa(i)
    end do

    nblxb(nbxb+1) = nblb

    do k = 1, nwb
      lb(:, k) = la(:, k)
      yb(:, k) = ya(:, k)
      pb(:, k) = pa(:, k)
      sb(k) = sa(k)
      db(k) = da(k)
    end do

    db0 = da0

    ib = pb(1, 1) + pb(2, 1) + sb(1)
    maxcb = sb(1)
    sign = 1.0

    if(db(1).lt.0.) then 
      sign = -1.0
    end if

    db0 = db0 * sign

    do k=1,nwb
      db(k) = db(k)*sign
      ib1 = pb(1,k) + pb(2,k) + sb(k)
      if(ib1 .gt. ib) ib = ib1
      if(sb(k) .gt. maxcb) maxcb = sb(k)
    end do

    do nbl1 = 1, nbla
      nb1 = nblka(nbl1) + 1
      nb2 = nblka(nbl1 + 1)
      iia(nbl1) = 0
      jja(nbl1) = 0
      kka(nbl1) = 0
      ijka(nbl1) = 0

      do i = nb1, nb2
        if(pa(1,i).gt.iia(nbl1)) then
          iia(nbl1) = pa(1,i)
        end if

        if(pa(2,i).gt.jja(nbl1))then 
          jja(nbl1) = pa(2,i)
        end if

        if(sa(i).gt.kka(nbl1)) then
          kka(nbl1) = sa(i)
        end if 

        j = pa(1,i) + pa(2,i) + sa(i)

        if(j.gt.ijka(nbl1))then 
          ijka(nbl1) = j
        end if
      end do
    end do

    do nbl1 = 1, nblb
      nb1 = nblkb(nbl1) + 1
      nb2 = nblkb(nbl1+1)
      iib(nbl1) = 0
      jjb(nbl1) = 0
      kkb(nbl1) = 0
      ijkb(nbl1) = 0

      do i = nb1, nb2
        if(pb(1,i).gt.iib(nbl1))then
          iib(nbl1) = pb(1,i)
        end if 

        if(pb(2,i).gt.jjb(nbl1))then 
          jjb(nbl1) = pb(2,i)
        end if

        if(sb(i).gt.kkb(nbl1)) then
          kkb(nbl1) = sb(i)
        end if 

        j = pb(1,i) + pb(2,i) + sb(i)
        
        if(j.gt.ijkb(nbl1)) then 
          ijkb(nbl1) = j
        end if
      end do
    end do

    dan = da0 / da(1)
    dbn = db0 / db(1)

    close(1,status='KEEP')

    call breit(hterm,htscr)

    hterm(:) = hterm(:) * 1048576
    htscr(:) = htscr(:) * 1048576
    htscr(:) = htscr(:) * 0.5 / (da0 * db0)

    write(*, fmt1) (i, htscr(i), hterm(i), i = 1, 3, 3)

    call routedOpen(1, file = matout, status = 'OLD')

    do while(errorChecker .ne. 0)
      if(errorChecker .gt. 0) then 
        print*, "Error Occured, program stopped"
        stop
      else if(errorChecker .lt. 0) then
        exit
      end if

      read(1,'(A14)',iostat = errorChecker) fnin
    end do
    
    backspace 1

    write(1,'(A14,1X,A14,1X,A47)')waveFn_A,waveFn_B,title
    write(1, fmt2) (hterm(i),i=1,3,3)
    write(*, fmt10) line
    write(*, fmt10)

    close(1, status = 'KEEP')

    !? RESET INPUT FILE FOR NEXT CALCULATION.
    lines(1) = line
    rewind 5
    errorChecker = 0

    do i=1,100
      ip = i-1
      read(5,'(A)', iostat = errorChecker) lines(i)
      if(errorChecker .ne. 0) then 
        if(errorChecker .gt. 0) then 
          print*, "Error Occured. Program stopped"
          stop
        else if(errorChecker .lt. 0) then 
          exit 
        end if
      end if
    end do

    rewind 5
    
    write(5,'(A)') lines(1)
    write(5,'(A40)') lines(2)

    do ii=4,ip
      write(5,'(A40)') lines(ii)
    end do

    write(5,'(A40)') lines(3)

    rewind 5

    read(5,'(/)')
  end do
end program del1


subroutine breit(hterm,htscr)

  !? Import CommonBlocks
  use a1a_block 
  use b1b_block 
  use c1c_block 
  use d1d_block 
  use f1f_block 
  use maxpow_block 

  !? Import Formatting
  use format

  !? Import OpenWave
  use wavExt

  implicit real*16 (a-h,o-z)
  
  dimension :: hterm(10),htscr(10)
  
  integer :: p,q,s,a1,b1,c1
  
  logical :: lgo,lcalc
  logical :: escapeCode = .false.


  !? FOR D STATES
  data isumd,mp2d,mq2d,ms2d,msld/-2,1,1,2,-10/

  !? FOR P STATES
  !? DATA ISUMD,MP2D,MQ2D,MS2D,MSLD/-2,1,1,2,-10/

  !? PARONIC VERSION.  REMOVE NEXT TWO EXECUTABLE STATEMENTS AND REVERSE THE
  !? SIGN OF SGN TO OBTAIN THE NORMAL VERSION..

  !? IST CHOOSES THE COORDINATES OF ELECTRON 1 OR 2 AND 'I' SPECIFIES WHICH
  !? SET OF COEFFICIENTS OF PL(COSTHETA12) IS TO BE USED.

  !? READ(*,*) ISUMD,MP2D,MQ2D,MS2D,MSLD
  !? NSPNA = 1 - NSPNA
  !? NSPNB = 1 - NSPNB

  ermac = 1.1_16

  do while(test.gt.1.1_16)
      ermac = ermac/2
      test = 1.1_16 + ermac
  end do
  
  kross = 9
  klog = -2
  kang = 0
  sq3 = sqrt(3.0_16)
  mla = lrgla
  mlb = lrglb
  md = - mla + mlb

  do i=1,10
      hterm(i) = 0.0
      htscr(i) = 0.0
  end do

  ss5 = (-1) ** (lrglb - mlb) * sqrt(15.0_16 / 8.0_16)
  fj = f3j(lrglb, -mlb, 2, md, lrgla, mla)

  if(fj.eq.0.0) then 
    write(4,fmt8)
    return
  end if

  ss5 = ss5 / fj

  !? CALCULATES THE REDUCED MATRIX ELEMENT (L'||L||L)(S'||S||S).  MULTIPLY RESULT
  !? BY F6J(J,1,S',L,L',S)*(-1)**(L+S'+J)*ALPHA**2.

  fj = f3j(lrglb, -mlb, 1, md, lrgla, mla) * (-1.) ** (lrglb - mlb)

  if(fj.eq.0.0) then 
    write(4,fmt8)
    return
  end if

  if((nspna + nspnb - 1) .lt. 0) then
    write(4,fmt7)
    return 
  else if((nspna + nspnb - 1) .eq. 0)then 
    so = -z * sq3 / fj / 2.0
    if(nspna.eq.0)then
      so = - so
    end if
    soo = so / z
  else if((nspna + nspnb - 1) .gt. 0) then 
    so = z * sq3 / fj / sqrt(2.0_16)
    soo = -3.0 * so / z
  end if

  !?   THE OUTER LOOP CORRESPONDS TO THE RIGHT HAND STATE AND IS DENOTED BY
  !?   'A' OR'1' OR '2' IN VARIOUS PLACES.

  !?   SUMMATION OVER RIGHT HAND STATE DIRECT AND EXCHANGE TERMS.

  do idx = 1, 2
    !? SUMMATION OVER BLOCKS OF TERMS WITH THE SAME NON-LINEAR PAREMETERS.
    do nhi = 1, nbxa
      nbx1 = nblxa(nhi) + 1
      nbx2 = nblxa(nhi+1)
      do nhj = 1, nbxb
        nbx3 = nblxb(nhj) + 1
        nbx4 = nblxb(nhj+1)

        !?   THE OUTER LOOP CORRESPONDS TO THE RIGHT HAND STATE AND IS DENOTED BY
        !?  'A' OR'1' OR '2' IN VARIOUS PLACES.

        do nbl1 = nbx1, nbx2
          nb1 = nblka(nbl1) + 1
          nb2 = nblka(nbl1 + 1)

          !?   THE INNER LOOP CORRESPONDS TO THE LEFT HAND STATE AND IS DENOTED BY
          !?   'B' OR '3' OR '4' IN VARIOUS PLACES (3 MEANS ELECTRON 1 AND 4 MEANS ELECTRON
          !?    2).

          do nbl3 = nbx3, nbx4
            nb3 = nblkb(nbl3) + 1
            nb4 = nblkb(nbl3 + 1)
            lgo = .false.

            if(nbl1.eq.1.and.nbl3.eq.1) then 
              lgo = .true.
            end if 

            if(izb.eq.0.and.nbl3.lt.nbl1)then 
              escapeCode = .true.
              exit
            end if

            !?   SUMMATION OVER LEFT HAND STATE DIRECT AND EXCHANGE TERMS.
            !?   THE ORDERING IS D-D, E-E, D-E AND E-D.

            do jdx = 1, 2
              lcalc = .false.

              if(nbl1.le.ncut.and.nbl3.le.ncut.and.(idx.eq.1.or.lrgla.eq.0))then 
                lcalc = .true.
              end if

              stmax = 0.
              ide =  jdx
              jde =  jdx
              
              if(idx.eq.2) then 
                ide = 3 - jde
              end if

              id = ide - 1
              jd = jde - 1
              idep = 3 - ide
              jdep = 3 - jde
              sgn1 = 1.

              if(nspna.eq.1) then 
                sgn1 = (-1.)**(ide-1)
              end if

              sgn = sgn1

              if(nspnb.eq.1) then 
                sgn = sgn1*(-1.)**(jde-1)
              end if

              !? FOR PARONIC STATES, USE THE FOLLOWING -
              !? SGN1 = 1.
              !? IF(NSPNA.EQ.0) SGN1 = (-1.)**(IDE-1)
              !? SGN = SGN1
              !? IF(NSPNB.EQ.0) SGN = SGN1*(-1.)**(JDE-1)

              l1 = la(ide,nb1)
              l2 = la(idep,nb1)
              l3 = lb(jde,nb3)
              l4 = lb(jdep,nb3)
              ll1 = l1*(l1+1)
              ll3 = l3*(l3+1)
              ll4 = l4*(l4+1)

              call cross(l3,l4,l1,l2,lrglb,mlb,lrgla,mla,ide,jde)

              mlp1 = mlt(1)
              lp1 = lt(1)

              do l = 1, lp1 + 1
                cplt(l,8) = cplt(l,11)
              end do

              do ll = mlp1, lp1, 2
                l = lp1 - ll + mlp1

                !? REDUCE cplt(L,8) = cplt(L,11) BY RECURSION RELATION.
                if(l.eq.1.and.abs(cplt(l+1,8)).gt.1.d-10)then 
                  stop 2
                end if

                cplt(l,8) = (2 * l - 1) * cplt(l + 1, 8)

                if(l.eq.1) then 
                  exit 
                end if

                cplt(l - 1, 8) = cplt(l - 1, 8) + cplt(l + 1, 8)

                if(abs(cplt(l - 1, 8)).lt.1.d-10)then 
                  cplt(l - 1, 8) = 0.
                end if
              end do

              mlt(8) = 0
              lt(8) = 0

              if(lt(11).ne.0)then 
                mlt(8) = mlt(11) + 1
                lt(8) = lt(11) - 1
              end if

              l8 = lt(8)

              if(kang.ne.0) then 
                write(4,fmt3) l3, l4, l1, l2, lrglb, mlb, lrgla, mla, ide, jde, nb3, nb4, nb1, nb2, izb
                write(4,fmt4)

                do ic = 1, 16
                  if(lt(ic).le.0) then
                    exit
                  end if
                  ic2 = lt(ic)
                  write(4,fmt4) ic, (cplt(ic3, ic), ic3 = 1, ic2)
                end do
              end if

              if(jdx.ne.2)then
                if(nbl1.le.nbx1.and.nbl3.le.nbx3) then 
                  y11 = ya(ide,nb1) + yb(jde,nb3)
                  y22 = ya(idep,nb1) + yb(jdep,nb3)
                  kex = 0

                  if(y22.gt.y11) then
                    y22 = ya(ide,nb1) + yb(jde,nb3)
                    y11 = ya(idep,nb1) + yb(jdep,nb3)
                    kex = 1
                  endif

                  lmax = lt(2)
                  mp1 = 4
                  ms1 = 2
                  msl = lt(1) + msld
                  imin = -3 + lrgla + lrglb
                  lsuma = 2*id
                  lsumb = l2*(iabs(l1-l2)/2)
                  ltota = lrgla + neiga - 1
                  ltotb = lrglb + neigb - 1

                  if((nb1.ne.1.and.nb3.ne.1).or.z.le.1) then 
                    isum = ijkb(nbl3) + ijka(nbl1) + isumd
                    mp2 = iib(nbl3) + (1-id)*iia(nbl1) + id*jja(nbl1) + mp1 + mp2d + 1 - id
                    mq2 = jjb(nbl3) + (1-id)*jja(nbl1) + id*iia(nbl1)+ mp1+mq2d + id
                    ms2 = ms1 + kkb(nbl3) + kka(nbl1) + 2*lt(1) + ms2d
                    escapeCode = .true.
                  end if


                  if((nb3.ne.1.or.nb1.ne.1) .and. (escapeCode .eqv. .false.)) then 
                    if(nb3.ne.1) then 

                      !? LIMITS FOR HYDROGENIC '1,N' TERMS.
                      isum = ltota + ijkb(nbl3) + isumd
                      mp2 = iib(nbl3) + id*ltota + mp1 + mp2d
                      mq2 = jjb(nbl3) + (1-id)*ltota + mp1 + mq2d
                      ms2 = ms1 + kkb(nbl3) + 2*lt(1) + ms2d
                      escapeCode = .true.
                    end if

                    if(escapeCode .eqv. .false.)then
                      !? LIMITS FOR HYDROGENIC 'N,1' TERMS.
                      isum = ltotb + ijka(nbl1) + isumd
                      mp2 = (1-id)*iia(nbl1) + id*(jja(nbl1)+ltotb) + mp1 + mp2d + 1
                      mq2 = (1-id)*(jja(nbl1)+ltotb) + id*iia(nbl1) + mp1 + mq2d + id

                      if(izb.ne.0) then 
                        mq2 = mq2 + id*(ltotb-1)
                      end if

                      ms2 = ms1 + maxca + 2*lt(1) + ms2d
                      escapeCode = .true.
                    end if
                  end if

                  if(escapeCode .eqv. .false.)then
                    isum = ltotb + ltota + isumd
                    mp2 = (1-id)*ltota + id*(ltotb+ltota) + mp1 + mp2d
                    mq2 = (1-id)*(ltotb+ltota) + id*ltota + mp1 + mq2d
                    ms2 = ms1 + 2*lt(1) + ms2d + id
                  end if

                  if(id.eq.0.and.(nbx4.gt.nbx3.or.nbx2.gt.nbx1)) then 
                    ms2 = ms2 +2
                  end if

                  if(id.eq.0.and.(nbx4.gt.nbx3.and.nbx2.gt.nbx1)) then 
                    ms2 = ms2 +2
                  end if

                  if(id.eq.1.and.(nbx2.gt.nbx1.or.nbx4.gt.nbx3).and.(l2.gt.0.or.l3.gt.0))then 
                    ms2 = ms2 +2
                  end if

                  if((nbx4.gt.nbx3.or.nbx2.gt.nbx1))then 
                    lmax =min0(lrgla +2,lmax+1)
                  end if

                  if((nbx4.gt.nbx3.and.nbx2.gt.nbx1))then 
                    lmax =min0(lrgla +2,lmax+1)
                  end if

                  if(msl.gt.12)then 
                    stop
                  end if

                  if(kex.eq.1) then
                    mtemp = mp2
                    mp2 = mq2
                    mq2 = mtemp
                  end if

                  if(nbl3.eq.nbl1) then
                    write(*,fmt5) kex,mp1,mp2,mq2,ms1,ms2,isum,y11,y22,l3,l4,l1,l2,nbl1,nbl3
                  end if

                  call genint(id,nbl3,nbl1,lmax,fac,0.q0)
                end if
              end if

              delt12 = 0.0
              jdxx = jdx

              if(kex.eq.1) then 
                jdxx = 3 - jdx
              end if

              m1 = mlt(1)
              m2 = lt(1)

              do j=m1,m2
                delt12 = delt12 + cplt(j,1)
              end do

              !? SUMMATION OVER RIGHT HAND STATE BASIS FUNCTIONS.
              do i=nb1,nb2
                a1 = pa(ide,i)
                b1 = pa(idep,i)
                c1 = sa(i)
                y1 = ya(ide,i)
                y2 = ya(idep,i)
                l9 = max(lt(11),lt(12))
                ml9 = min(mlt(11),mlt(12))

                if(ml9.le.0)then
                  ml9 = max(mlt(11),mlt(12))
                end if

                if(ml9.gt.0)then
                  do l=ml9,l9,2
                    cplt(l,9) = a1*cplt(l,11) - cplt(l,12)
                  end do

                  if(abs(cplt(l9,9)).lt.1.d-08)then 
                    l9 = l9 - 2
                  end if
                end if

                mlt(9) = ml9
                lt(9) = l9

                !? SUMMATION OVER LEFT HAND STATE BASIS FUNCTIONS.
                do j=nb3,nb4
                  if(izb.eq.0.and.j.lt.i)then
                    exit
                  end if

                  psi2 = da(i)*db(j)*sgn

                  if(izb.eq.0.and.j.ne.i)then
                    psi2 = psi2*2.
                  end if

                  y3 = yb(jde,j)
                  y4 = yb(jdep,j)
                  p = pa(ide,i) + pb(jde,j)
                  q = pa(idep,i) + pb(jdep,j)
                  s = sa(i) + sb(j)
                  y = y1 + y3
                  yp = y2 + y4
                  c14 = sa(i)*sb(j)
                
                  !? CALCULATION OF SPIN-DEPENDENT STONE CORRECTION.
                  sum = 0.

                  if(l9.gt.0)then 
                    sum = spl(p-1,q-2,s,9)
                  end if

                  if(l8.gt.0) then 
                    sum = sum - y1*spl(p-1,q-3,s+2,8)/(s+2)

                    if(c1.eq.0) then 
                      sum = sum + c1*spl(p,q-3,s,8)/s
                    end if
                  end if

                  hterm(1) = hterm(1) - 2. * sum * so * psi2

                  if(lgo) then 
                    htscr(1) = htscr(1) - 2.*sum*so*psi2
                  end if
                end do
              end do
            end do
          end do
          if(escapeCode) then 
              exit
          end if
        end do
        if(escapeCode) then 
          exit
        end if
      end do
      if(escapeCode) then 
          exit
      end if
    end do
    if(escapeCode) then 
      exit
    end if
  end do

  return

end subroutine breit