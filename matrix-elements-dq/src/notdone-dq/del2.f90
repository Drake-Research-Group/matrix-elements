!* PROGRAM FOR STONE DEL2 SPIN-INDEPENDENT TERM *!

include 'cross1.f90'
include 'spin68.f90'

program del2

      !? Import CommonBlocks
      use a1a_block
      use h1h_block
      use f1f_block
      use b1b_block
      use maxpow_block

      !? Import Formatting
      use format

      !? Import File Management
      use wavExt

      implicit real*16 (a-h,o-z)
  
      character :: cx*1
      character :: lines(100)

      integer :: endOfFile = 0
      integer :: errorChecker = 0

      !? More Initializations
      dimension :: hterm(10),htscr(10),nn(50)
      character(len = 16) :: matout
      character title*52,nama*12,namb*12,line*51,fmt3a(3)*12, eigen*5, dotdmp*4

      logical :: lscreen
      logical :: escapeCode = .false.

      !? Open matl.dat for reading and open matl.out for writing
      call routedOpen(5, file = 'matl.dat', status = 'UNKNOWN')
      call routedOpen(4, file = 'matl.out' ,status = 'UNKNOWN')

      !? Variable Initializations
      n1 = 42
      n2 = 42
      na = n1 + 1
      nc = n2 + 1
      n3 = n1 - 2
      nfa = 1
      itt = 0
      ih = 2
      da0 = 1
      db0 = 1
      amm = 0.
      dotdmp = '.DMP'

      read(5, fmt10) line

      write(*,fmt10) line
      write(4,fmt10) line

      read(5,*) matout,kontrol

      do while(endOfFile .eq. 0)
            read(5,*, iostat = endOfFile) waveFn_A,waveFn_B
            if(endOfFile .gt. 0) then
                  print*, "Error Occured. Program Stopped."
                  stop

            else if(endOfFile .lt. 0) then 
                  stop         
            end if


            if(waveFn_A.eq.'EXIT') then 
                  stop
            end if

            !? right state reading routine. 
            if(waveFn_A(1:2).eq.'A:') then
                  nfa = 3
            end if

            !? Calls the subroutine routedOpen to open the waveFunction Files.
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

            errorChecker = 0

            do iff = 1, 3
                  if = iff

                  if(if.gt.1) then 
                        backspace(1)
                  end if

                  read(1, fmt3a(if), iostat = errorChecker), wa(1), wa(2), amm

                  if(errorChecker .lt. 0) then 
                        write(*,*) 'A Format match not found'
                        stop

                  else if(errorChecker .eq. 0) then
                        exit
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
                  end if
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

                  if(abs(ya(1, 1) - 1.q0 / (linc + 1)).gt.1.d-06)then 
                        lincs = 1
                  end if 

                  ya(2, 1) = (iza - 1.q0) / ((n - lincs) * iza)
                  da0 = da(neiga)
                  xx= 1
                  
                  do m = 1, neiga
                        do k = 1, 2
                        la(k,m) = la(k,neiga)
                        ya(k,m) = ya(k,1)
                        end do
                        pa(1,m) = linc
                        pa(2,m) = lrgla + m - 1
                        sa(m) = 0
                        da(m) = da(neiga) * xx * 2 ** (l + lrgla + 2) * (iza - 1.q0) ** (lrgla + 1) * sqrt((iza - 1) * fac(n + lrgla + 1) / (iza * fac(n - lrgla) * fac(2 * l + 2))) / (fac(2 * lrgla + 2) * n * (l + 1.q0) ** (l + 2) * (n * z) ** (lrgla + 1))
                        xx = -xx * 2 * (iza - 1) * (neiga - m) / ((2 * lrgla + m + 1) * m * n * iza) 
                  end do
            end if

            maxca = sa(1)
            ia = pa(1, 1) + pa(2, 1) + sa(1)
            sign = 1.

            if(da(1) .lt. 0.) then 
                  sign = -1.
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

            !?   left state reading routine - if left state is the same as the right state,
            !?   use the file name 'SAME'.

            izb = 0

            if(waveFn_B .ne. 'SAME') then 
            nfb = 1

            if(waveFn_B(1:2).eq.'A:') then 
                  nfb = 3
            end if

            !? Calls the subroutine routedOpen to open the waveFunction Files.
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

            !?   find the right format
            fmt3a(1) = '(1X,3D26.19)'
            fmt3a(2) = '(1X,3D30.23)'
            fmt3a(3) = '(1X,3D38.31)'

            do iff=1,3
                  if = iff
                  if(if.gt.1) then 
                        backspace(1)
                  end if 

                  read(1,fmt3a(if), iostat = errorChecker) wb(1), wb(2), ammb

                  if(errorChecker .ne. 0) then
                        cycle 
                  end if

                  escapeCode = .true.
                  exit

            end do

            if(.not. escapeCode) then
                  write(*,*) 'B Format match not found'
                  stop
            end if 

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
                                    do k=1,2
                                          lb(k, m) = lb(k, neigb)
                                          yb(k, m) = yb(k, 1)
                                    end do

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
                  do i = 1, 2
                        lb(i, k) = la(i, k)
                        yb(i, k) = ya(i, k)
                        pb(i, k) = pa(i, k)
                  end do
                  sb(k) = sa(k)
                  db(k) = da(k)
            end do

            db0 = da0

            ib = pb(1, 1) + pb(2, 1) + sb(1)
            maxcb = sb(1)
            sign = 1.

            if(db(1).lt.0.) then 
                  sign = -1.
            end if

            db0 = db0 * sign

            do k = 1, nwb
                  db(k) = db(k)*sign
                  ib1 = pb(1,k) + pb(2,k) + sb(k)

                  if(ib1 .gt. ib) then 
                        ib = ib1
                  end if

                  if(sb(k) .gt. maxcb)then 
                        maxcb = sb(k)
                  end if
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
            leigb = .false.
            errorChecker = 0

            read(1,'(1X,A5)',iostat = errorChecker) eigen

            if(eigen.eq.'EIGEN') then
                  read(1, fmt16) (eigb(i), i = 1, nwb - neigb + 1)
                  leigb = .true.
            end if
      
            close(1,status='KEEP')

            !? OPEN .DMP FILE WITH SQUARE ARRAY OF WAVEFUNCTION COEFFICIENTS.
            if(leigb) then
                  write(*,'(2A)') 'DEL2 OPENING ',waveFn_B(1:7)//dotdmp
                  call routedOpen(8, file = waveFn_B(1:7) //dotdmp, form = 'UNFORMATTED', status = 'OLD')
            endif

            close(1,status='KEEP')

            call breit(hterm,htscr,amm)

            do i=1,3
                  hterm(i) = hterm(i)*1048576
                  htscr(i) = htscr(i)*1048576
                  htscr(i) = htscr(i)*0.5/(da0*db0)
            end do

            hterm(1) = hterm(1) + 1.
            htscr(1) = (2*lrgla+1)*(htscr(1) + 1.)
            hterm(2) = hterm(2) + 1.

            write(*, fmt1) (i,htscr(i),hterm(i),i=1,3,3)
            call routedOpen(1, file = matout, status = 'OLD')

            errorChecker = 0
            
            do while(errorChecker .eq. 0)
                  read(1,'(A14)', iostat = errorChecker) fnin
                  if(errorChecker .gt. 0)then 
                        print*, "An error has occurred. Stopping the program."
                        stop 
                  
                  else if(errorChecker .lt. 0)then 
                        exit 
                  end if
            end do

            backspace 1

            write(1,'(A14,1X,A14,1X,A51)')waveFn_A,waveFn_B,title
            write(1, fmt2) (hterm(i),i=1,3,3)
            write(*, fmt10) line
            write(*, fmt10)

            close(1, status = 'KEEP')
      end do

      rewind 5
      errorChecker = 0

      do i=1,100
            ip = i-1
            read(5,'(A)', iostat = errorChecker) lines(i)
            if(errorChecker .gt. 0) then 
                  print*, "An Error Occurred. Stopping the program."
            
            else if(errorChecker .lt. 0)then
                  exit
            end if
      end do

      rewind 5

      write(5,'(A)') lines(1)
      write(5,'(A16)') lines(2)

      do ii=4,ip
            write(5,'(A32)') lines(ii)
      end do

      write(5,'(A32)') lines(3)

      rewind 5

      read(5,'(/)')
end program del2

subroutine breit(hterm,htscr,amm)

      !?Import CommonBlock
      use a1a_block
      use b1b_block
      use c1c_block
      use d1d_block
      use f1f_block
      use h1h_block
      use maxpow_block

      !?Import Formatting
      use format

      implicit real*16 (a-h,o-z)

      dimension hterm(10),htscr(10),ipow(30),cpltsto(12),trb(3601),pert(3601),tr(3601,3601),temp(3601,3601)

      integer :: p
      integer :: q
      integer :: s
      integer :: a1,b1,c1
      integer :: a4,b4,c4
      integer :: pam
      integer :: qam

      logical :: lgo, lcalc
      logical :: escapeCode = .false.

      character op*3,title*15,dotdmp*4


      !? FOR D STATES
      data isumd,mp2d,mq2d,ms2d,msld/-1,2,2,1,-10/

      !? FOR P STATES
      !? DATA ISUMD,MP2D,MQ2D,MS2D,MSLD/-1,1,2,1,-10/

      data ipow/-2,-1,0, -2,1,-2, 0,-1,-2, -2,0,0, -1,-1,0, 0,0,-2,-1,1,-2, -1,0,0, -1,1,-2, -1,0,-2/

      !? IST CHOOSES THE COORDINATES OF ELECTRON 1 OR 2 AND 'I' SPECIFIES WHICH
      !? SET OF COEFFICIENTS OF PL(COSTHETA12) IS TO BE USED.

      ermac = 1.0d0

      do while(test .gt. 1.0d0)
            ermac = ermac/2
            test = 1.0d0 + ermac
      end do 

      !? USE .DML TO DENOTE THE FINITE MASS CASE.
      dotdmp = '.DMP'

      if(amm.ne.0) then 
            dotdmp = '.DML'
      end if 

      kross = 9
      klog = -2
      kang = 1
      kang = 0
      sq3 = dsqrt(3.d0)
      mla = lrgla
      mlb = lrglb
      md = - mla + mlb

      do i=1,10
            hterm(i) = 0.0
            htscr(i) = 0.0
      end do

      tr = 0
      trb = 0

      !? THE OUTER LOOP CORRESPONDS TO THE RIGHT HAND STATE AND IS DENOTED BY
      !? 'A' OR'1' OR '2' IN VARIOUS PLACES.
      !? SUMMATION OVER RIGHT HAND STATE DIRECT AND EXCHANGE TERMS.

      do idx=1,2
            !? SUMMATION OVER BLOCKS OF TERMS WITH THE SAME NON-LINEAR PAREMETERS.
            do nhi=1,nbxa
                  nbx1 = nblxa(nhi) + 1
                  nbx2 = nblxa(nhi+1)
                  do nhj=1,nbxb
                        nbx3 = nblxb(nhj) + 1
                        nbx4 = nblxb(nhj+1)

                        !? THE OUTER LOOP CORRESPONDS TO THE RIGHT HAND STATE AND IS DENOTED BY
                        !? 'A' OR'1' OR '2' IN VARIOUS PLACES.

                        do nbl1=nbx1,nbx2
                              nb1 = nblka(nbl1) + 1
                              nb2 = nblka(nbl1+1)

                              !? THE INNER LOOP CORRESPONDS TO THE LEFT HAND STATE AND IS DENOTED BY
                              !?'B' OR '3' OR '4' IN VARIOUS PLACES (3 MEANS ELECTRON 1 AND 4 MEANS ELECTRON
                              !? 2).

                              do nbl3=nbx3,nbx4
                                    nb3 = nblkb(nbl3) + 1
                                    nb4 = nblkb(nbl3 + 1)
                                    lgo = .false.

                                    if(nbl1.eq.1.and.nbl3.eq.1) then 
                                          lgo = .true.
                                    end if 

                                    if(izb.eq.0.and.nbl3.lt.nbl1) then 
                                          cycle
                                    end if 

                                    !? SUMMATION OVER LEFT HAND STATE DIRECT AND EXCHANGE TERMS.
                                    !? THE ORDERING IS D-D, E-E, D-E AND E-D.

                                    do jdx=1,2
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

                                          l1 = la(ide,nb1)
                                          l2 = la(idep,nb1)
                                          l3 = lb(jde,nb3)
                                          l4 = lb(jdep,nb3)
                                          ll1 = l1*(l1+1)
                                          ll3 = l3*(l3+1)
                                          ll4 = l4*(l4+1)

                                          call cross(l3,l4,l1,l2,lrglb,mlb,lrgla,mla,ide,jde)

                                          mlp2 = mlt(2)
                                          mlp3 = mlt(3)
                                          lp2 = lt(2)
                                          lp3 = lt(3)

                                          do ll=mlp2,lp2,2
                                                cplt(ll,15) = 0.
                                                cplt(ll,16) = 0.
                                                cplt(ll,4) = cplt(ll,9)
                                                cplt(ll,5) = cplt(ll,10)
                                          end do

                                          do ll=mlp3,lp3,2
                                                cpltsto(ll) = cplt(ll,7)
                                                cplt(ll,7) = 0.
                                          end do

                                          do n=4,14
                                                mlt(n) = mlp2
                                                lt(n) = lp2
                                          end do

                                          mlt(15) = mlp3
                                          lt(15) = lp3
                                          mlt(16) = mlp3
                                          lt(16) = lp3

                                          if(kang.ne.0) then
                                                write(4, fmt3) l3,l4,l1,l2,lrglb,mlb,lrgla,mla,ide,jde,nb3,nb4,nb1,nb2,izb
                                                write(4, fmt4)

                                                do ic=1,16
                                                      if(lt(ic).le.0) then 
                                                            cycle 
                                                      end if 
                                                      ic2 = lt(ic)
                                                      write(4,fmt4) ic,(cplt(ic3,ic), ic3=  1,ic2)
                                                end do
                                          end if 

                                          if((nbl1.le.nbx1.and.nbl3.le.nbx3) .and. (jdx.le.2)) then 

                                                y11 = ya(ide,nb1) + yb(jde,nb3)
                                                y22 = ya(idep,nb1) + yb(jdep,nb3)
                                                kex = 0

                                                if(y22.gt.y11) then
                                                      y22 = ya(ide,nb1) + yb(jde,nb3)
                                                      y11 = ya(idep,nb1) + yb(jdep,nb3)
                                                      kex = 1
                                                endif

                                                lmax = lt(3)
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
                                                      mp2 = iib(nbl3) + (1-id)*iia(nbl1) + id*jja(nbl1) + mp1 + mp2d+ 2 - id
                                                      mq2 = jjb(nbl3) + (1-id)*jja(nbl1) + id*iia(nbl1)+ mp1+mq2d- 1 + 2*id
                                                      ms2 = ms1 + kkb(nbl3) + kka(nbl1) + 2*lt(1) + ms2d

                                                      escapeCode = .true.
                                                end if 

                                                if((nb3.ne.1.or.nb1.ne.1) .and. (.not. escapeCode)) then 
                                                      
                                                      !? LIMITS FOR HYDROGENIC '1,N' TERMS.
                                                      if(nb3.ne.1) then 
                                                            isum = ltota + ijkb(nbl3) + isumd
                                                            mp2 = iib(nbl3) + id*ltota + mp1 + mp2d + 1
                                                            mq2 = jjb(nbl3) + (1-id)*ltota + mp1 + mq2d - 1 + id
                                                            ms2 = ms1 + kkb(nbl3) + 2*lt(1) + ms2d
                                                            
                                                            escapeCode = .true.
                                                      end if 
                                                      

                                                      !? LIMITS FOR HYDROGENIC 'N,1' TERMS.
                                                      if(.not. escapeCode) then
                                                            isum = ltotb + ijka(nbl1) + isumd
                                                            mp2 = (1-id)*iia(nbl1) + id*(jja(nbl1)+ltotb) + mp1 + mp2d + 1
                                                            mq2 = (1-id)*(jja(nbl1)+ltotb) + id*iia(nbl1) + mp1 + mq2d - 1 +id+1
                                                            ms2 = ms1 + kka(nbl1) + 2*lt(1) + ms2d + 1
                                                            
                                                            escapeCode = .true.
                                                      end if 
                                    
                                                end if 

                                                if(.not. escapeCode) then
                                                      isum = ltotb + ltota + isumd
                                                      mp2 = (1-id)*ltota + id*(ltotb+ltota) + mp1 + mp2d + 1 - id
                                                      mq2 = (1-id)*(ltotb+ltota) + id*ltota + mp1 + mq2d - 1 + id
                                                      ms2 = ms1 + 2*lt(1) + ms2d
                                                end if 

                                                if(id.eq.0.and.(nbx4.gt.nbx3.or.nbx2.gt.nbx1)) then 
                                                      ms2 = ms2 +lrgla
                                                end if 

                                                if(id.eq.0.and.(nbx4.gt.nbx3.and.nbx2.gt.nbx1)) then 
                                                      ms2 = ms2 +lrgla
                                                end if 

                                                if(id.eq.1.and.nbx2.gt.nbx1.and.l3.gt.0) then 
                                                      ms2 = ms2 +lrgla
                                                end if 

                                                if((nbx4.gt.nbx3.or.nbx2.gt.nbx1)) then 
                                                      lmax =min0(lrgla +3,lmax+1)
                                                end if 

                                                if((nbx4.gt.nbx3.and.nbx2.gt.nbx1)) then 
                                                      lmax =min0(lrgla +3,lmax+1)
                                                end if 

                                                if(msl.gt.12) then 
                                                      stop 180
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
                                                ar1 = 2.*a1*a1 - l1*(l1+1)
                                                ar2 = -2.*y1*(2.*a1 + 1.)
                                                ar3 = -2.*c1*y1
                                                ar4 = 2.*y1*y1

                                                do ll=mlp3,lp3,2
                                                      cplt(ll,15) = c1*y2*(cplt(ll,1) + cplt(ll,3))
                                                      cplt(ll,16) = c1*(-b1*cplt(ll,3) -cpltsto(ll) + (2.*a1-b1)*cplt(ll,1))
                                                end do

                                                do ll=mlp2,lp2,2
                                                      cplt(ll,7) = 2.*a1*b1*cplt(ll,2) + b1*cplt(ll,5) + 2.*a1*cplt(ll,4) + cplt(ll,6)
                                                      cplt(ll,8) = -c1*(2.*a1*cplt(ll,2) + cplt(ll,5))
                                                      cplt(ll,9) = 2.*c1*(b1*cplt(ll,2) + cplt(ll,4))
                                                      cplt(ll,10) = -y2*(2.*a1*cplt(ll,2) + cplt(ll,5))
                                                      cplt(ll,11) = -2.*y1*(b1*cplt(ll,2) + cplt(ll,4))
                                                      cplt(ll,12) = -2.*c1*y2*cplt(ll,2)
                                                      cplt(ll,13) = 2.*c1*y1*cplt(ll,2)
                                                      cplt(ll,14) = 2.*y1*y2*cplt(ll,2)
                                                end do

                                                do n=7,14
                                                      lt(n) = lp2
                                                      do ll=mlp2,lp2,2
                                                            lpn = lt(n)
                                                            if(abs(cplt(lpn,n)).lt.1.d-06) then 
                                                                  lt(n) = lt(n) - 2
                                                            end if 
                                                      end do
                                                end do

                                                do n=15,16
                                                      lt(n) = lp3
                                                      do ll=mlp3,lp3,2
                                                            lpn = lt(n)
                                                            if(abs(cplt(lpn,n)).lt.1.d-06) lt(n) = lt(n) - 2
                                                      end do
                                                end do

                                                !? SUMMATION OVER LEFT HAND STATE BASIS FUNCTIONS.
                                                do j=nb3,nb4
                                                      if(izb.eq.0.and.j.lt.i) then 
                                                            cycle
                                                      end if 

                                                      psi2 = da(i)*db(j)*sgn

                                                      if(izb.eq.0.and.j.ne.i) then 
                                                            psi2 = psi2*2.
                                                      end if 

                                                      y3 = yb(jde,j)
                                                      y4 = yb(jdep,j)
                                                      p = pa(ide,i) + pb(jde,j)
                                                      q = pa(idep,i) + pb(jdep,j)
                                                      s = sa(i) + sb(j)

                                                      !? CALCULATION OF STONE DELTA(2) CORRECTION.
                                                      kkr = -2
                                                      sum = 0.

                                                      do itr=7,16
                                                            kkr = kkr + 3
                                                            if(lt(itr).gt.0) then 
                                                                  ip = p + ipow(kkr)
                                                                  iq = q + ipow(kkr+1)
                                                                  is = s + ipow(kkr+2)
                                                                  sum = sum + spl(ip,iq,is,itr)
                                                            end if
                                                      end do

                                                      if(ar1.ne.0.) then 
                                                            sum = sum + ar1*spl(p-3,q,s,1)
                                                      end if 

                                                      if(ar3.ne.0.) then 
                                                            sum = sum + ar3*spl(p,q,s-2,1)
                                                      end if 

                                                      sum = sum + ar2*spl(p-2,q,s,1)+ ar4*spl(p-1,q,s,1)
                                                      trb(j) = trb(j) + sum*psi2/db(j)
                                                      tr(j,i) = tr(j,i) + sum*psi2/(da(i)*db(j))
                                                      hterm(1) = hterm(1) + sum*psi2

                                                      if(lgo) then 
                                                            htscr(1) = htscr(1) + sum*psi2
                                                      end if 
                                                end do
                                          end do
                                    end do 
                              end do
                        end do
                  end do
            end do
      end do

      if(leigb) then
            op = 'ST2'
            write(*,'(2A)') 'ST2 OPENING ',waveFn_A(1:4)//waveFn_B(1:4)//op//dotdmp
            open(9,file=waveFn_A(1:4)//waveFn_B(1:4)//op//dotdmp,form='UNFORMATTED',status='UNKNOWN')
            nar = 1

            !? ADJUST NEIGB AND NEIGA TO ALLOW FOR NEIG = n-1 FOR TRIPLET-S STATES

            neffb = neigb
            neffa = neiga

            if(lrglb.eq.0.and.nspnb.eq.1) then 
                  neffb = neigb-1
            end if 

            if(lrgla.eq.0.and.nspna.eq.1) then 
                  neffa = neiga-1
            end if 

            if(leiga) then 
                  nar = nwa - neiga + 1
            end if 

            nbr = nwb - neigb + 1
            eiga(1) = wa(2)

            do i=1,nbr
                  call input_b
                  sumss = 0
                  do j=1,nwb
                        sumss = sumss + trb(j)*db(j)
                  end do
                  pert(i) = 1048576*sumss*z**3
            end do

            write(*,'(A,1PD24.16)') 'Check Ratio =',pert(neffb)/(hterm(1)*1048576*z**3)
            rewind(8)

            do i=1,nbr
                  call input_b

                  do j=1,nwa
                        sum = 0
                        do k=1,nwb
                              sum = sum + db(k)*tr(k,j)
                        end do

                        temp(i,j) = sum
                  end do 
            end do

            rewind(7)
            rewind(8)
            rewind(9)

            do j=1,nar
                  if(nar.gt.1) then
                        call input_a
                  end if
                  do i=1,nbr
                        sum = 0

                        do k=1,nwa
                              sum = sum + temp(i,k)*da(k)
                        end do

                        tr(i,j) = 1048576*sum*z**3
                  end do
            end do

            write(9) amm,nar,nbr,waveFn_A(1:7)//'-'//waveFn_B(1:7)
            write(9) (eiga(i), i=1,nar)
            write(9) (eigb(i), i=1,nbr)

            do j=1,nar
                  write(9) (tr(i,j),i=1,nbr)
            end do

            write(*,'(A,I1,A,I1,A,D20.12)') '    tr(',neffb,',',neffa,') =',tr(neffb,neffa)
            rewind 9

            close(9)
            close(7)
            close(8)

            write(*,'(A,1PD24.16)') 'Check Ratio =',tr(neffb,neffa)/hterm(1)/1048576/z**3
      end if

      return
end subroutine breit


subroutine input_b

      !?Import CommonBlock
      use a1a_block
      use b1b_block

      !?Import Formatting
      use format

      implicit real*16(a-h,o-z)

      call diskr(db(neigb),nwb-neigb+1)

      !? CALCULATION OF SCREENED HYDROGENIC WAVEFUNCTIONS.
      if(izb.ne.0) then
            n = lrglb + neigb
            db0 = db(neigb)
            xx = 1
            l = 0

            do m=1,neigb
                  db(m) = db(neigb)*xx*2**(l+lrglb+2)*(izb-1)**(lrglb+1)*sqrt((izb-1)*fac(n+lrglb+1)/(izb*fac(n-lrglb)*fac(2*l+2)))/(fac(2*lrglb+2)*n*(l+1.q0)**(l+2)*(n*z)**(lrglb+1))
                  xx = -xx*2*(izb-1)*(neigb-m)/((2*lrglb+m+1)*m*n*izb)
            end do
      end if 

      sign = 1.

      if(db(1).lt.0.) then 
            sign = -1.
      end if 

      db0 = db0*sign

      do k=1,nwb
            db(k) = db(k)*sign
      end do

      dbn = db0/db(1)

      return
end subroutine input_b


subroutine diskr(db,n)
      implicit real*16(a-h,o-z)

      dimension db(n)

      read(8) (db(i),i=1,n)

      return
end subroutine diskr


subroutine input_a

      !?Import CommonBlock
      use a1a_block
      use b1b_block

      !?Import Formatting
      use format

      implicit real*16(a-h,o-z)

      call diskra(da(neiga),nwa-neiga+1)

      !? CALCULATION OF SCREENED HYDROGENIC WAVEFUNCTIONS.
      if(iza.ne. 0) then
            n = lrgla + neiga
            da0 = da(neiga)
            xx = 1
            l = 0

            do m=1,neiga
                  da(m) = da(neiga)*xx*2**(l+lrgla+2)*(iza-1)**(lrgla+1)*sqrt((iza-1)*fac(n+lrgla+1)/(iza*fac(n-lrgla)*fac(2*l+2)))/(fac(2*lrgla+2)*n*(l+1.q0)**(l+2)*(n*z)**(lrgla+1))
                  xx = -xx*2*(iza-1)*(neiga-m)/((2*lrgla+m+1)*m*n*iza)
            end do
      end if

      sign = 1.

      if(da(1).lt.0.) then 
            sign = -1.
      end if 

      da0 = da0*sign

      do k=1,nwa
            da(k) = da(k)*sign
      end do

      dan = da0/da(1)

      return
end subroutine input_a


subroutine diskra(da,n)
      implicit real*16(a-h,o-z)

      dimension da(n)

      read(7) (da(i),i=1,n)
      
      return
end subroutine diskra


