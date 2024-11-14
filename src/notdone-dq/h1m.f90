include 'cross2.f90'
include 'spin68.f90'

program h1mass

   !? PROGRAM TO EVALUATE THE FINITE MASS CORRECTION TO H1.  THE RESULTS ARE TO
   !? BE ADDED DIRECTLY TO THE OUTPUT FROM H1.

   !? Import Common Blocks
   use a1a_block
   use b1b_block
   use f1f_block
   use h1h_block
   use maxpow_block

   !? Import Formatting 
   use format

   !? Import routedOpen routine
   use wavExt 

   implicit real*16 (a-h,o-z)

   dimension hterm(10),htscr(10),nn(50)

   character :: cx
   character :: lines(100)
   character(len = 16) :: matout
   character(len = 52) :: title
   character(len = 12) :: nama,namb
   character(len = 51) :: line
   character(len = 12) :: fmt3a(3)
   character(len = 5) :: eigen
   character(len = 4) :: dotdmp

   logical :: lscreen

   integer :: endOfFile = 0
   integer :: errorChecker = 0

!? Open matl.dat for reading and open matl.out for writing
   call routedOpen(5,file='matl.dat',status='UNKNOWN')
   call routedOpen(4,file='matl.out',status='UNKNOWN')

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
            if(if.gt.1) backspace(1)
            read(1,fmt3a(if), iostat = errorChecker) wb(1), wb(2), ammb

            if(errorChecker .ne. 0) then
               cycle 
            end if

            go to 43
         end do

         write(*,*) 'B Format match not found'
         stop

      43 if(amm .ne. ammb) then
            write(*,'(2D20.10,A)') amm,ammb,' Masses not equal.'
            stop
         endif
      !
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

      do i=1,3
         hterm(i) = 2.*hterm(i)*amm
         htscr(i) = 2.*htscr(i)*amm
      end do

      do i=1,4
         hterm(i) = hterm(i)*1048576
         htscr(i) = htscr(i)*1048576
         htscr(i) = htscr(i)*0.5/(da0*db0)
      end do

      htscr(4) = htscr(4) - 1.
      hterm(4) = hterm(4) - 1.
      
      write(*, fmt1) (i,htscr(i),hterm(i),i=1,4,3)

      call routedOpen (1,file=matout,status='OLD')

      do while(endOfFile .eq. 0) 
         read(1,'(A14)',iostat = errorChecker) fnin

         if(errorChecker .gt. 0) then 
            print*, "An error has occurred, the program will now stop."
            stop 
         end if 
      end do 

      backspace 1

      write(1,'(A14,1X,A14,1X,A47)')waveFn_A,waveFn_B,title
      write(1, fmt39) (hterm(i),i=1,4,3)

      close(1,status='KEEP')

      write(*, fmt10) line
      write(*, fmt10)

      !? RESET INPUT FILE FOR NEXT CALCULATION.

      rewind 5

      do i=1,100
         ip = i-1
         read(5,'(A)',iostat = errorChecker) lines(i)

         if(errorChecker .lt. 0) then 
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
   end do
end program h1mass
!


subroutine breit(hterm,htscr)

   use a1a_block
   use b1b_block
   use c1c_block
   use d1d_block
   use f1f_block
   use maxpow_block

   implicit real*16 (a-h,o-z)

   dimension hterm(10),htscr(10),d1d2r(12),d1r2r(12),r1d2r(12),d1d2l(12),d1r2l(12),r1d2l(12),cplt6(12),cplt7(12),cplt8(12)

   integer p,q,s,a1,b1,c1,a4,b4,c4

   logical lgo,lcalc

!   FOR D, F STATES
!      DATA ISUMD,MP2D,MQ2D,MS2D,MSLD/0,2,2,1,1/
!   FOR D, F STATES

   data isumd,mp2d,mq2d,ms2d,msld/0,2,2,1,1/

!   FOR S STATES
!      DATA ISUMD,MP2D,MQ2D,MS2D,MSLD/0,1,2,0,0/
!   FOR P STATES
!      DATA ISUMD,MP2D,MQ2D,MS2D,MSLD/0,2,2,0,1/
!   POWERS OF R1 AND R12.  POWER OF R2 = 0.
!
!    IST CHOOSES THE COORDINATES OF ELECTRON 1 OR 2 AND 'I' SPECIFIES WHICH
!    SET OF COEFFICIENTS OF PL(COSTHETA12) IS TO BE USED.

      ermac = 1.1q0
 1000 ermac = ermac/2
      test = 1.1q0 + ermac
      if(test.gt.1.1q0) go to 1000
      kross = 9
      klog = -2
      kang = 1
      kang = 0
      mla = lrgla
      mlb = lrglb
      md = - mla + mlb
      do 6 i=1,10
      hterm(i) = 0.0
      htscr(i) = 0.0
    6 continue
      zz1 = 1./z
      e = wa(2)/z**2
!   SUMMATION OVER RIGHT HAND STATE DIRECT AND EXCHANGE TERMS.
      do 620 idx=1,2
!   SUMMATION OVER BLOCKS OF TERMS WITH THE SAME NON-LINEAR PAREMETERS.
      do 20 nhi=1,nbxa
      nbx1 = nblxa(nhi) + 1
      nbx2 = nblxa(nhi+1)
      do 20 nhj=1,nbxb
      nbx3 = nblxb(nhj) + 1
      nbx4 = nblxb(nhj+1)
!
!   THE OUTER LOOP CORRESPONDS TO THE RIGHT HAND STATE AND IS DENOTED BY
!   'A' OR'1' OR '2' IN VARIOUS PLACES.
!
      do 20 nbl1=nbx1,nbx2
      nb1 = nblka(nbl1) + 1
      nb2 = nblka(nbl1+1)
!
!   THE INNER LOOP CORRESPONDS TO THE LEFT HAND STATE AND IS DENOTED BY
!   'B' OR '3' OR '4' IN VARIOUS PLACES (3 MEANS ELECTRON 1 AND 4 MEANS ELECTRON
!    2).
!
      do 20 nbl3=nbx3,nbx4
      nb3 = nblkb(nbl3) + 1
      nb4 = nblkb(nbl3 + 1)
      lgo = .false.
      if(nbl1.eq.1.and.nbl3.eq.1) lgo = .true.
      if(izb.eq.0.and.nbl3.lt.nbl1) go to 20
!
!   SUMMATION OVER LEFT HAND STATE DIRECT AND EXCHANGE TERMS.
!   THE ORDERING IS D-D, E-E, D-E AND E-D.
!
!     DO 220 JDX=1,2
      jdx = 1
      ide =  jdx
      jde =  jdx
      if(idx.eq.2) ide = 3 - jde
      id = ide - 1
      jd = jde - 1
      idep = 3 - ide
      jdep = 3 - jde
      sgn1 = 1.
      if(nspna.eq.1) sgn1 = (-1.)**(ide-1)
      sgn = sgn1
      if(nspnb.eq.1) sgn = sgn1*(-1.)**(jde-1)
      l1 = la(ide,nb1)
      l2 = la(idep,nb1)
      l3 = lb(jde,nb3)
      l4 = lb(jdep,nb3)
      ll1 = l1*(l1+1)
      ll2 = l2*(l2+1)
      ll3 = l3*(l3+1)
      ll4 = l4*(l4+1)
      l13 = l1 + l3
      call cross(l3,l4,l1,l2,lrgla,ide,jde)
      mlp1 = mlt(1)
      lp1 = lt(1)
      lp4 = lt(4)
      mlp5 = mlt(5)
      lp5 = lt(5)
      llp1 = lp1*(lp1-1)
      ldel1 = ll3 - ll1 - llp1
      ldel2 = ll4 - ll2 - llp1
      ldel3 = ll1 - ll3 - llp1
      ldel4 = ll2 - ll4 - llp1
      do 72 l=mlp5,lp5,2
      d1d2r(l) = cplt(l,6)
      d1r2r(l) = cplt(l,7)
      r1d2r(l) = cplt(l,8)
      d1d2l(l) = cplt(l,9)
      d1r2l(l) = cplt(l,10)
   72 r1d2l(l) = cplt(l,11)
      if(kang.eq.0) go to 110
      write(*,991) l3,l4,l1,l2,lrgla,ide,nb3,nb4,nb1,nb2,izb
  991 format(20i3)
      write(*,111)
      do 210 ic=1,11
      if(lt(ic).le.0) go to 210
      ic2 = lt(ic)
      write(*,111) ic,(cplt(ic3,ic), ic3=  1,ic2)
  111 format(i4,1p,4d18.11/4x,4d18.11)
  210 continue
  110 if(jdx.eq.2) go to 77
      if(nbl1.gt.nbx1.or.nbl3.gt.nbx3) go to 77
      y11 = ya(ide,nb1) + yb(jde,nb3)
      y22 = ya(idep,nb1) + yb(jdep,nb3)
      kex = 0
      if(y22.gt.y11) then
        y22 = ya(ide,nb1) + yb(jde,nb3)
        y11 = ya(idep,nb1) + yb(jdep,nb3)
        kex = 1
      endif
      lmax = lp5
      mp1 = 4
      ms1 = 3
      msl = lt(1) + msld
      imin = -3 + lrgla + lrglb
      lsuma = 2*id
      lsumb = l2*(iabs(l1-l2)/2)
      ltota = lrgla + neiga - 1
      ltotb = lrglb + neigb - 1
      if((nb1.eq.1.or.nb3.eq.1).and.z.ne.1.) go to 78
      isum = ijkb(nbl3) + ijka(nbl1) + isumd
      mp2 = iib(nbl3) + (1-id)*iia(nbl1) + id*jja(nbl1) + mp1 + mp2d+ 1 - id
      mq2 = jjb(nbl3) + (1-id)*jja(nbl1) + id*iia(nbl1)+ mp1+mq2d-1+id
      ms2 = ms1 + kkb(nbl3) + kka(nbl1) + 2*lt(1) + ms2d
      go to 80
   78 if(nb3.eq.1.and.nb1.eq.1) go to 179
      if(nb3.eq.1) go to 79
      isum = ltota + ijkb(nbl3) + isumd
      mp2 = iib(nbl3) + id*ltota + mp1 + mp2d
      mq2 = jjb(nbl3) + (1-id)*ltota + mp1 + mq2d-1
      ms2 = ms1 + kkb(nbl3) + 2*lt(1) + ms2d + 1
      go to 80
   79 isum = ltotb + ijka(nbl1) + isumd
      mp2 = (1-id)*iia(nbl1) + id*(jja(nbl1)+ltotb) + mp1 + mp2d
      mq2 = (1-id)*(jja(nbl1)+ltotb) + id*iia(nbl1) + mp1 + mq2d-1
      ms2 = ms1 + kka(nbl1) + 2*lt(1) + ms2d
      go to 80
  179 isum = ltotb + ltota + isumd
      mp2 = (1-id)*ltota + id*(ltotb+ltota) + mp1 + mp2d
      mq2 = (1-id)*(ltotb+ltota) + id*ltota + mp1 + mq2d
      ms2 = ms1 + 2*lt(1) + ms2d
   80 continue
      if(id.eq.0.and.(nbx4.gt.nbx3.or.nbx2.gt.nbx1)) ms2 = ms2 + 2
      if(id.eq.0.and.(nbx4.gt.nbx3.and.nbx2.gt.nbx1)) ms2 = ms2 + 2
      if(id.eq.1.and.nbx2.gt.nbx1.and.l3.gt.0) ms2 = ms2 + 2
      if(id.eq.0.and.(nbx4.gt.nbx3.or.nbx2.gt.nbx1)) msl = msl+1
      if((nbx4.gt.nbx3.or.nbx2.gt.nbx1))lmax =min0(lrgla +2,lmax+1)
      if((nbx4.gt.nbx3.and.nbx2.gt.nbx1))lmax =min0(lrgla +2,lmax+1)
      if(msl.gt.12) stop 180
      if(kex.eq.1) then
        mtemp = mp2
        mp2 = mq2
        mq2 = mtemp
      endif
      if(nbl3.eq.nbl1) write(*,180) kex,mp1,mp2,mq2,ms1,ms2,isum,y11,y22,l3,l4,l1,l2,nbl1,nbl3
  180 format(' GENINT',i1,6i3,2f10.6,4i2,3x,2i2,$)
      call genint(id,nbl3,nbl1,lmax,fac,0.q0)
   77 delt12 = 0.0
      jdxx = jdx
      if(kex.eq.1) jdxx = 3 - jdx
!
!   SUMMATION OVER RIGHT HAND STATE BASIS FUNCTIONS.
!
      do 21 i=nb1,nb2
      a1 = pa(ide,i)
      b1 = pa(idep,i)
      c1 = sa(i)
      y1 = ya(ide,i)
      y2 = ya(idep,i)
      do 73 l=mlp5,lp5,2
      cplt6(l) = a1*b1*cplt(l,5) + a1*r1d2r(l) + b1*d1r2r(l) + d1d2r(l)
      cplt7(l) = -y2*(a1*cplt(l,5) + d1r2r(l))
   73 cplt8(l) = -y1*(b1*cplt(l,5) + r1d2r(l))
      do 68 k=6,11
   68 mlt(k) = mlp5
!
!   SUMMATION OVER LEFT HAND STATE BASIS FUNCTIONS.
!
  224 do 21 j=nb3,nb4
      if(izb.eq.0.and.j.lt.i) go to 21
      psi2 = 2.*da(i)*db(j)*sgn
      if(izb.eq.0.and.j.ne.i) psi2 = psi2*2.
      a4 = pb(jde,j)
      b4 = pb(jdep,j)
      c4 = sb(j)
      y3 = yb(jde,j)
      y4 = yb(jdep,j)
      p = pa(ide,i) + pb(jde,j)
      q = pa(idep,i) + pb(jdep,j)
      s = sa(i) + sb(j)
      y_ = y1 + y3
      yp_ = y2 + y4
      y1234 = y1*y2 + y3*y4
      cc1 = c1*(c1+1) + c4*(c4+1)
      do 74 l=mlp5,lp5,2
      cplt(l,6) = cplt6(l) + a4*b4*cplt(l,5) + a4*r1d2l(l) + b4*d1r2l(l)+ d1d2l(l)
      cplt(l,7) = cplt7(l) - y4*(a4*cplt(l,5) + d1r2l(l))
      cplt(l,8) = cplt8(l) - y3*(b4*cplt(l,5) + r1d2l(l))
      cplt(l,9) = e*cplt(l,6) + cplt(l,7) + cplt(l,8)
      cplt(l,10) = e*cplt(l,7) + y1234*cplt(l,5)
   74 cplt(l,11) = e*cplt(l,8) + y1234*cplt(l,5)
      do 69 k=6,11
      lt(k) = lp5
   69 if(abs(cplt(lp5,k)).lt.1.d-10) lt(k) = lp5 - 2
      if(abs(cplt(lt(6),6)).lt.1.d-10) lt(6) = lt(6) - 2
      if(abs(cplt(lt(6),6)).lt.1.d-10) lt(6) = lt(6) - 2
      if(abs(cplt(lt(6),6)).lt.1.d-10) lt(6) = lt(6) - 2
!
      sum = 0.
      sumz = 0.
      sum4 = 0.
      cjp = c1
      cip = c4
      cjl = c1
      cil = c4
      is = s - 1
      if(is.eq.0) go to 67
      cjl = cjl/is
      cil = cil/is
   67 if(s.eq.0) go to 66
      cjp = cjp/s
      cip = cip/s
   66 qj = 0.5*(a1*(p+1)*2 + ldel1)
      qi = 0.5*(a4*(p+1)*2 + ldel3)
      qd11 = cjp*qj + cip*qi
      qd12 = qd11 - a1*cjp - a4*cip
      qd14 = cjl*qj + cil*qi
      qj = 0.5*(b1*(q+1)*2 + ldel2)
      qi = 0.5*(b4*(q+1)*2 + ldel4)
      qe11 = cjp*qj + cip*qi
      qe13 = qe11 - b1*cjp - b4*cip
      qe14 = cjl*qj + cil*qi
      qj = -(a1*y_ + y1*(p+2))
      qi = -(a4*y_ + y3*(p+2))
      qd21 = cjp*qj + cip*qi
      qd22 = qd21 + cjp*y1 + cip*y3
      qd24 = cjl*qj + cil*qi
      qj = -(b1*yp_ + y2*(q+2))
      qi = -(b4*yp_ + y4*(q+2))
      qe21 = cjp*qj + cip*qi
      qe23 = qe21 + cjp*y2 + cip*y4
      qe24 = cjl*qj + cil*qi
      qd31 = y_*(y1*cjp + y3*cip)
      qd34 = y_*(y1*cjl + y3*cil)
      qe31 = yp_*(y2*cjp + y4*cip)
      qe34 = yp_*(y2*cjl + y4*cil)
!
      ov = spl(p,q,s,1)
      sum = qd12*spl(p-3,q,s,1) + qd11*spl(p-2,q-1,s,1) + qe13*spl(p,q-3,s,1) + qe11*spl(p-1,q-2,s,1)
      sum = sum + (e*qd11 + qd22)*spl(p-2,q,s,1) + (e*qe11 + qe23)*spl(p,q-2,s,1) + (qd21 + qe21)*spl(p-1,q-1,s,1) + (e*qd21 + qd31 + qe31)*spl(p-1,q,s,1) + (e*qe21 + qd31 + qe31)*spl(p,q-1,s,1) + e*(qd31 + qe31)*ov
      if(s.eq.0) go to 75
      sumz = qd14*splg(p-2,q,is,1,is) + qe14*splg(p,q-2,is,1,is) + qd24*splg(p-1,q,is,1,is) + qe24*splg(p,q-1,is,1,is) + (qd34 + qe34)*splg(p,q,is,1,is) - cc1*spl(p,q,s-3,1)
      sum = sum - cc1*(spl(p-1,q,s-2,1) + spl(p,q-1,s-2,1) + e*spl(p,q,s-2,1))
   75 sumz = sumz + y1234*spl(p,q,is,5) + spl(p-1,q-1,is,6) + spl(p-1,q,is,7) + spl(p,q-1,is,8)
      if(abs(sumz).le.1.d70) go to 555
      go to 555
      write(*,132)3,i,j,p,q,s,sumz
      do 556 kk=5,8
      write(*,557)kk,mlt(kk),lt(kk),(cplt(ll,kk),ll=mlt(kk),lt(kk),2)
  557 format(3i3,6d11.4)
  556 continue
  555 continue
      if(s.eq.0.and.mlp5.gt.1) go to 76
      sum = sum + spl(p-2,q-1,s,6) + spl(p-1,q-2,s,6) + spl(p-2,q,s,7)   + spl(p,q-2,s,8)
      sum = sum + spl(p-1,q-1,s,9) + spl(p-1,q,s,10) + spl(p,q-1,s,11) + e*y1234*spl(p,q,s,5)
   76 if(lp4.eq.0.or.s.eq.0) go to 65
      sum4 = e*(spl(p-2,q,s,4) + spl(p,q-2,s,4)) + spl(p-3,q,s,4) + spl(p-2,q-1,s,4) + spl(p,q-3,s,4) + spl(p-1,q-2,s,4)
      sumz = sumz - (cil+cjl)*(splg(p-2,q,is,4,is) + splg(p,q-2,is,4,is))
      if(abs(sumz).gt.1.d70) write(*,132)7,i,j,p,q,s,sumz
   65 sum = sum - sum4 - zz1*sumz
  132 format(' TRAP', 6i4,2d12.5)
      hterm(1) = hterm(1) + sum*psi2
      if(lgo) htscr(1) = htscr(1) + sum*psi2
  888 format(2i3,2x,3i3,d20.12)
      hterm(4) = hterm(4) + ov*psi2
      if(lgo) htscr(4) = htscr(4) + ov*psi2
   21 continue
  220 continue
   20 continue
  620 continue
  100 return
      end

