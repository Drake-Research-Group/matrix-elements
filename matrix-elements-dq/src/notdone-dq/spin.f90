!?   MODIFIED MARCH 28/21 TO USE KONTROL TO SPECIFY WHICH OPERATOR FOR OUTPUT.
!?   KONTROL = 1 FOR delta(r1)
!?   KONTROL = 2 FOR delta(r12)
!?   KONTROL = 5 FOR SPIN-SPIN

include 'cross1.f90'
include 'spin68.f90'

program spin

   !? Import Common Blocks
   use a1a_block
   use b1b_block
   use f1f_block
   use h1h_block
   use maxpow_block

   !? Import Formatting
   use format 

   !? Import OpenWavefunctions routine
   use wavExt

   implicit real*16 (a-h,o-z)

   logical lwrite

   !? More Initializations
   dimension hterm(20),htscr(20),nn(50)

   character*16 matout
   character title*52,nama*12,namb*12,line*51,fmt3a(3)*12, eigen*5, dotdmp*4
   character(len = 57) :: lines(100)

   logical lscreen

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

      leigb = .false.

      read(1,'(1X,A5)',end=28) eigen

      if(eigen.eq.'EIGEN') then 
         read(1, fmt3a(if)) (eigb(i),i=1,nwb-neigb+1)

         leigb = .true.
      end if 


   28 close(1,status='KEEP')

      ! OPEN .DMP FILE WITH SQUaRE ARRAY OF WAVEFUNCTION COEFFICIENTS.
      if(leigb) then
         write(*,'(2A)') 'spin OPENING ',waveFn_B(1:7)//dotdmp
         call routedOpen(8,file=waveFn_B(1:7)//dotdmp,form='UNFORMATTED',status='OLD')
      endif

   151 lpar = 0

      close(1,status='KEEP')

      call breit(hterm,htscr,amm,lpar,kontrol)

      do i=1,20
         hterm(i) = hterm(i)*1048576
         htscr(i) = htscr(i)*1048576
         htscr(i) = htscr(i)*0.5/(da0*db0)
      end do 

      if(hterm(1).ne.0.) then 
         hterm(1) = 2.*hterm(1) - 1.
      end if 

      if(hterm(16).ne.0.) then 
         hterm(16) = 2.*hterm(16) - 1.
      end if 

      i1 = 1
      i2 = 18

      if(lrgla.ne.lrglb) then
        i1 = 13
        i2 = 18
      endif

      write(*,20) ((i+2)/3,hterm(i),i=i1,i2,3)
   20 format(' TERM',i2,' =',d20.12)

      open(1,file=matout,status='OLD',err=226)

      go to 26

  226 open(1,file=matout,status='NEW')

      write(1, fmt10) line
      write(1,'(a,a12,a)') '''',matout,''''

   26 read(1,'(A14,1X,A14)',end=21) fnwv1,fnwv2

      go to 26

   21 backspace 1

      if(lpar.eq.0) then 
         write(1,'(A14,1X,A14,1X,A52)') waveFn_A,waveFn_B,title
      end if 

      if(lpar.eq.1) then 
         write(1,'(A14,1X,A14,'' PARONIC'')')waveFn_A,waveFn_B
      end if 

      write(1,23) (hterm(i),i=1,17,3)
   23 format(1x,1p,d18.11)

  521 close(1,status='KEEP')

      write(*, fmt10) line
      write(*, fmt10)

      ! RESET INPUT FILE FOR NEXT CALCULATION.
      cycle
      rewind 5

      do i=1,100
         ip = i-1
         read(5,'(A)',end=901) lines(i)
      end do 

  901 rewind 5

      write(5,'(A)') lines(1)
      write(5,'(A16)') lines(2)

      do ii=4,ip
         write(5,'(A32)') lines(ii)
      end do 

      write(5,'(A32)') lines(3)

      rewind 5

      read(5,'(/)')

   end do
end program spin

subroutine breit(hterm,htscr,amm,lpar,kontrol)

      !? Import Common Block
      use a1a_block
      use b1b_block
      use c1c_block
      use d1d_block
      use f1f_block
      use h1h_block
      use maxpow_block

      implicit real*16 (a-h,o-z)

      dimension hterm(20),htscr(20),d1d2(12),r1d2(12),r2d1(12),cpltm(12,3),trb(3601,5),pert(3601,5),tr(3601,3601,5),temp(3601,3601)
      integer p,q,s,a1,b1,c1,a4,b4,c4,pam,qam
      logical lgo,lcalc
      character op*3,title*15,dotdmp*4

      !? FOR P-STATES
      !? DATA ISUMD,MP2D,MQ2D,MS2D,MSLD/-2,4,3,0,1/

      !? FOR D-STATES
      data isumd,mp2d,mq2d,ms2d,msld/-2,3,2,0,2/

      !? PARONIC VERSION.  REMOVE NEXT TWO EXECUTABLE STATEMENTS AND REVERSE THE
      !? SIGN OF SGN TO OBTAIN THE NORMAL VERSION..
      !? IST CHOOSES THE COORDINATES OF ELECTRON 1 OR 2 AND 'I' SPECIFIES WHICH
      !? SET OF COEFFICIENTS OF PL(COSTHETA12) IS TO BE USED.

      ermac = 1.1d0
 1000 ermac = ermac/2

      !? USE .DML TO DENOTE THE FINITE MASS CASE.
      dotdmp = '.DMP'
      if(amm.ne.0) dotdmp = '.DML'
      test = 1.1d0 + ermac
      if(test.gt.1.1d0) go to 1000
      if(lpar.eq.1) nspna = 1 - nspna
      if(lpar.eq.1) nspnb = 1 - nspnb
      kross = 9
      klog = -2
      kang = 0
      sq3 = dsqrt(3.d0)
      mla = lrgla
      mlb = lrglb
      mla = min0(lrgla,lrglb)
      mlb = min0(lrglb,lrgla)
      md = - mla + mlb
      do 6 i=1,18
      hterm(i) = 0.0
      htscr(i) = 0.0
    6 continue
      do 7 k=1,5
      do 7 i=1,nwa
      do 7 j=1,nwb
      trb(j,k) = 0
    7 tr(j,i,k) = 0
      if(lrgla.eq.0.and.lrglb.eq.0) go to 71
      ss5 = (-1)**(lrglb-mlb)*dsqrt(15.d0/8.)
      fj = f3j(lrglb,-mlb,2,md,lrgla,mla)
      if(fj.eq.0.0) go to 103
      ss5 = ss5/fj
      if(lrgla.eq.0.or.lrglb.eq.0) go to 71

      !? CALCULATES THE REDUCED MATRIX ELEMENT (L'||L||L)(S'||S||S).  MULTIPLY RESULT
      !? BY F6J(J,1,S',L,L',S)*(-1)**(L+S'+J)*ALPHA**2.
   72 continue
      fj = f3j(lrglb,-mlb,1,md,lrgla,mla)*(-1.)**(lrglb-mlb)
      if(fj.eq.0.0) go to 103
      if(nspna+nspnb-1) 76,75,76
   75 so =-z*sq3/fj/2.
      if(nspna.eq.0) so = - so
      soo = so/z
      go to 71
   76 so = z*sq3    /fj/dsqrt(2.d0)
      soo = -3.*so/z

      !? THE OUTER LOOP CORRESPONDS TO THE RIGHT HAND STATE AND IS DENOTED BY
      !? 'A' OR'1' OR '2' IN VARIOUS PLACES.

   !   SUMMATION OVER RIGHT HAND STATE DIRECT AND EXCHANGE TERMS.
   71 do 620 idx=1,2
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
      do 220 jdx=1,2
   !      JDX = 1
      ide =  jdx
      jde =  jdx
      if(idx.eq.2) ide = 3 - jde
      id = ide - 1
      jd = jde - 1
      idep = 3 - ide
      jdep = 3 - jde
   !      SGN1 = 1.
   !      IF(NSPNA.EQ.1) SGN1 = (-1.)**(IDE-1)
   !      SGN = SGN1
   !      IF(NSPNB.EQ.1) SGN = SGN1*(-1.)**(JDE-1)
      sgn1 = 1.
      if(nspna.eq.1-lpar) sgn1 = (-1.)**(ide-1)
      sgn = sgn1
      if(nspnb.eq.1-lpar) sgn = sgn1*(-1.)**(jde-1)
      sgn2 = sgn
      if(nspna.ne.nspnb.and.nspna.eq.1-lpar) sgn2 = sgn*(-1.)**(ide-1)
      if(nspna.ne.nspnb.and.nspnb.eq.1-lpar) sgn2 = sgn*(-1.)**(jde-1)
      l1 = la(ide,nb1)
      l2 = la(idep,nb1)
      l3 = lb(jde,nb3)
      l4 = lb(jdep,nb3)
      ll1 = l1*(l1+1)
      ll3 = l3*(l3+1)
      ll4 = l4*(l4+1)
      call cross(l3,l4,l1,l2,lrglb,mlb,lrgla,mla,ide,jde)
   !      CALL CROSS(L3,L4,L1,L2,LRGLB,0  ,LRGLA,0  ,IDE,JDE)
      mlp2 = mlt(2)
      lp2 = lt(2)
      lll = max0(lt(14),lt(15),lt(16))

   !      lt(4) = 0
   !      lt(5) = 0
   !      lt(6) = 0
   !      DO 555 L=1,LLL
   !      cplt(L,4) = cplt(L,14)
   !      IF(DABS(cplt(L,4)).GT.1.D-10) lt(4) = L
   !      cplt(L,5) = cplt(L,15)
   !      IF(DABS(cplt(L,5)).GT.1.D-10) lt(5) = L
   !      cplt(L,6) = cplt(L,16)
   !  555 IF(DABS(cplt(L,6)).GT.1.D-10) lt(6) = L
   !      mlt(4) = mlt(14)
   !      mlt(5) = mlt(15)
   !      mlt(6) = mlt(16)
   !
   !   SAVE DEL1.DEL2, R2.DEL1 AND R1.DEL2 ANGULAR COEFS.
      do 54 l=1,lp2
      d1d2(l) = cplt(l,6)
      r2d1(l) = cplt(l,10)
   54 r1d2(l) = cplt(l,9)
      mlp3 = mlt(3) - 1
      if(mlp3.lt.1) mlp3 = mlp3 + 2
      lp3 = lt(3) - 1
      mlt(3) = mlp3
      lt(3) = lp3
      do 50 ll=mlp3,lp3,2
      l = lp3 - ll + mlp3
   !   FORM COS**2(THETA) - 1 AND REDUCE BY RECURSION RELATION.
      if(l.eq.1.and.abs(cplt(l+1,3)).gt.1.d-10) stop 1
      cplt(l,3) = (2*l-1)*cplt(l+1,3)
      if(l.eq.1) go to 50
      cplt(l-1,3) = cplt(l-1,3) + cplt(l+1,3) - cplt(l-1,1)
      if(abs(cplt(l-1,3)).lt.1.d-10) cplt(l-1,3) = 0.
   50 continue
      do 52 ll=mlp3,lp3,2
      l = lp3 - ll + mlp3
      if(cplt(l,3).ne.0.) go to 53
   52 lt(3) = lt(3) - 2
   53 l10 = max0(lt(11)-1,lt(13))
      mlp1 = mlt(1)
      lp1 = lt(1)
      do 60 ll=mlp1,lp1,2
      l = lp1 - ll + mlp1
   !   REDUCE cplt(L,11) BY RECURSION RELATION.
      if(l.eq.1.and.abs(cplt(l+1,11)).gt.1.d-10) stop 2
      cplt(l,11) = (2*l-1)*cplt(l+1,11)
      if(l.eq.1) go to 60
      cplt(l-1,11) = cplt(l-1,11) + cplt(l+1,11)
      if(abs(cplt(l-1,11)).lt.1.d-10) cplt(l-1,11) = 0.
   60 continue
      ml10 = 0
      if(mlt(11).gt.0) ml10 = mlt(11) + 1
      if(mlt(13).gt.0) ml10 = min0(mlt(13),ml10)
      if(mlt(13).gt.0.and.ml10.eq.0) ml10 = mlt(13)
      mlt(10) = ml10
      lt(10) = l10
      if(lt(11).eq.0) go to 44
      mlt(11) = mlt(11) + 1
      lt(11) = lt(11) -1
   44 l11 = lt(11)
      lss = max0(lt(14)-1,lt(15)-1,lt(16)-2)
      if(lss.eq.1) lss = 3
      mlt(7) = 0
      lt(7) = 0
      mlt(14) = 0
      lt(14) = 0
      mlt(15) = 0
      lt(15) = 0
      mlt(16) = mlt(16) - 1
      lt(16) = lt(16) - 1
      cplt16 = cplt(1,16)
   !   CORRECT FOR NON-ZERO cplt(1,16).
      if(cplt16.ne.0) then
        cplt(1,16) = 0
        cplt(2,14) = cplt(2,14) - 1.5d0*cplt16
        cplt(2,15) = cplt(2,15) - 1.5d0*cplt16
        cplt(3,16) = cplt(3,16) - 2*cplt16
        mlt(16) = 2
      endif
      do 62 l=1,lss+1
      cplt(l,14) = -3.*(cplt(l,14) - 0.5*cplt(l+1,16))
      cplt(l,15) = -3.*(cplt(l,15) - 0.5*cplt(l+1,16))
      cplt(l,7) = 3.*(l-1)*cplt(l+1,16)
      cplt(l,16) = -1.5*cplt(l+1,16) + cplt(l,1)
      if(l.eq.2) cplt(l,16) = cplt(l,16) - 9*cplt16/2
      if(abs(cplt(l,7)).gt.1.d-10.and.mlt(7).eq.0) mlt(7) = l
      if(abs(cplt(l,7)).gt.1.d-10) lt(7) = l
      if(abs(cplt(l,14)).gt.1.d-10.and.mlt(14).eq.0) mlt(14) = l
      if(abs(cplt(l,14)).gt.1.d-10) lt(14) = l
      if(abs(cplt(l,15)).gt.1.d-10.and.mlt(15).eq.0) mlt(15) = l
      if(abs(cplt(l,15)).gt.1.d-10) lt(15) = l
   62 continue
      mlp14 = mlt(14) - 1
      if(mlp14.lt.1) mlp14 = mlp14 + 2
      lp14 = lt(14) - 1
      do 63 ll=mlp14,lp14,2
      l = lp14 - ll + mlp14
   !   REDUCE cplt(L,14) BY RECURSION RELATION.
      if(l.eq.1.and.abs(cplt(l+1,14)).gt.1.d-10) stop 63
      cplt(l,14) = (2*l-1)*cplt(l+1,14)
      if(l.eq.1) go to 63
      cplt(l-1,14) = cplt(l-1,14) + cplt(l+1,14)
      if(abs(cplt(l-1,14)).lt.1.d-10) cplt(l-1,14) = 0.
   63 continue
      if(lp14.gt.0) lt(14) = lp14
      if(lp14.gt.0) mlt(14) = mlp14
      mlp15 = mlt(15) - 1
      if(mlp15.lt.1) mlp15 = mlp15 + 2
      lp15 = lt(15) - 1
      do 64 ll=mlp15,lp15,2
      l = lp15 - ll + mlp15
   !   REDUCE cplt(L,15) BY RECURSION RELATION.
      if(l.eq.1.and.abs(cplt(l+1,15)).gt.1.d-10) stop 64
      cplt(l,15) = (2*l-1)*cplt(l+1,15)
      if(l.eq.1) go to 64
      cplt(l-1,15) = cplt(l-1,15) + cplt(l+1,15)
      if(abs(cplt(l-1,15)).lt.1.d-10) cplt(l-1,15) = 0.
   64 continue
      if(lp15.gt.0) lt(15) = lp15
      if(lp15.gt.0) mlt(15) = mlp15
   66 if(kang.eq.0) go to 110
      write(4,991) lss,l3,l4,l1,l2,lrglb,mlb,lrgla,mla,ide,jde,nb3,nb4,nb1,nb2,izb
  991 format(20i3)
      write(4,111)
      do 210 ic=1,16
      if(lt(ic).le.0) go to 210
      ic2 = lt(ic)
      write(4,111) ic,mlt(ic),(cplt(ic3,ic), ic3=  1,ic2)
  111 format(i2,i2,1p,4d18.11/4x,4d18.11)
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
      lmax = lt(2)
      mp1 = 4
      ms1 = 4
      msl = lt(1) + msld
      imin = -3 + lrgla + lrglb
      lsuma = 0
      lsumb = l2*(iabs(l1-l2)/2)
      ltota = lrgla + neiga - 1
      ltotb = lrglb + neigb - 1
      if((nb1.eq.1.or.nb3.eq.1).and.z.gt.1) go to 78
      isum = ijkb(nbl3) + ijka(nbl1) + isumd
      mp2 = iib(nbl3) + (1-id)*iia(nbl1) + id*jja(nbl1) + mp1 + mp2d + 1 - id
      mq2 = jjb(nbl3) + (1-id)*jja(nbl1) + id*iia(nbl1)+ mp1+mq2d + id
      ms2 = ms1 + kkb(nbl3) + kka(nbl1) + 2*lt(1) + ms2d + id
      go to 80
   78 if(nb3.eq.1.and.nb1.eq.1) go to 179
      if(nb3.eq.1) go to 79
   !   LIMITS FOR HYDROGENIC '1,N' TERMS.
      isum = ltota + ijkb(nbl3) + isumd
      mp2 = iib(nbl3) + id*ltota + mp1 + mp2d
      mq2 = jjb(nbl3) + (1-id)*ltota + mp1 + mq2d
      ms2 = ms1 + kkb(nbl3) + 2*lt(1) + ms2d + 1
      go to 80
   !   LIMITS FOR HYDROGENIC 'N,1' TERMS.
   79 isum = ltotb + ijka(nbl1) + isumd
      mp2 = (1-id)*iia(nbl1) + id*(jja(nbl1)+ltotb) + mp1 + mp2d
      mq2 = (1-id)*(jja(nbl1)+ltotb) + id*iia(nbl1) + mp1 + mq2d+1
      ms2 = ms1 + kka(nbl1) + 2*lt(1) + ms2d
      if(izb.ne.0) mq2 = mq2 + id*(ltotb)
      go to 80
  179 isum = ltotb + ltota + isumd
      mp2 = (1-id)*ltota + id*(ltotb+ltota) + mp1 + mp2d
      mq2 = (1-id)*(ltotb+ltota) + id*ltota + mp1 + mq2d
      ms2 = ms1 + 2*lt(1) + ms2d + id
   80 continue
      if(id.eq.0.and.(nbx4.gt.nbx3.or.nbx2.gt.nbx1)) ms2 = ms2 + 2
      if(id.eq.0.and.(nbx4.gt.nbx3.and.nbx2.gt.nbx1)) ms2 = ms2 + 2
      if(id.eq.1.and.(nbx2.gt.nbx1.or.nbx4.gt.nbx3).and.(l2.gt.0.or.l3.gt.0)) ms2 = ms2 +2
      if(id.eq.0.and.(nbx4.gt.nbx3.or.nbx2.gt.nbx1)) msl = msl + 1
      if((nbx4.gt.nbx3.or.nbx2.gt.nbx1))lmax =min0(lrgla +2,lmax+1)
      if((nbx4.gt.nbx3.and.nbx2.gt.nbx1))lmax =min0(lrgla +2,lmax+1)
      if(amm.ne.0.) ms2 = ms2 + 1
      if(msl.gt.12) stop 180
      if(kex.eq.1) then
        mtemp = mp2
        mp2 = mq2
        mq2 = mtemp
      endif
  180 format(' GENINT',i1,6i3,2f10.6,4i2,3x,2i2,$)
      mq2 = mq2 + 2
      mp2 = mp2 + 2
      isum = isum + 2
      if(kontrol.eq.1.or.kontrol.eq.2) go to 77
      if(nbl3.eq.nbl1) write(*,180) kex,mp1,mp2,mq2,ms1,ms2,isum,y11,y22,l3,l4,l1,l2,nbl1,nbl3
      call genint(id,nbl3,nbl1,lmax,fac,0.q0)
   77 delt12 = 0.0
      jdxx = jdx
      if(kex.eq.1) jdxx = 3 - jdx
      m1 = mlt(1)
      m2 = lt(1)
      do 59 j=m1,m2
   59 delt12 = delt12 + cplt(j,1)
   !
   !   SUMMATION OVER RIGHT HAND STATE BASIS FUNCTIONS.
   !
      do 21 i=nb1,nb2
      a1 = pa(ide,i)
      b1 = pa(idep,i)
      c1 = sa(i)
      fc1 = c1
      y1 = ya(ide,i)
      y2 = ya(idep,i)
      cc1 = c1*(c1 - 2)
      do 55 l=mlp2,lp2,2
      cplt(l,9) = b1*cplt(l,2) + r1d2(l)
      cplt(l,4) = b1*r2d1(l) + d1d2(l) - cplt(l,9)*(1. + 0.5*c1)
      cpltm(l,2) = -y2*r2d1(l) + y2*cplt(l,2)
   55 cpltm(l,3) = c1*(r2d1(l) + 0.5*cplt(l,9) - cplt(l,2))
      mlt(4) = mlp2
      lt(4) = lp2
      mlt(9) = mlp2
      lt(9) = lp2
      if(abs(cplt(lp2,9)).lt.1.d-08) lt(9) = lt(9) - 2
      if(abs(cplt(lp2,4)).lt.1.d-08) lt(4) = lt(4) - 2

      !? SUMMATION OVER LEFT HAND STATE BASIS FUNCTIONS.
      do 21 j=nb3,nb4
      if(izb.eq.0.and.j.lt.i) go to 21
      psi2 = da(i)*db(j)*sgn
      if(izb.eq.0.and.j.ne.i) psi2 = psi2*2.
      y3 = yb(jde,j)
      y4 = yb(jdep,j)
      p = pa(ide,i) + pb(jde,j)
      q = pa(idep,i) + pb(jdep,j)
      s = sa(i) + sb(j)
      y_ = y1 + y3
      yp_ = y2 + y4
      c14 = sa(i)*sb(j)
      if(lrgla.ne.lrglb) go to 5

      !? CALCULATION OF PI*DELTA(R1) MATRIX ELEMENT
      sum = 0.
      if(l1+l3.ne.0) go to 81
      if(p.ne.0) go to 81
      ipow = q + s + 3
      sum = fac(ipow)/(yp_**ipow*4.*1048576)
      hterm(16) = hterm(16) + sum*psi2
      trb(j,1) = trb(j,1) + sum*psi2/db(j)
      tr(j,i,1) = tr(j,i,1) + sum*psi2/(da(i)*db(j))

   81 if(nspna.ne.nspnb) go to 3
      tx1 = spl(p-2,q,s,1)
      tx2 = spl(p+1,q,s-3,1)
      tx3 = spl(p,q+1,s-3,2)
      tx23 = tx2-tx3
      is = s - 2
      cc = c14
      if(is.ne.0) cc = c14/is
      tx4 = 0.
      if(cc.ne.0.) tx4 = - cc*splg(p-2,q+1,is,3,is)
      cm = sa(i) - sb(j)
      if(s.gt.0) cm = cm/s
      cc = 0.5*(ll1*(1.-cm) + ll3*(1.+cm))
      tx5 = 0.
      if(cc.ne.0.) tx5 = tx5 + cc*spl(p-3,q,s,1)
      sum = tx1 - (tx23)/z - tx4 - tx5
   !
   !   ADD MASS POLARIZATION CONTRIBUTION IF AMM # 0.
      if(amm.eq.0.) go to 49
      do 43 k=2,3
      do 43 l=mlp2,lp2,2
   43 cplt(l,k+3) = cpltm(l,k)
      if(c1.eq.0) go to 56
      do 57 l=mlp2,lp2,2
      cplt(l,5) = cplt(l,5) - y2*c1*cplt(l,3)/s
   57 if(is.ne.0) cplt(l,6) = cplt(l,6) + cc1*cplt(l,3)/is
   56 do 58 k=5,6
      mlt(k) = mlp2
      lt(k) = lp2
   58 if(abs(cplt(lp2,k)).lt.1.d-08) lt(k) = lt(k) - 2
      txm1 = spl(p-2,q-1,s,4)
      txm2 = spl(p-2,q,s,5)
      txm = txm1 + txm2
  668 if(c1.gt.0) txm = txm + spl(p-2,q+1,is,6) + c1*(0.5*spl(p,q-1,is,9) - (b1+1)*spl(p-1,q,is,1))
      if(is.eq.0.and.cc1.ne.0) txm = txm + cc1*fpl(p-2,q+1,0,3,0)
      sum = sum + txm*amm
      hterm(17) = hterm(17) + amm*txm*psi2*0.5
   !      IF(LGO) HTSCR(7) = HTSCR(7) + AMM*TXM*PSI2*0.5
   !
   49 continue
   !      IF(abs(SUM).GT.1.D70) WRITE(*,666)6,I,J,P,Q,S,SUM
   !      SUM = SUM*SGN2
      if(lgo) htscr(1) = htscr(1) + sum*psi2*0.5
      hterm(1) = hterm(1) + sum*psi2*0.5
   !   81 CONTINUE  81 appears above.
   !
   !   CALCULATION OF PI*DELTA(R12) MATRIX ELEMENT
   !
    2 sum = 0.
      if(s.ne.0) go to 98
      ipow = p + q + 3
      sum = delt12*fac(ipow)/((y_+yp_)**ipow*2.*1048576)
   !      IF(LGO) HTSCR(2) = HTSCR(2) + SUM*PSI2
   97 hterm(4) = hterm(4) + sum*psi2
      trb(j,2) = trb(j,2) + sum*psi2/db(j)
      tr(j,i,2) = tr(j,i,2) + sum*psi2/(da(i)*db(j))
   98 dr12 = sum
      if(lrgla.eq.0) go to 28
   !
   !   CALCULATION OF SPIN-ORBIT MATRIX ELEMENT.
   !
    3 sum = 0.0
      if(l10.lt.1) go to 4
      rs = 0.
      if(s.gt.0) rs = fc1/s
      do 61 ll=ml10,l10,2
   61 cplt(ll,10) = cplt(ll,13) - rs*cplt(ll,11)
      sum = spl(p-3,q,s,10)
      hterm(7) = hterm(7) + sum*so*psi2
   !      IF(abs(SUM).GT.1.D70) WRITE(*,666)3,I,J,P,Q,S,SUM
   !
   !   CALCULATION OF SPIN-OTHER-ORBIT MATRIX ELEMENT.
   !
    4 sum = 0.0
      if(nspna+nspnb.eq.0) go to 28
      if(lt(13).gt.0) sum = spl(p,q,s-3,13)
      if(l11.eq.0) go to 45
      is = s - 1
      rs = 1.
      if(is.ne.0) rs = 1.d0/is
      sum = sum - y1*rs*splg(p-1,q,is,11,is)
      if(a1.ne.0) sum = sum + a1*rs*splg(p-2,q,is,11,is)
   45 if(lt(12).gt.0) sum = sum - spl(p-1,q+1,s-3,12)
      hterm(10) = hterm(10) + sum*soo*psi2
   !      IF(abs(SUM).GT.1.D70) WRITE(*,666)4,I,J,P,Q,S,SUM
   !
   !   CALCULATION OF MAGNETIC SPIN-SPIN INTERACTION.
   !
    5 sum = 0.0
      if(jdx.eq.2) go to 28
      if(nspna.ne.nspnb) go to 28
   !
   !   32 SUM = SPL(P,Q,S-3,1)
   !      SUM1 = SPL(P+2,Q,S-5,4) + SPL(P,Q+2,S-5,5)
   !     1    - SPL(P+1,Q+1,S-5,6)
   !      STOTX = 2.*(SUM - 3.*SUM1)*SS5*PSI2
   !      STOTX = 2.*(SUM - 3.*SUM1)*SS5*1048576
   !
      iss = s - 3
      css = 1.
      if(iss.ne.0) css = 1.d0/iss
      sum = spl(p,q,iss,16) + css*(splg(p+1,q-1,iss,14,iss) + splg(p-1,q+1,iss,15,iss))
      sum = sum + css*splg(p,q,iss,7,iss)
      stot = 2.*sum*ss5*psi2
      trb(j,5) = trb(j,5) - stot/db(j)
      tr(j,i,5) = tr(j,i,5) - stot/(da(i)*db(j))
  556 hterm(13) = hterm(13) - stot
  666 format(' T =',6i4,1p3d15.7/28x,3d15.7)
   28 continue
   21 continue
  220 continue
   20 continue
  620 continue
      if(leigb) then
         op = ' '
         if(kontrol.eq.1) then
            op='dl1'
            kterm = 16
         elseif(kontrol.eq.2) then
            op = 'd12'
            kterm = 4
         elseif(kontrol.eq.5) then
            op = 'ss5'
            kterm = 13
         else
            write(*,*) 'ERROR FROM SPIN, KONTROL = ',kontrol
            stop
         endif
         write(*,'(2A)') 'spin OPENING ', waveFn_A(1:4)//waveFn_B(1:4)//op//dotdmp
         open(9,file=waveFn_A(1:4)//waveFn_B(1:4)//op//dotdmp, form='UNFORMATTED',status='UNKNOWN')
   !----------------------------------------------------------------------------
      nar = 1
   !   ADJUST NEIGB AND NEIGA TO ALLOW FOR NEIG = n-1 FOR TRIPLET-S STATES
      neffb = neigb
      neffa = neiga
      if(lrglb.eq.0.and.nspnb.eq.1) neffb = neigb-1
      if(lrgla.eq.0.and.nspna.eq.1) neffa = neiga-1
      if(leiga) nar = nwa - neiga + 1
      nbr = nwb - neigb + 1
      eiga(1) = wa(2)
      do 12 i=1,nbr
      call inputb
      sumss = 0
      do 13 j=1,nwb
      sumss = sumss + trb(j,kontrol)*db(j)

   13 continue
      pert(i,kontrol) = 1048576*sumss*z**3
   12 continue
      if(hterm(kterm).eq.0) hterm(kterm) = 1.d-40
      write(*,'(A,1PD24.16)') 'Check Ratio =',pert(neffb,kontrol) / (hterm(kterm)*1048576*z**3)
      rewind(8)
      do 90 i=1,nbr
      call inputb
      do 90 j=1,nwa
      sum = 0
      do 91 k=1,nwb
      sum = sum + db(k)*tr(k,j,kontrol)
   91 continue
      temp(i,j) = sum
   90 continue

      rewind(7)
      rewind(8)
      rewind(9)
      do 92 j=1,nar
      if(nar.gt.1) call inputa
      do 92 i=1,nbr
      sum = 0
      do 93 k=1,nwa
      sum = sum + temp(i,k)*da(k)
   93 continue
   92 tr(i,j,kontrol) = 1048576*sum*z**3
      write(9) amm,nar,nbr,waveFn_A(1:7)//'-'//waveFn_B(1:7)
      write(9) (eiga(i), i=1,nar)
      write(9) (eigb(i), i=1,nbr)
      do 94 j=1,nar
   94 write(9) (tr(i,j,kontrol),i=1,nbr)
      write(*,'(A,I1,A,I1,A,D20.12)') '    tr(',neffb,',',neffa, ') =',tr(neffb,neffa,kontrol)
      write(*,*) tr(1,1,kontrol),tr(2,2,kontrol),neffb,neffa
      rewind 9

      close(9)
      close(7)
      close(8)
      write(*,'(A,1PD24.16)') 'Check Ratio =',tr(neffb,neffa,kontrol) / hterm(kterm)/1048576/z**3

      endif
  100 return
  103 write(4,104)
  104 format(' RETURN CALLED DUE TO VANISHING 3-J SYMBOL')
      return
end subroutine breit


subroutine inputb
      
   !? Import Common Blocks
   use a1a_block
   use b1b_block
   use f1f_block

   implicit real*16(a-h,o-z)

   if(nwb.gt.3601)  then 
      write(*,*) 'NWB too big ',nwb,neigb
   end if 

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

   if(db(1).lt.0.)then 
      sign = -1.
   end if 

   db0 = db0*sign
   
   do k=1,nwb
      db(k) = db(k)*sign
   end do

   dbn = db0/db(1)

   return
end subroutine inputb


subroutine diskr(db,n)

   implicit real*16(a-h,o-z)

   dimension db(n)

   read(8) (db(i),i=1,n)

   return
end subroutine diskr


subroutine inputa

   use a1a_block
   use b1b_block
   use f1f_block

   implicit real*16(a-h,o-z)

   if(nwa.gt.3601) then 
      write(*,*) 'NWA too big ',nwa,neiga
   end if 

   call diskra(da(neiga),nwa-neiga+1)

   !? CALCULATION OF SCREENED HYDROGENIC WAVEFUNCTIONS.
   if(iza.ne.0) then
      n = lrgla + neiga
      da0 = da(neiga)
      xx = 1
      l = 0

      do m=1,neiga
         da(m) = da(neiga)*xx*2**(l+lrgla+2)*(iza-1)**(lrgla+1) *sqrt((iza-1)*fac(n+lrgla+1)/(iza*fac(n-lrgla)*fac(2*l+2))) /(fac(2*lrgla+2)*n*(l+1.q0)**(l+2)*(n*z)**(lrgla+1))
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
end subroutine inputa


subroutine diskra(da,n)

   implicit real*16(a-h,o-z)

   dimension da(n)

   read(7) (da(i),i=1,n)

   return
end
