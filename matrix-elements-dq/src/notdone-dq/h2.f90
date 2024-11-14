include 'cross1.f90'
include 'spin68.f90'

program h2

   !? Import Common Blocks
   use a1a_block
   use b1b_block
   use h1h_block
   use f1f_block
   use maxpow_block

   !? Import RoutedOpen Subroutine 
   use wavExt

   !? Import Formatting
   use format

   implicit real*16 (a-h,o-z)

   !? Initializations
   dimension hterm(10),htscr(10),nn(50)

   character(len = 57) :: lines(100)
   character cx*1,eigen*5,dotdmp*4
   character*16 matout
   character title*52,nama*12,namb*12,line*51,fmt3a(3)*12

   integer :: errorChecker = 0
   integer :: endOfFile = 0

   logical lscreen

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

      read(1,'(1X,A5)',iostat = errorChecker) eigen

      if(errorChecker .eq. 0) then
         if(eigen.eq.'EIGEN') then
            read(1, fmt3a(if)) (eigb(i),i=1,nwb-neigb+1)
            leigb = .true.
         end if 
      end if 

      close(1,status='KEEP')
      
      ! open .dmp file with square array of wavefunction coefficients.
      if(leigb) then
         write(*,'(2A)') 'H2 OPENING ',waveFn_B(1:7)//dotdmp
         call routedOpen(8,file=waveFn_B(1:7)//dotdmp,form='UNFORMATTED',status='OLD')
      endif

      close(1,status='KEEP')

      call breit(hterm,htscr,amm)

      do i=1,3
         hterm(i) = hterm(i)*1048576
         htscr(i) = htscr(i)*1048576
         htscr(i) = htscr(i)*0.5/(da0*db0)
      end do 

      write(*, fmt1) (i,htscr(i),hterm(i),i=1,3,3)

      open(1,file=matout,status='OLD')


      do while(errorChecker .eq. 0)
         read(1,'(A14)', iostat = errorChecker) fnin
      end do 

      backspace 1

      write(1,'(A14,1X,A14,1X,A47)')waveFn_A,waveFn_B,title
      write(1, fmt2) (hterm(i),i=1,3,3)

      close(1,status='KEEP')

      write(*, fmt10) lines(1)
      write(*, fmt10)      

      ! reset input file for next calculation.
      cycle
      rewind 5

      do i=1,100
         ip = i-1
         read(5,'(A)',iostat = errorChecker) lines(i)

         if(errorChecker .ne. 0) then 
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
end program h2


subroutine breit(hterm,htscr,amm)

      !? Import Common Blocks
      use a1a_block
      use b1b_block
      use c1c_block
      use d1d_block
      use f1f_block
      use h1h_block
      use maxpow_block

      implicit real*16 (a-h,o-z)

      dimension hterm(10),htscr(10),ipow(30),trb(3601),pert(3601),tr(3601,3601),temp(3601,3601)
      integer p,q,s,a1,b1,c1,a4,b4,c4
      logical lgo,lcalc
      character op*3,title*15,dotdmp*4

      !for p states . . .
      !     data isumd,mp2d,mq2d,ms2d,msld/-1,3,3,-1,1/
      !for d states . . .
      data isumd,mp2d,mq2d,ms2d,msld/-1,4,4,-1,2/
      data ipow/-1,-1,-1, -1,0,-1, 0,-1,-1, 0,1,-3, 1,0,-3, -1,1,-3, 1,-1,-3, 2,-1,-3, -1,2,-3, 1,1,-3/
   !
   !   ist chooses the coordinates of electron 1 or 2 and 'I' specifies which
   !   set of coefficients of pl(costheta12) is to be used.
   !
      ermac = 1.1q0
 1000 ermac = ermac/2
      test = 1.1q0 + ermac
      if(test.gt.1.1q0) go to 1000
!  use .dml to denote the finite mass case.
      dotdmp = '.DMP'
      if(amm.ne.0) dotdmp = '.DML'
      kross = 5
      klog = -2
      kang = 1
      kang = 0
      sq3 = sqrt(3.q0)
      mla = lrgla
      mlb = lrglb
      md = - mla + mlb
      do 6 i=1,10
      hterm(i) = 0.0
      htscr(i) = 0.0
    6 continue
      tr = 0
      trb = 0
!
!  the outer loop corresponds to the right hand state and is denoted by
!  'A' or'1' or '2' in various places.
!
!  summation over right hand state direct and exchange terms.
      do 620 idx=1,2
!  summation over blocks of terms with the same non-linear paremeters.
      do 20 nhi=1,nbxa
      nbx1 = nblxa(nhi) + 1
      nbx2 = nblxa(nhi+1)
      do 20 nhj=1,nbxb
      nbx3 = nblxb(nhj) + 1
      nbx4 = nblxb(nhj+1)
!
!  the outer loop corresponds to the right hand state and is denoted by
!  'A' or'1' or '2' in various places.
!
      do 20 nbl1=nbx1,nbx2
      nb1 = nblka(nbl1) + 1
      nb2 = nblka(nbl1+1)
!
!  the inner loop corresponds to the left hand state and is denoted by
!  'B' or '3' or '4' in various places (3 means electron 1 and 4 means electron
!   2).
!
      do 20 nbl3=nbx3,nbx4
      nb3 = nblkb(nbl3) + 1
      nb4 = nblkb(nbl3 + 1)
      lgo = .false.
      if(nbl1.eq.1.and.nbl3.eq.1) lgo = .true.
      if(izb.eq.0.and.nbl3.lt.nbl1) go to 20
!
!  summation over left hand state direct and exchange terms.
!  the ordering is d-d, e-e, d-e and e-d.
!
!    do 220 jdx=1,2
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
      ll3 = l3*(l3+1)
      ll4 = l4*(l4+1)
      call cross(l3,l4,l1,l2,lrglb,mlb,lrgla,mla,ide,jde)
      mlp2 = mlt(2)
      mlp3 = mlt(3)
      lp2 = lt(2)
      lp3 = lt(3)
      do 120 ll=mlp2,lp2,2
      cplt(ll,4) = cplt(ll,9)
  120 cplt(ll,5) = cplt(ll,10)
      do 121 n=4,15
      mlt(n) = mlp2
  121 lt(n) = lp2
      mlp3 = mlt(3) - 1
      if(mlp3.lt.1) mlp3 = mlp3 + 2
      lp3 = lt(3) - 1
      mlt(3) = mlp3
      lt(3) = lp3
      do 50 ll=mlp3,lp3,2
      l = lp3 - ll + mlp3
!  form cos**2(theta) - 1 and reduce by recursion relation.
      if(l.eq.1.and.cplt(l+1,3).gt.1.d-10) stop
      cplt(l,3) = (2*l-1)*cplt(l+1,3)
      if(l.eq.1) go to 50
      cplt(l-1,3) = cplt(l-1,3) + cplt(l+1,3) - cplt(l-1,1)
      if(abs(cplt(l-1,3)).lt.1.d-10) cplt(l-1,3) = 0.
   50 continue
      do 52 ll=mlp3,lp3,2
      l = lp3 - ll + mlp3
      if(cplt(l,3).ne.0.) go to 53
   52 lt(3) = lt(3) - 2
   53 mlt(16) = mlp3
      lt(16) = lt(3)
      if(kang.eq.0) go to 110
      write(4,991) lss,l3,l4,l1,l2,lrglb,mlb,lrgla,mla,ide,jde,nb3,nb4,nb1,nb2,izb
  991 format(20i3)
      write(4,111)
      do 210 ic=1,16
      if(lt(ic).le.0) go to 210
      ic2 = lt(ic)
      write(4,111) ic,(cplt(ic3,ic), ic3=  1,ic2)
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
      lmax = lt(2)
      mp1 = 2
      ms1 = 4
      msl = lt(1) + msld
      imin = -3 + lrgla + lrglb
      lsuma = 2*id
      lsumb = l2*(iabs(l1-l2)/2)
      ltota = lrgla + neiga - 1
      ltotb = lrglb + neigb - 1
      if((nb1.eq.1.or.nb3.eq.1).and.z.gt.1) go to 78
      isum = ijkb(nbl3) + ijka(nbl1) + isumd
      mp2 = iib(nbl3) + (1-id)*iia(nbl1) + id*jja(nbl1) + mp1 + mp2d + 1 - id
      mq2 = jjb(nbl3) + (1-id)*jja(nbl1) + id*iia(nbl1)+ mp1+mq2d - 1 + id + id
      ms2 = ms1 + kkb(nbl3) + kka(nbl1) + 2*lt(1) + ms2d + 1
      go to 80
   78 if(nb3.eq.1.and.nb1.eq.1) go to 179
      if(nb3.eq.1) go to 79
      isum = ltota + ijkb(nbl3) + isumd
      mp2 = iib(nbl3) + id*ltota + mp1 + mp2d + 1 - id
      mq2 = jjb(nbl3) + (1-id)*ltota + mp1 + mq2d - 1 + id 
      ms2 = ms1 + kkb(nbl3) + 2*lt(1) + ms2d + 1
      go to 80
   79 isum = ltotb + ijka(nbl1) + isumd
      mp2 = (1-id)*iia(nbl1) + id*(jja(nbl1)+ltotb) + mp1 + mp2d+1-id
      mq2 = (1-id)*(jja(nbl1)+ltotb) + id*iia(nbl1) + mp1 + mq2d-1+id+2*id
      ms2 = ms1 + kka(nbl1) + 2*lt(1) + ms2d + id
      go to 80
  179 isum = ltotb + ltota + isumd
      mp2 = (1-id)*ltota + id*(ltotb+ltota) + mp1 + mp2d + 1 - id
      mq2 = (1-id)*(ltotb+ltota) + id*ltota + mp1 + mq2d - 1 + id
      ms2 = ms1 + 2*lt(1) + ms2d
   80 continue
      if(id.eq.0.and.(nbx4.gt.nbx3.or.nbx2.gt.nbx1)) ms2 = ms2 +2
      if(id.eq.0.and.(nbx4.gt.nbx3.and.nbx2.gt.nbx1)) ms2 = ms2 +2
      if(id.eq.1.and.nbx2.gt.nbx1.and.l3.gt.0) ms2 = ms2 +2
      if(id.eq.0.and.(nbx4.gt.nbx3.or.nbx2.gt.nbx1)) msl = msl + 1
      if((nbx4.gt.nbx3.or.nbx2.gt.nbx1))lmax =min0(lrgla +2,lmax+1)
      if((nbx4.gt.nbx3.and.nbx2.gt.nbx1))lmax =min0(lrgla +2,lmax+1)
      if(msl.gt.12) stop 180
      if(kex.eq.1) then
        mtemp = mp2
        mp2 = mq2
        mq2 = mtemp
      endif
      if(nbl3.eq.nbl1) then
         write(*,180) kex,mp1,mp2,mq2,ms1,ms2,isum,y11,y22,l3,l4,l1,l2,nbl1,nbl3
      end if 
  180 format(' GENINT',i1,6i3,2f10.6,4i2,3x,2i2,$)
      call genint(id,nbl3,nbl1,lmax,fac,0.q0)
   77 delt12 = 0.0
      jdxx = jdx
      if(kex.eq.1) jdxx = 3 - jdx
!
!  summation over right hand state basis functions.
!
      do 21 i=nb1,nb2
      a1 = pa(ide,i)
      b1 = pa(idep,i)
      c1 = sa(i)
      y1 = ya(ide,i)
      y2 = ya(idep,i)
      ar1 = a1*b1 + 2*c1*(a1+b1+c1) + 0.5q0*(lrgla*(lrgla+1)-l1*(l1+1)-l2*(l2+1))
      ar2 = y2*(a1+2*c1)
      ar3 = y1*(b1+2*c1)
      ar4 = 2.*y1*y2
      do 122 ll=mlp2,lp2,2
      cplt(ll,7) = 1.5q0*(a1*b1*cplt(ll,2) + a1*cplt(ll,4) + b1*cplt(ll,5) + cplt(ll,6))
      cplt(ll,8) = a1*cplt(ll,2) + cplt(ll,5)
      cplt(ll,9) = b1*cplt(ll,2) + cplt(ll,4)
      cplt(ll,10) = y1*(-2.*c1*cplt(ll,2) - 0.5*b1*cplt(ll,2) + 0.5*cplt(ll,4))
      cplt(ll,11) = y2*(-2.*c1*cplt(ll,2) - 0.5*a1*cplt(ll,2) + 0.5*cplt(ll,5))
      cplt(ll,12) = 2.*c1*cplt(ll,8) + 0.5*(a1*b1*cplt(ll,2) - a1*cplt(ll,4) + b1*cplt(ll,5) - cplt(ll,6))
      cplt(ll,13) = 2.*c1*cplt(ll,9) + 0.5*(a1*b1*cplt(ll,2) + a1*cplt(ll,4) - b1*cplt(ll,5) - cplt(ll,6))
      cplt(ll,14) = -0.5*y1*cplt(ll,9)
      cplt(ll,15) = -0.5*y2*cplt(ll,8)
      cplt(ll,8) = -1.5q0*y2*cplt(ll,8)
      cplt(ll,9) = -1.5q0*y1*cplt(ll,9)
  122 continue
      do 125 ll=mlp3,lp3,2
  125 cplt(ll,16) = y1*y2*cplt(ll,3)
      do 123 n=7,15
      lt(n) = lp2
      do 123 ll=mlp2,lp2,2
      lpn = lt(n)
      if(abs(cplt(lpn,n)).lt.1.d-06) lt(n) = lt(n) - 2
  123 continue
!
!  summation over left hand state basis functions.
!
      do 221 j=nb3,nb4
      if(izb.eq.0.and.j.lt.i) go to 221
      psi2 = da(i)*db(j)*sgn*2.
      if(izb.eq.0.and.j.ne.i) psi2 = psi2*2.
      y3 = yb(jde,j)
      y4 = yb(jdep,j)
      p = pa(ide,i) + pb(jde,j)
      q = pa(idep,i) + pb(jdep,j)
      s = sa(i) + sb(j)
      if(s+l1+l3.eq.0) go to 221
      y = y1 + y3
      yp = y2 + y4
!
      kkr = -2
      sum = 0.
      do 124 itr=7,15
      kkr = kkr + 3
      if(lt(itr).le.0) go to 124
      ip = p + ipow(kkr)
      iq = q + ipow(kkr+1)
      is = s + ipow(kkr+2)
      sum = sum + spl(ip,iq,is,itr)
      if(abs(sum).lt.1.d70) go to 124
!     write(*,131) 1,i,j,itr,p,q,s,ip,iq,is,sum
  131 format(' TRAP',i2, 9i4,2d12.5)
  124 continue
      is = s - 1
      cc = 1.
      if(is.ne.0) cc = cc/is
      sum = sum + cc*splg(p,q,is,16,is)
      if(ar1.ne.0.) sum = sum - ar1*spl(p,q,s-3,1)
      if(ar2.ne.0.) sum = sum + ar2*spl(p,q+1,s-3,1)
      if(ar3.ne.0.) sum = sum + ar3*spl(p+1,q,s-3,1)
      sum = sum + ar4*spl(p,q,s-1,2)
!     if(abs(sum).gt.1.d70) write(*,132) i,j,p,q,s,sum
  132 format(5i4,d20.12,i4)
      hterm(1) = hterm(1) + sum*psi2
      trb(j) = trb(j) + sum*psi2/db(j)
      tr(j,i) = tr(j,i) + sum*psi2/(da(i)*db(j))
      if(lgo) htscr(1) = htscr(1) + sum*psi2
   28 continue
  221 continue
   21 continue
  220 continue
   20 continue
  620 continue
      if(leigb) then
         op = 'H_2'
         write(*,'(2A)') 'H_2 OPENING ', waveFn_A(1:4)//waveFn_B(1:4)//op//dotdmp
         open(9,file=waveFn_A(1:4)//waveFn_B(1:4)//op//dotdmp, form='UNFORMATTED',status='UNKNOWN')
         nar = 1
   !  adjust neigb and neiga to allow for neig = n-1 for triplet-s states
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
      sumss = sumss + trb(j)*db(j)
   13 continue
      pert(i) = 1048576*sumss*z**3
   12 continue
      write(*,'(A,1PD24.16)') 'Check Ratio =',pert(neffb)/(hterm(1) *1048576*z**3)
      rewind(8)
      do 90 i=1,nbr
      call inputb
      do 90 j=1,nwa
      sum = 0
      do 91 k=1,nwb
   91 sum = sum + db(k)*tr(k,j)
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
   93 sum = sum + temp(i,k)*da(k)
   92 tr(i,j) = 1048576*sum*z**3
      write(9) amm,nar,nbr,waveFn_A(1:7)//'-'//waveFn_B(1:7)
      write(9) (eiga(i), i=1,nar)
      write(9) (eigb(i), i=1,nbr)
      do 94 j=1,nar
   94 write(9) (tr(i,j),i=1,nbr)
      write(*,'(A,I1,A,I1,A,D20.12)') '    tr(',neffb,',',neffa, ') =',tr(neffb,neffa)
      rewind 9
      close(9)
      close(7)
      close(8)
      write(*,'(A,1PD24.16)') 'Check Ratio =',tr(neffb,neffa)/hterm(1)/1048576/z**3
      endif
  100 return
end subroutine breit


subroutine inputb

      !? Import Common Blocks
      use a1a_block
      use b1b_block

      implicit real*16(a-h,o-z)

      call diskr(db(neigb),nwb-neigb+1)

   !  calculation of screened hydrogenic wavefunctions.
   18 if(izb.eq.0) go to 19
      n = lrglb + neigb
      db0 = db(neigb)
      xx = 1
      l = 0
      do 612 m=1,neigb
      db(m) = db(neigb)*xx*2**(l+lrglb+2)*(izb-1)**(lrglb+1)*sqrt((izb-1)*fac(n+lrglb+1)/(izb*fac(n-lrglb)*fac(2*l+2)))/(fac(2*lrglb+2)*n*(l+1.q0)**(l+2)*(n*z)**(lrglb+1))
  612 xx = -xx*2*(izb-1)*(neigb-m)/((2*lrglb+m+1)*m*n*izb)
   19 continue
      sign = 1.
      if(db(1).lt.0.) sign = -1.
      db0 = db0*sign
      do 15 k=1,nwb
      db(k) = db(k)*sign
   15 continue
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

   !? Import Common Blocks
   use a1a_block
   use b1b_block

   implicit real*16(a-h,o-z)

   call diskra(da(neiga),nwa-neiga+1)
   !  calculation of screened hydrogenic wavefunctions.

   18 if(iza.eq.0) go to 19
      n = lrgla + neiga
      da0 = da(neiga)
      xx = 1
      l = 0
      do 612 m=1,neiga
      da(m) = da(neiga)*xx*2**(l+lrgla+2)*(iza-1)**(lrgla+1) *sqrt((iza-1)*fac(n+lrgla+1)/(iza*fac(n-lrgla)*fac(2*l+2))) /(fac(2*lrgla+2)*n*(l+1.q0)**(l+2)*(n*z)**(lrgla+1))
  612 xx = -xx*2*(iza-1)*(neiga-m)/((2*lrgla+m+1)*m*n*iza)
   19 continue
      sign = 1.
      if(da(1).lt.0.) sign = -1.
      da0 = da0*sign
      do 15 k=1,nwa
      da(k) = da(k)*sign
   15 continue
      dan = da0/da(1)
      return
end subroutine inputa


subroutine diskra(da,n)
   implicit real*16(a-h,o-z)

   dimension da(n)

   read(7) (da(i),i=1,n)

   return
end
