include 'cross1.f90'
include 'spin68.f90'

program p4

    !? VERSION OF P4.FOR FOR S- and P-STATES WITH diVided BASIS SETS. Feb. 16/94.
    !? CONVERTED BACK TO p4.f April 5, 2021
    !? - UPDATED March 29/94.
    !* Converted to Fortran90 July 2023 by Evan Petrimoulx

    !? Import Common Blocks
    use a1a_block
    use b1b_block
    use f1f_block
    use h1h_block
    use maxpow_block

    !? Import Format
    use format 

    !? Import OpenWave
    use wavExt

    implicit real*16 (a-h,o-z)

    dimension hterm(10),htscr(10),nn(50)

    character :: cx, lines(100)
    character(len = 16) :: matout
    character(len = 52) :: title
    character(len = 12) :: nama, namb
    character(len = 51) :: line
    character(len = 12) :: fmt3a(3)
    character(len = 5) :: eigen
    character(len = 4) :: dotdmp

    integer :: endOfFile = 0
    integer :: errorChecker = 0
    
    logical :: lscreen

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
            if(ib1 .gt. ib) then 
                ib = ib1
            end if 

            if(sb(k) .gt. maxcb) then 
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

        read(1,'(1X,A5)',end=28) eigen

        if(eigen.eq.'EIGEN') then
            read(1,fmt3a(if)) (eigb(i),i=1,nwb-neigb+1)
            leigb = .true.
        end if

    28 close(1,status='KEEP')

        if(leigb) then
            !? OPEN .DMP FILE WITH SQUARE ARRAY OF WAVEFUNCTION COEFFICIENTS.
            write(*,'(2A)') 'p^4 OPENING ',waveFn_B(1:7)//dotdmp
            call routedOpen(8,file=waveFn_B(1:7)//dotdmp,form='UNFORMATTED',status='OLD')
        end if

        close(1,status='KEEP')

        write(*,*)

        if(nwa.gt.2401.or.nwb.gt.2401) then
            write(*,*) 'NWA, NWB TOO BIG ',nwa,nwb
            stop
        end if

        call breit(hterm,htscr,ks,xp,amm)

        write(*, fmt42) (hterm(i),i=1,3)

        call routedOpen(1,file=matout,status='OLD')

        errorChecker = 0
        do while(errorChecker .eq. 0)
            read(1,'(A14)',iostat = errorChecker) fnin
        end do

        backspace 1

        write(1,'(A14,1X,A14,1X,A47)')waveFn_A,waveFn_B,title
        write(1, fmt42) (hterm(i),i=1,3)
        write(1, fmt42) (hterm(i),i=4,ks+3)

        close(1,status='KEEP')

        write(*, fmt10) line
        write(*, fmt10)
        
        !? RESET INPUT FILE FOR NEXT CALCULATION.
        cycle

        rewind 5

        do i=1,100
            ip = i-1
            read(5,'(A)',iostat = errorChecker) lines(i)

            if(errorChecker .gt. 0) then
                print*, "An error has occurred, stopping the program."

            else if(errorChecker .lt. 0) then 
                exit
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

end program p4


subroutine breit(hterm,htscr,ks,xp,amm)

      !? Import Common Blocks
      use at2_block
      use a1a_block
      use b1b_block
      use c1c_block
      use d1d_block
      use f1f_block
      use h1h_block
      use maxpow_block

      implicit real*16 (a-h,o-z)

      dimension hterm(10),htscr(10),ipow(12),ar(6),al(6),xc(6,6),dar(6),dal(6),a(6,6),d(6),drs(6),trb(2401),pert(2401),tr(2401,2401),temp(2401,2401)
      integer p,q,s,a1,b1,c1,a4,b4,c4,pam,qam
      logical lgo,lcalc
      character op*3,title*15,dotdmp*4
    !   FOR D STATES
    !      DATA ISUMD,MP2D,MQ2D,MS2D,MSLD/0,-1,1,-2,1/
    !   FOR S STATES
        data isumd,mp2d,mq2d,ms2d,msld/0,2,3,0,2/
    !   POWERS OF R1 AND R12.  POWER OF R2 = 0.
        data ipow/-2,0, -1,0, 0,0, 0,-2, 1,-2 ,0,-1/,xc/36*0.q0/
    !
    !    IST CHOOSES THE COORDINATES OF ELECTRON 1 OR 2 AND 'I' SPECIFIES WHICH
    !    SET OF COEFFICIENTS OF PL(COSTHETA12) IS TO BE USED.
    !
    !
    !      READ(*,*) ISUMD,MP2D,MQ2D,MS2D,MSLD
    !   SET KROSS = 5 OR 9 FOR USE IN 'CROSS'.
        ermac = 1.1q0
    1000 ermac = ermac/2
        test = 1.1q0 + ermac
        if(test.gt.1.1q0) go to 1000
    !   USE .DML TO DENOTE THE FINITE MASS CASE.
        dotdmp = '.DMP'
        if(amm.ne.0) dotdmp = '.DML'
        kross = 5
        klog = -3
        kang = 0
        pam = 0
        qam = 0
        do 1 i=1,nwa
        if(pa(1,i).gt.pam) pam=pa(1,i)
        1 if(pa(2,i).gt.qam) qam=pa(2,i)
        mla = lrgla
        mlb = lrglb
        md = - mla + mlb
        do 6 i=1,10
        hterm(i) = 0.0
        htscr(i) = 0.0
        6 continue
        do 7 j=1,nwb
    !      TOV(J) = 0
        7 trb(j) = 0
        tr = 0
    !
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
        if(nbl1.eq.1.and.nbl3.eq.1.and.z.ne.1.) lgo = .true.
        if(izb.eq.0.and.nbl3.lt.nbl1) go to 20
    !
    !   SUMMATION OVER LEFT HAND STATE DIRECT AND EXCHANGE TERMS.
    !   THE ORDERING IS D-D, E-E, D-E AND E-D.
    !
        do 220 jdx=1,2
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
        l13 = l1 + l3
        call cross(l3,l4,l1,l2,lrglb,mlb,lrgla,mla,ide,jde)
    !
        lam = lt(1) - 1
        lam1 = lam + 1
    !   NOTE THAT LT, LTP AND LS ARE TWICE T, T' AND S ON p. 3 OF NOTES.
        lt =  -(l3*(l3+1) - l1*(l1+1) - lam*(lam+1))
        ltp = -(l1*(l1+1) - l3*(l3+1) - lam*(lam+1))
        ls =  -(l1*(l1+1) + l3*(l3+1) - lam*(lam+1))
        mlt(10) = 1
        lt(10) = 1
        cplt(1,10) = 0.5
        mlt(11) = 3
        lt(11) = 3
        cplt(1,11) = 0
    !
        mlp2 = mlt(2)
        lp2 = lt(2)
        mlp3 = mlt(3)
        lp3 = lt(3)
        do 54 l=1,lp3
    54 cplt(l,16) = cplt(l,3)
        mlp16 = mlp3 - 1
        if(mlp16.lt.1) mlp16 = mlp16 + 2
        lp16 = lp3 - 1
        mlt(16) = mlp16
        lt(16) = lp16
        do 50 ll=mlp16,lp16,2
        l = lp16 - ll + mlp16
    !   FORM COS**2(THETA) - 1 AND REDUCE BY RECURSION RELATION.
        if(l.eq.1.and.cplt(l+1,16).gt.1.d-10) stop 1
        cplt(l,16) = (2*l-1)*cplt(l+1,16)
        if(l.eq.1) go to 50
        cplt(l-1,16) = cplt(l-1,16) + cplt(l+1,16) - cplt(l-1,1)
        if(abs(cplt(l-1,16)).lt.1.d-10) cplt(l-1,16) = 0.
    50 continue
        do 52 ll=mlp16,lp16,2
        l = lp16 - ll + mlp16
        if(cplt(l,16).ne.0.) go to 53
    52 lt(16) = lt(16) - 2
    53 if(kang.eq.0) go to 110
        write(6,991) lss,l3,l4,l1,l2,lrglb,mlb,lrgla,mla,ide,jde,nb3,nb4,nb1,nb2,izb
    991 format(20i3)
        write(6,111)
        do 210 ic=1,16
        if(lt(ic).le.0) go to 210
        ic2 = lt(ic)
        write(6,111) ic,(cplt(ic3,ic), ic3=  1,ic2)
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
        npj = 100
        nqj = 100
        nsj = 100
        do 51 j=nb3,nb4
        if(nsj.gt.sb(j)) nsj = sb(j)
        if(npj.gt.pb(jde,j)) npj = pb(jde,j)
    51 if(nqj.gt.pb(jdep,j)) nqj = pb(jdep,j)
        mp1 = 3
        ms1 = 3
        if(nsj.lt.0) ms1 = 5
        if(npj.lt.l3.or.nqj.lt.l4) mp1 = 5
        msl = lt(1) + msld
        imin = -2 + lrgla + lrglb
        if(nsj.lt.0.or.npj.lt.l3.or.nqj.lt.l4) imin = imin - 2
        lsuma = 2*id
        lsumb = l2*(iabs(l1-l2)/2)
        ltota = lrgla + neiga - 1
        ltotb = lrglb + neigb - 1
    !      IF((NB1.EQ.1.OR.NB3.EQ.1).AND.NEIGA.NE.0) GO TO 78
        isum = ijkb(nbl3) + ijka(nbl1) + isumd
        mp2 = iib(nbl3) + (1-id)*iia(nbl1) + id*jja(nbl1) + mp1 + mp2d+ 3*(1-id)      
        mq2 = jjb(nbl3) + (1-id)*jja(nbl1) + id*iia(nbl1)+ mp1+mq2d - 1 + id
        ms2 = ms1 + kkb(nbl3) + kka(nbl1) + 2*lt(1) + ms2d
        go to 80
    78 if(nb3.eq.1.and.nb1.eq.1) go to 179
    !   LIMITS FOR HYDROGENIC '1,N' TERMS.
        if(nb3.eq.1) go to 79
        isum = ltota + ijkb(nbl3) + isumd
        mp2 = iib(nbl3) + id*ltota + mp1 + mp2d + 1 - id 
        mq2 = jjb(nbl3) + (1-id)*ltota + mp1 + mq2d 
        ms2 = ms1 + kkb(nbl3) + 2*lt(1) + ms2d
        go to 80
    79 isum = ltotb + ijka(nbl1) + isumd
    !   LIMITS FOR HYDROGENIC 'N,1' TERMS.
        mp2 = (1-id)*iia(nbl1) + id*(jja(nbl1)+ltotb) + mp1 + mp2d+1-id
        mq2 = (1-id)*(jja(nbl1)+ltotb) + id*iia(nbl1) + mp1 + mq2d+2*id
        ms2 = ms1 + kka(nbl1) + 2*lt(1) + ms2d
        go to 80
    179 isum = ltotb + ltota + isumd
        mp2 = (1-id)*ltota + id*(ltotb+ltota) + mp1 + mp2d
        mq2 = (1-id)*(ltotb+ltota) + id*ltota + mp1 + mq2d
        ms2 = ms1 + 2*lt(1) + ms2d
    80 continue
        lmax1 = min0(l1+l3,l2+l4) + 2
        lmax2 = min0(l2+l3,l1+l4) + 2
        lmax = lmax1
        if((nbx4.gt.nbx3.or.nbx2.gt.nbx1))lmax = max0(lmax1,lmax2)
        ms2 = ms2 + 2*(lmax - lt(1) - 1) - 2
        msl = msl + lmax - lt(1) - 1
        if(nbl3.eq.nbl1) write(*,180) mp1,mp2,mq2,ms1,ms2,isum,y11,y22,l3,l4,l1,l2,nbl1,nbl3
    180 format(' GENINT ',6i3,2f10.6,4i2,3x,2i2,$)
        if(kex.eq.1) then
            mtemp = mp2
            mp2 = mq2
            mq2 = mtemp
        endif
        if(id.eq.1) then
        mp2 = mp2 + 1
        mq2 = mq2 + 1
        endif
    !      ISUM = ISUM + 1
        call genint(id,nbl3,nbl1,lmax,fac,2.q0)
    !      WRITE(*,'(5H DONE)')
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
        ar(1)= a1*(a1+1) - l1*(l1+1)
        ar(2)= -2.*y1*(a1+1)
        ar(3) = y1*y1
        ar(4) = c1*(c1+1) + 2*a1*c1
        ar(5) = - 2.*y1*c1
        ar(6) = - 2.*c1
    !
    !   SUMMATION OVER LEFT HAND STATE BASIS FUNCTIONS.
    !
    224 do 21 j=nb3,nb4
        if(izb.eq.0.and.j.lt.i) go to 21
        psi2 = da(i)*db(j)*sgn
        if(izb.eq.0.and.j.ne.i) psi2 = psi2*2.
        psia = psi2/db(j)*1048576
        a4 = pb(jde,j)
        b4 = pb(jdep,j)
        c4 = sb(j)
        y3 = yb(jde,j)
        y4 = yb(jdep,j)
        p = pa(ide,i) + pb(jde,j)
        q = pa(idep,i) + pb(jdep,j)
        s = sa(i) + sb(j)
        y_ = y1 + y3
        yp = y2 + y4
    !
    !   <p**4> VERSION.
    !
    250 sum = 0.
        ov = spl(p,q,s,1)
        al(1)= a4*(a4+1) - l3*(l3+1)
        al(2)= -2.*y3*(a4+1)
        al(3) = y3*y3
        al(4) = c4*(c4+1) + 2*a4*c4
        al(5) = - 2.*y3*c4
        al(6) = - 2.*c4
        do 141 itr=1,6
        do 140 jtr=itr,6
    140 a(itr,jtr) = al(itr)*ar(jtr) + al(jtr)*ar(itr)
        a(itr,itr) = al(itr)*ar(itr)
        drs(itr) = 0.5*(al(itr)*ar(6)*(2*a1*lam1 - lt)+ al(6)*ar(itr)*(2*a4*lam1 - ltp))
    141 d(itr) = al(itr)*ar(6)*a1 + al(6)*ar(itr)*a4
        rs = 0.
        if(s.ne.0) rs = 1.q0/s
        if(a(1,1).ne.0.) sum = a(1,1)*spl(p-4,q,s,1)
        if(c.ne.0.) sum = sum + c*spl(p-3,q,s,1)
        sum = sum + (a(1,3) + a(2,2) + y_*rs*a(1,5) - (p+1)*rs*a(2,5)) *spl(p-2,q,s,1) + (a(2,3) + y_*rs*a(2,5) - (p+2)*rs*a(3,5)) *spl(p-1,q,s,1) + (a(3,3) - 2.5 + y_*rs*a(3,5))*ov
        if(c1.eq.0.and.c4.eq.0) go to 142
        if(s.eq.0) then
            sum = sum - p*a(1,5) *splg(p-3,q,0,1,0) + (y_*a(1,5) - (p+1)*a(2,5))*splg(p-2,q,0,1,0) + (y_*a(2,5) - (p+2)*a(3,5))*splg(p-1,q,0,1,0) +  y_*a(3,5) *splg(p,q,0,1,0)
        endif
    !
        is = s - 2
        if(a(1,4).ne.0) sum = sum + a(1,4)*spl(p-2,q,is,1)
        rs = 0.
        if(is.ne.0) rs = 1.q0/is
        sum = sum + (a(2,4) - (p+2)*rs*a(4,5))*spl(p-1,q,is,1) + (a(3,4) + y_*rs*a(4,5) + a(5,5))*spl(p,q,is,1)
        if(is.eq.0) sum = sum - (p+2)*a(4,5)*splg(p-1,q,is,1,0) + y_*a(4,5)*splg(p,q,is,1,0)
        if(is.eq.0) rs = 1.
        sum = sum + rs*a(5,5)*splg(p-1,q+1,is,16,is)
        if(a(4,4).ne.0) sum = sum + a(4,4)*spl(p,q,s-4,1)
    132 format(' TRAP ',i2,5i4,2d12.5)
    !
        if(lam.eq.0) then
            if(d(1).ne.0.) sum = sum + d(1)*spl(p-3,q+1,is,2)
            if(d(2).ne.0.) sum = sum + d(2)*spl(p-2,q+1,is,2)
                            sum = sum + d(3)*spl(p-1,q+1,is,2)
            if(d(4).ne.0.) sum = sum + d(4)*spl(p-1,q+1,s-4,2)
            if(d(5).ne.0.) sum = sum + rs*d(5)*(-(p+1)*splg(p-2,q+1,is,2,is) + y_*splg(p-1,q+1,is,2,is))
            sum = sum + a(6,6)*(2*rs*(a1*a4+ls/4)*splg(p-3,q+1,is,2,is) + a1*a4*spl(p-2,q+2,s-4,1))
        else if(lam.eq.1) then
            rs0 = 1
            if(s.ne.0) rs0 = 1.q0/s
            if(d(1).ne.0.) sum = sum + d(1)*spl(p-3,q+1,is,10)
            if(drs(1).ne.0.) sum = sum + drs(1)*rs0*splg(p-4,q,s,1,s)
            if(d(2).ne.0.) sum = sum + d(2)*spl(p-2,q+1,is,10)
            if(drs(2).ne.0.) sum = sum + drs(2)*rs0*splg(p-3,q,s,1,s)
            sum = sum + d(3)*spl(p-1,q+1,is,10) + drs(3)*rs0*splg(p-2,q,s,1,s)
            if(d(4).ne.0.) sum = sum + d(4)*spl(p-1,q+1,s-4,10)
            if(drs(4).ne.0.) sum = sum + drs(4)*rs*splg(p-2,q,s-2,1,s-2)
            if(d(5).ne.0.) sum = sum + rs*d(5)* (-(p+1)*splg(p-2,q+1,is,10,is) + y_*splg(p-1,q+1,is,10,is))
            if(drs(5).ne.0.) sum = sum + 0.5*drs(5)*(rs*(splg(p-1,q,is,1,is) - splg(p-3,q+2,is,1,is)) + rs0*splg(p-3,q,s,1,s))
                cplt(3,11) = (4*a1-lt)*(4*a4-ltp)/16.q0
                cd6 = ((4*a1-lt)*(4*a4-ltp) + 2*(2*a1+lt)*(2*a4+ltp))/24.q0
                sum = sum + a(6,6)*(cd6*spl(p-2,q+2,s-4,1) + rs*splg(p-3,q+1,is,11,is))
        else
            write(*,*) 'LAM TOO BIG FOR THIS VERSION ',lam
            stop
        endif
    142 hterm(1) = hterm(1) + sum*psi2
        trb(j) = trb(j) + sum*psi2/db(j)
        tr(j,i) = tr(j,i) + sum*psi2/(da(i)*db(j))
        if(abs(sum).gt.1.d70) write(*,132)6,i,j,p,q,s,sum
        hterm(2) = hterm(2) + ov*psi2
    21 continue
    220 continue
    20 continue
    620 continue
        if(leigb) then
            op = 'p^4'
            write(*,'(2A)') 'p^4 OPENING ', waveFn_A(1:4)//waveFn_B(1:4)//op//dotdmp
            open(9,file=waveFn_A(1:4)//waveFn_B(1:4)//op//dotdmp, form='UNFORMATTED',status='UNKNOWN')

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
        sumss = sumss + trb(j)*db(j)
    13 continue
        pert(i) = 1048576*sumss*z**4
    12 continue
        write(*,'(A,1PD24.16)') 'Check Ratio =',pert(neffb)/(hterm(1)*1048576*z**4)
        rewind(8)
        do 90 i=1,nbr
        call inputb
        do 90 j=1,nwa
        sum = 0
        do 91 k=1,nwb
    91 sum = sum + db(k)*tr(k,j)
        temp(i,j) = sum
    90 continue
    !
        rewind(7)
        rewind(8)
        rewind(9)
        do 92 j=1,nar
        if(nar.gt.1) call inputa
        do 92 i=1,nbr
        sum = 0
        do 93 k=1,nwa
    93 sum = sum + temp(i,k)*da(k)
    92 tr(i,j) = 1048576*sum*z**4
        write(9) amm,nar,nbr,waveFn_A(1:7)//'-'//waveFn_B(1:7)
        write(9) (eiga(i), i=1,nar)
        write(9) (eigb(i), i=1,nbr)
        do 94 j=1,nar
    94 write(9) (tr(i,j),i=1,nbr)
        write(*,'(A,I1,A,I1,A,D20.12)') '    tr(',neffb,',',neffa,') =',tr(neffb,neffa)
        rewind 9
        close(9)
        close(7)
        close(8)
        write(*,'(A,1PD24.16)') 'Check Ratio =',tr(neffb,neffa) /hterm(1)/1048576/z**4
        endif
    100 return
    end subroutine breit


    subroutine inputb

        !? Import COmmon Blocks
        use a1a_block
        use b1b_block

        implicit real*16(a-h,o-z)

        call diskr(db(neigb),nwb-neigb+1)

        !? CALCULATION OF SCREENED HYDROGENIC WAVEFUNCTIONS.
        if(izb.ne.0) then 
                n = lrglb + neigb
                db0 = db(neigb)
                xx = 1

                l = 0

                do m=1,neigb
                    db(m) = db(neigb)*xx*2**(l+lrglb+2)*(izb-1)**(lrglb+1) *sqrt((izb-1)*fac(n+lrglb+1)/(izb*fac(n-lrglb)*fac(2*l+2))) /(fac(2*lrglb+2)*n*(l+1.q0)**(l+2)*(n*z)**(lrglb+1))
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

        !? CALCULATION OF SCREENED HYDROGENIC WAVEFUNCTIONS.
        
        if(iza.ne.0) then
                n = lrgla + neiga
                da0 = da(neiga)
                xx = 1
                l = 0

                do m=1,neiga
                    da(m) = da(neiga)*xx*2**(l+lrgla+2)*(iza-1)**(lrgla+1) *sqrt((iza-1)*fac(n+lrgla+1)/(iza*fac(n-lrgla)*fac(2*l+2)))/(fac(2*lrgla+2)*n*(l+1.q0)**(l+2)*(n*z)**(lrgla+1))
                    xx = -xx*2*(iza-1)*(neiga-m)/((2*lrgla+m+1)*m*n*iza)
                end do
        end if

        sign = 1.

        if(da(1).lt.0.)then
                sign = -1.
        end if

        da0 = da0*sign

        do k=1,nwa
                da(k) = da(k)*sign
        end do

        dan = da0/da(1)

        return
    end subroutine inputa
    !
    !
    subroutine diskra(da,n)

        implicit real*16(a-h,o-z)

        dimension da(n)

        read(7) (da(i),i=1,n)

        return
end subroutine diskra