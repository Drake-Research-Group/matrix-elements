!?   PROGRAM TO CALCULATE REDUCED MATRIX ELEMENTS <psiB || T(k) || psiA>
!?   AS DEFINED BY EDMONDS, FOR AN IRREDUCIBLE TENSOR OPERATOR T(k) OF RANK k.
!?   WRITTEN BY G.W.F. DRAKE, DEPARTMENT OF PHYSICS, UNIVERSITY OF WINDSOR,
!?              WINDSOR, ONTARIO, CANADA N9B 3P4.
!?
!?   Updated to Fortran90 and optimized by Evan Petrimoulx July 2023

include 'cross3.f90'
include 'spin68.f90'

program hypr

	!? Include Common Blocks
	use at2_block
	use a1a_block
	use b1b_block
	use h1h_block
	use f1f_block
	use maxpow_block

	!? Include Formatting
	use format 

	!? Include routedOpen routine 
	use wavExt

	implicit real*16 (a-h,o-z)

	dimension :: hterm(10),htscr(10),nn(50)
	dimension :: fab(2)

	character :: lines(100)
	character(len = 16) :: matout
	character(len = 52) :: title
	character(len = 12) :: nama, namb
	character(len = 51) :: line
	character(len = 12) :: fmt3a(3)
	character(len = 5) :: eigen
	character(len = 4) :: dotdmp

	logical :: lscreen
	logical :: exitCode = .false.

	integer :: errorChecker = 0
	integer :: endOfFile = 0

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

			!? find the right format
			fmt3a(1) = '(1X,3D26.19)'
			fmt3a(2) = '(1X,3D30.23)'
			fmt3a(3) = '(1X,3D38.31)'

			do iff=1,3
				if = iff
				if(if.gt.1) backspace(1)
				read(1,fmt3a(if), iostat = errorChecker) wb(1), wb(2), ammb

				if(errorChecker .ne. 0) then
					if(iff .eq. 3) then 
						write(*,*) 'B Format match not found'
						stop
					end if 
					cycle 
				end if

				exit
			end do

			if(amm .ne. ammb) then
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
		

		!? MATEL CALCULATES A COMPLETE REDUCED MATRIX ELEMENT, AS DEFINED BY EDMONDS.
		!? THE RESULT IS RETURNED IN HTERM WITH HTSCR BEING THE SCREENED HYDROGENIC
		!? PART. KT = MULTIPOLE MOMENT OF TRANSITION OPERATOR TK(R1) + TK(R2).
		!? KT = 1

		kt = 2
		ncut = 2

		close(1,status='KEEP')

		call matel(hterm,htscr,kt,ncut)

		do i=1,4
			hterm(i) = hterm(i)*1048576
			htscr(i) = htscr(i)*1048576*0.5/(da0*db0)
		end do 

		write(*, fmt1) (i,htscr(i),hterm(i),i=1,1)

		!? ADD THE RESULT OF THE PRESENT CALCULATION TO THE RESULTS OF PREVIOUS
		!? CALCULATIONS IN THE CUMULATIVE FILE MATOUT.

		call routedOpen(1,file=matout,status='OLD')

		do while(errorChecker .eq. 0)
			read(1,'(A14)',iostat = errorChecker) fnwv1,fnwv2

			if(errorChecker .gt. 0) then 
				print*, "Error Occurred, stopping the program"
				stop 
			end if 
		end do 

		backspace 1
		write(1,'(A14,1X,A14,1X,A47)')waveFn_A,waveFn_B,title
		write(1, fmt40) (hterm(i),i=1,1)
		write(*, fmt10) line
		write(*, fmt10)

		close(1,status='KEEP')

		!? RESET INPUT FILE FOR NEXT CALCULATION.
		exitCode = .true.

		if(.not. exitCode) then 
			rewind 5

			do i=1,100
				ip = i-1
				read(5,'(A)',iostat = errorChecker) lines(i)
				if(errorChecker .gt. 0) then 
					print*, "An error has occurred, the program will now stop."
					stop 
				else if(errorChecker .lt. 0) then 
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
		end if 
    end do 
end program hypr 


subroutine matel(hterm,htscr,kt,ncut)

	!? Import Common Blocks
	use a1a_block
	use b1b_block
	use c1c_block
	use d1d_block
	use f1f_block
	use h1h_block
	use maxpow_block

	!? Import Formatting
	use format

	implicit real*16 (a-h,o-z)

	dimension hterm(30),htscr(30),delt12(6)

	integer p,q,s,a1,b1,c1,a4,b4,c4,pam,qam
	logical lgo,lcalc
	logical :: escapeCode = .false.

	data isumd,mp2d,mq2d,ms2d,msld/-3,0,1,-1,-100/

	ermac = 1.1q0

	do while(test .gt. 1.1q0)
		ermac = ermac/2
		test = 1.1q0 + ermac
	end do

	kross = 9
	klog = -2

	!? SET KANG = 1 FOR DETAILED PRINTOUT OF ANGULAR COEFFICIENTS.
	kang = 0
	sq3 = sqrt(3.q0)
	mla = lrgla
	mlb = lrglb
	lab = max0(lrgla,lrglb)
	md = - mla + mlb

	do i=1,10
		hterm(i) = 0.0
		htscr(i) = 0.0
	end do

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
					!? 'B' OR '3' OR '4' IN VARIOUS PLACES (3 MEANS ELECTRON 1 AND 4 MEANS ELECTRON 2).

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

							lcalc = .false.

							if(nbl1.le.ncut.and.nbl3.le.ncut.and.(idx.eq.1.or.lrgla.eq.0)) then 
								lcalc = .true.
							end if 

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

							sgnx = 1
							
							if(nspna.ne.nspnb) then 
								sgnx = -1
							end if 

							ya1 = ya(ide,nb1)
							ya2 = ya(idep,nb1)
							yb1 = yb(jde,nb3)
							yb2 = yb(jdep,nb3)
							y11 = ya1 + yb1
							y22 = ya2 + yb2
							y11m = ya1 - yb1
							y22m = ya2 - yb2
							y1122 = y11 + y22
							l1 = la(ide,nb1)
							l2 = la(idep,nb1)
							l3 = lb(jde,nb3)
							l4 = lb(jdep,nb3)
							ll1 = l1*(l1+1)
							ll3 = l3*(l3+1)
							ll4 = l4*(l4+1)

							!? CALCULATION OF ANGULAR COEFFICIENTS.  SEE 'CROSS' FOR FURTHER DETAILS.
							call cross(l3,l4,l1,l2,lrglb,lrgla,kt,ide,jde)
							
							write(4,*) 'ok1',cplt(1,1),cplt(2,1)
							
							mlp1 = mlt(1)
							lp1 = lt(1)
							mlp2 = mlt(2)
							lp2 = lt(2)
							
							if(kang.ne.0) then 
								write(6, fmt3) l3,l4,l1,l2,lrglb,lrgla,kt,ide,jde,nb3,nb4,nb1,nb2
								write(6, fmt4)
								
								do ic=1,6
									if(lt(ic).gt.0) then 
										ic2 = lt(ic)
										write(6, fmt4) ic,(cplt(ic3,ic), ic3=  1,ic2)
									end if  
								end do
							end if 

							if(jdx.ne.2) then 

								y11 = ya(ide,nb1) + yb(jde,nb3)
								y22 = ya(idep,nb1) + yb(jdep,nb3)
								kex = 0
								
								if(y22.gt.y11) then
									y22 = ya(ide,nb1) + yb(jde,nb3)
									y11 = ya(idep,nb1) + yb(jdep,nb3)
									kex = 1
								end if

								!?   SET MP1 = 1 - (LOWEST POWER OF R1 OR R2 NEEDED)
								!?   SET MS1 = 1 - (LOWEST POWER OF R12 NEEDED)
								!?   SET IMIN = LOWEST SUM OF POWERS OF R1, R2 AND R12 NEEDED.
								lmax = max0(lp1,lp2)
								mp1 = 6
								ms1 = 1
								imin = -3
								msl = lt(1) + msld
								imin = imin + lrgla + lrglb
								ltota = lrgla + neiga - 1
								ltotb = lrglb + neigb - 1
								lp12 = max0(lp1,lp2)
								
								if(nb1.ne.1.and.nb3.ne.1) then 
									isum = ia + ib + isumd

									!? CALCULATE THE HIGHEST POWERS OF R1, R2 AND R12 NEEDED.
									isum = ijkb(nbl3) + ijka(nbl1) + isumd
									mp2 = iib(nbl3) + (1-id)*iia(nbl1) + id*jja(nbl1) + mp1 + mp2d+ 1 + id
									mq2 = jjb(nbl3) + (1-id)*jja(nbl1) + id*iia(nbl1)+ mp1+mq2d
									ms2 = ms1 + kkb(nbl3) + kka(nbl1) + 2*lp12 + ms2d
									
									escapeCode = .true.
								end if 

								if((nb3.eq.1.and.nb1.eq.1) .and. (.not. escapeCode)) then 
									isum = ltotb + ltota + isumd +2
									mp2 = (1-id)*ltota + id*(ltotb+ltota) + mp1 + mp2d + 1 - id +2
									mq2 = (1-id)*(ltotb+ltota) + id*ltota + mp1 + mq2d + id + 1 +2
									ms2 = ms1 + 2*lp12 + ms2d + 2*id +2
								end if 

								if((nb3.ne.1) .and. (.not. escapeCode)) then 
									isum = ltota + ijkb(nbl3) + isumd +2
									mp2 = iib(nbl3) + id*ltota + mp1 + mp2d + 1 +2
									mq2 = jjb(nbl3) + (1-id)*ltota + mp1 + mq2d +2
									ms2 = ms1 + kkb(nbl3) + 2*lp12 + ms2d +2
									
									escapeCode = .true.
								end if 

								if(.not. escapeCode) then
									isum = ltotb + ijka(nbl1) + isumd +2
									mp2 = (1-id)*iia(nbl1) + id*(jja(nbl1)+ltotb) + mp1 + mp2d + 1 -id +2
									mq2 = (1-id)*(jja(nbl1)+ltotb) + id*iia(nbl1) + mp1 + mq2d + id +1 +2
									ms2 = ms1 + kka(nbl1) + 2*lp12 + ms2d +2
								end if 

								if(id.eq.0.and.(nbx4.gt.nbx3.or.nbx2.gt.nbx1)) then 
									ms2 = ms2 + 2
								end if 

								if(id.eq.0.and.(nbx4.gt.nbx3.and.nbx2.gt.nbx1)) then 
									ms2 = ms2 + 2
								end if 

								if(id.eq.1.and.nbx2.gt.nbx1.and.l3.gt.0) then 
									ms2 = ms2 + 2
								end if 

								lmax = lmax + 2

								if(kex.eq.1) then
									mtemp = mp2
									mp2 = mq2
									mq2 = mtemp
								end if

								if(nbl3.eq.nbl1) then 
									write(*, fmt5) kex,mp1,mp2,mq2,ms1,ms2,isum,y11,y22,l3,l4,l1,l2,nbl1,nbl3
								end if 


								!? GENERATE ARRAYS OF ALL RADIAL INTEGRALS NEEDED.
								call genint(id,nbl3,nbl1,lmax,fac,z)

							end if 

							jdxx = jdx

							if(kex.eq.1) then 
								jdxx = 3 - jdx
							end if 

							!? SUMMATION OVER RIGHT HAND STATE BASIS FUNCTIONS.
							do i=nb1,nb2
								ipa = pa(ide,i)
								iqa = pa(idep,i)
								isa = sa(i)

								!? SUMMATION OVER LEFT HAND STATE BASIS FUNCTIONS.
								do j=nb3,nb4
									
									!? MULTIPLY PSI2 BY 2 IF SUM OVER JDX IS OMITTED FOR SYMMETRIZED OPERATORS.
									psi2 = da(i)*db(j)*sgn
									
									!? TAKE ADVANTAGE OF SYMMETRY FOR DIAGONAL MATRIX ELEMENTS.
									if(izb.ne.0.or.j.ge.i) then 

										if(izb.eq.0.and.j.ne.i) then 
											psi2 = psi2*2.
										end if 

										ipb = pb(jde,j)
										iqb = pb(jdep,j)
										isb = sb(j)
										p = ipa + ipb
										q = iqa + iqb
										s = isa + isb

										!? THE NEXT IS THE CRUCIAL LINE WHERE THE MATRIX ELEMENT IS ASSEMBLED.
										!? 'SPL' COMBINES RADIAL WITH ANGULAR INTEGRALS AND SUMS OVER LAMBDA AS SHOWN
										!? IN EQ. (13) OF DRAKE, PHYS. REV. A 18, 820 (1978).
										!? THE FIRST THREE ARGUMENTS OF 'SPL' ARE THE TOTAL POWERS OF R1, R2 AND R12.
										!? THE FOURTH PICKS WHICH SET OF ANGULAR COEFFICIENTS FROM 'CROSS' TO USE.
										t1 = spl(p+kt-5,q,s,1)
										t2 = spl(p,q+kt-5,s,2)
										sum = spl(p+kt-5,q,s,1) + sgnx*spl(p,q+kt-5,s,2)

										if(abs(t1).gt.1.d70.or.abs(t2).gt.1.d70) then
											write(*, fmt41) i,j,p,q,s,sum,t1,t2
											pause
										endif

										hterm(1) = hterm(1) + sum*psi2

										if(lgo) then 
											htscr(1) = htscr(1) + sum*psi2
										end if 
									end if 
								end do
							end do
						end do
					end do
				end do
			end do
		end do
	end do

	return 
end subroutine matel

