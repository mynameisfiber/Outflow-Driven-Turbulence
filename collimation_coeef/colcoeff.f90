PROGRAM colcoeff
	REAL :: sum, osoft=.1, mu,  PI = 3.1415926535897932384626433832795029, cosphi, collen, r
	INTEGER :: i,j,k, or=6
	REAL, DIMENSION(3) :: col = (/ 1.0, 0.0, 0.0 /)
	REAL :: oi = 10, oj = 10, ok = 10
	
	collen = (col(1)**2 + col(2)**2 + col(3)**2)**.5
	sum = 0
  DO k=ok-or,ok+or
    DO j=oj-or,oj+or
      DO i=oi-or,oi+or
			  if (i .ne. 0 .and. j .ne. 0 .and. k .ne. 0) then
			    r = (i**2 + j**2 + k**2)**.5
			    cosphi = DOT_PRODUCT((/ i-oi,j-oj,k-ok /),col) / (r*collen)
			    !if we are at the asymptote, move over the smallest discrete angle
			    IF (cosphi .ge. 1) THEN
			      print*,"+++",cosphi
			      cosphi = 1-atan(1.0/or)
			    ELSE IF (cosphi .le. -1) THEN 
  			    print*,"---",cosphi
  		      cosphi = atan(1.0/or) - 1
			    END IF
          mu = ACOS(cosphi) * 2.0 / PI - 1
          sum = sum + 1 / (1 + osoft**2 - mu**2)
          print*,1 / (1 + osoft**2 - mu**2), cosphi
			    
        !           dot = DOT_PRODUCT((/ i-or,j-or,k-or /),(/ 1.0,0.0,0.0 /)) / ((i**2 + j**2 + k**2)**.5 * (3*.5**2)**.5 )
        !           mu = 1.0
        !           if (abs(dot) .lt. 1.0) THEN
        !             mu = ACOS(dot) * 2.0 / PI - 1
        !           end if
        !   sum = sum + 1 / (1 + osoft**2 - mu**2)
        ! end if
        ! if (j .eq. or) then
        !   print*,1 / (1 + osoft**2 - mu**2)
  			end if
			end do
		end do
	end do
  print*,1.0/sum
END PROGRAM colcoeff

! import numpy, pylab as py
! a = file("testdata.dat").readlines()
! a = [float(x) for x in a]
! b = numpy.array(a)
! b = b.reshape((12,12))
! py.pcolor(b)
! py.show()
