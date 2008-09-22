MODULE analysis
  IMPLICIT NONE
  INTEGER, PARAMETER :: bins=100
  REAL, DIMENSION(bins) :: sigma
  
  CONTAINS
  
  SUBROUTINE analyse(u,n,ghost)
    INTEGER :: n, ghost, k,j,i
    REAL :: kinetic=0, compressional=0, meanrho=0, meanvel = 0
    REAL, DIMENSION(n,n,n,4) :: u
  
    !First we get some values ready for the main calculations
    meanrho = SUM(u(:,:,:,1))/(n*n*n)
    
    do k=ghost+1,n-ghost
      do j=ghost+1,n-ghost
        do i=ghost+1,n-ghost
        
          meanvel = (u(i,j,k,2)**2 + u(i,j,k,3)**2 + u(i,j,k,4)**2)**(.5)
        
          !Calculate kinetic energy
          kinetic = kinetic + .5 * meanvel / u(i,j,k,1)
              
          !Calculate compressional energy
          compressional = compressional + u(i,j,k,1)*log(u(i,j,k,1)/meanrho)
          
          !Calculate sigma
        end do
      end do
    end do
  
  END SUBROUTINE analyse
  
  SUBROUTINE outputone(fileprefix, filesuffix, value)
    CHARACTER*64 fileprefix, filename
    INTEGER :: filesuffix
    REAL :: value
    
    IF (filesuffix .eq. -1) THEN
      WRITE(filename,800) fileprefix
      800 format('output-',S)
      OPEN(UNIT=2,FILE=TRIM(filename),ACCESS='APPEND')
    ELSE
      WRITE(filename,801) fileprefix, filesuffix
      801 format('output-',S,'-',I3.3)
      OPEN(UNIT=2,FILE=TRIM(filename))
    END IF
    WRITE(1,"(E15.6)") value
    CLOSE(2)
  
  END SUBROUTINE outputone

END MODULE analysis
