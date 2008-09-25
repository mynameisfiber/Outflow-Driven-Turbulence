MODULE analysis
  IMPLICIT NONE
  
  INTEGER, PRIVATE :: binsv, numbins, node
  REAL, PRIVATE :: charv
  REAL, DIMENSION(:), ALLOCATABLE, PRIVATE :: sigma, sigmavx
  
  CONTAINS
  
  SUBROUTINE analysis_init(inumbins,ibinsv,icharv,inode)
    INTEGER :: inumbins, ibinsv, inode, pos
    REAL :: icharv
    
    node = inode
    numbins = inumbins
    binsv = ibinsv
    charv = icharv
    ALLOCATE(sigma(numbins*binsv))
    ALLOCATE(sigmavx(2*numbins*binsv))
  
  END SUBROUTINE analysis_init
  
  SUBROUTINE analysis_calc(u,n,ghost,t,nstep)
    INTEGER :: n, ghost, k,j,i, nstep, pos
    REAL :: t, kinetic=0, compressional=0, meanrho=0, meanmom = 0
    REAL, DIMENSION(n,n,n,4) :: u
    
    !PRINT*,"Analyisis: node-",node,"nstep-",nstep,"t-",t
  
    !First we get some values ready for the main calculations
    sigma = 0
    sigmavx = 0
    meanrho = SUM(u(:,:,:,1))/(n*n*n)
    
    !$OMP PARALLEL DO SCHEDULE(STATIC) &
    !$OMP shared(u,n,ghost,kinetic,compressional,sigma,sigmavx,charv,binsv,numbins) &
    !$OMP PRIVATE(i,j,k,pos,meanmom,meanrho) DEFAULT(none)
    do k=ghost+1,n-ghost
      do j=ghost+1,n-ghost
        do i=ghost+1,n-ghost
        
          meanmom = (u(i,j,k,2)**2 + u(i,j,k,3)**2 + u(i,j,k,4)**2)**(.5)
        
          !Calculate kinetic energy
          kinetic = kinetic + .5 * meanmom**2 / u(i,j,k,1)
              
          !Calculate compressional energy
          compressional = compressional + u(i,j,k,1)*log(u(i,j,k,1)/meanrho)
          
          !Update sigma
          pos = INT(meanmom / (charv * u(i,j,k,1)) * binsv)+1
          IF (pos .LT. numbins*binsv .AND. pos .GT. 0) THEN
            sigma(pos) = sigma(pos) + u(i,j,k,1)
          ELSE IF (pos .GT. 0) THEN
            pos = numbins*binsv
            sigma(pos) = sigma(pos) + u(i,j,k,1)
          END IF
          
          !update sigmavx
          pos = INT(u(i,j,k,2) / (charv * u(i,j,k,1)) * binsv + charv*binsv)+1
          IF (pos .LT. 2*numbins*binsv .AND. pos .GT. 0) THEN
            sigmavx(pos) = sigmavx(pos) + u(i,j,k,1)
          ELSE IF (pos .GT. 0) THEN
            pos = 2*numbins*binsv
            sigmavx(pos) = sigmavx(pos) + u(i,j,k,1)
          END IF
          
        end do
      end do
    end do
    
    kinetic = kinetic / (n*n*n)
    compressional = compressional / (n*n*n)
    sigma = sigma
    
    call outputone("kinetic",-1,kinetic,t,nstep)
    call outputone("compress",-1,compressional,t,nstep)
    call outputn(numbins*binsv,"sigma",nstep,sigma,t,nstep)
    call outputn(2*numbins*binsv,"sigmavx",nstep,sigmavx,t,nstep)
  
  END SUBROUTINE analysis_calc
  
  SUBROUTINE outputone(fileprefix, filesuffix, value, t, nstep)
    CHARACTER*54 filename
    CHARACTER*(*) fileprefix
    INTEGER :: filesuffix, nstep
    REAL :: value, t
    
    if (node .eq. 0) then
      WRITE(filename,800) TRIM(fileprefix)
      800 format('output-times-',A)
      OPEN(UNIT=2, FILE=TRIM(filename), ACCESS='APPEND')
      WRITE(2,805) t,nstep
      805 FORMAT(E15.6,' ',I10.10)
      CLOSE(2)
    end if
    
    IF (filesuffix .eq. -1) THEN
      WRITE(filename,810) TRIM(fileprefix), node
      810 format('output-',A,'-',I3.3)
      OPEN(UNIT=2,FILE=TRIM(filename),ACCESS='APPEND')
    ELSE
      WRITE(filename,820) TRIM(fileprefix), filesuffix, node
      820 format('output-',A,'-',I8.8,'-',I3.3)
      OPEN(UNIT=2,FILE=TRIM(filename))
    END IF
    
    WRITE(2,"(E15.6)") value
    CLOSE(2)
  
  END SUBROUTINE outputone
  
  SUBROUTINE outputn(size, fileprefix, filesuffix, value, t, nstep)
    CHARACTER*54 filename
    CHARACTER*(*) fileprefix
    INTEGER :: filesuffix, nstep, size, i
    REAL :: t
    REAL, DIMENSION(size) :: value
    
    if (node .eq. 0) then
      WRITE(filename,900) TRIM(fileprefix)
      900 format('output-times-',A)
      OPEN(UNIT=2, FILE=TRIM(filename), ACCESS='APPEND')
      WRITE(2,905) t,nstep
      905 FORMAT(E15.6,' ',I10.10)
      CLOSE(2)
    end if
    
    IF (filesuffix .eq. -1) THEN
      WRITE(filename,910) TRIM(fileprefix), node
      910 format('output-',A,'-',I3.3)
      OPEN(UNIT=2,FILE=TRIM(filename),ACCESS='APPEND')
    ELSE
      WRITE(filename,920) TRIM(fileprefix), filesuffix, node
      920 format('output-',A,'-',I8.8,'-',I3.3)
      OPEN(UNIT=2,FILE=TRIM(filename))
    END IF
    
    DO i=1,size
      WRITE(2,"(E15.6)") value(i)
    END DO
    
    CLOSE(2)
  
  END SUBROUTINE outputn

END MODULE analysis
