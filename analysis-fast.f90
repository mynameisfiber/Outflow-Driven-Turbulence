MODULE analysis
  IMPLICIT NONE
  
  INTEGER :: binsv, numbins, node
  REAL :: charv, meanrho, totalrho, kinetic, compressional
  REAL, DIMENSION(:), ALLOCATABLE :: sigma, sigmavx
  CHARACTER*256 :: aodir
  LOGICAL, DIMENSION(5) :: whichAnalysis

  CONTAINS
  
  
  SUBROUTINE analysis_init(inumbins,ibinsv,icharv,inode,iaodir,iwhichAnalysis)
    INTEGER :: inumbins, ibinsv, inode
    REAL :: icharv
    CHARACTER(LEN=*) :: iaodir
    LOGICAL, DIMENSION(5) :: iwhichAnalysis
    
    aodir = iaodir
    node = inode
    numbins = inumbins
    binsv = ibinsv
    charv = icharv
    whichAnalysis = iwhichAnalysis
    if (whichAnalysis(4)) ALLOCATE(sigma(numbins*binsv))
    if (whichAnalysis(5)) ALLOCATE(sigmavx(2*numbins*binsv))
  
  END SUBROUTINE analysis_init
  
  SUBROUTINE analysis_final()
    if (whichAnalysis(4)) DEALLOCATE(sigma)
    if (whichAnalysis(5)) DEALLOCATE(sigmavx)
  END SUBROUTINE analysis_final
  
  SUBROUTINE analysis_calc_init()
    
    !Get some values ready for the main calculations
    if (whichAnalysis(4)) sigma = 0
    if (whichAnalysis(5)) sigmavx = 0
    compressional = 0
    kinetic = 0
!    if (whichAnalysis(2) .or. whichAnalysis(3) ) THEN
!      !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(u,n,ghost) PRIVATE(k) DEFAULT(none) REDUCTION(+:totalrho)
!      do k=ghost+1,n-ghost
!        totalrho = totalrho + SUM(u(k,ghost+1:n-ghost,ghost+1:n-ghost,1))
!      end do
!      meanrho = totalrho / (n + 2*ghost)**3
!    END IF
  END SUBROUTINE analysis_calc_init
  
  SUBROUTINE analysis_calc_cell(u)
    REAL, DIMENSION(4) :: u
    REAL :: meanmom
    INTEGER :: pos
    
    meanmom = (u(2)**2 + u(3)**2 + u(4)**2)**(.5)
  
    !Calculate kinetic energy
    if (whichAnalysis(1) ) kinetic = kinetic + .5 * meanmom**2 / u(1)
        
    !Calculate compressional energy
    if (whichAnalysis(3) ) compressional = compressional + u(1)*log(u(1)/meanrho)
    
    !Update sigma
    if (whichAnalysis(4) ) THEN
      pos = INT(meanmom / (charv * u(1)) * binsv)+1
      IF (pos .LT. numbins*binsv .AND. pos .GT. 0) THEN
        sigma(pos) = sigma(pos) + u(1)
      ELSE IF (pos .GT. 0) THEN
        pos = numbins*binsv
        sigma(pos) = sigma(pos) + u(1)
      END IF
    END IF
    
    !update sigmavx
    if (whichAnalysis(5) ) THEN
      pos = INT(u(2) / (charv * u(1)) * binsv + charv*binsv)+1
      IF (pos .LT. 2*numbins*binsv .AND. pos .GT. 0) THEN
        sigmavx(pos) = sigmavx(pos) + u(1)
      ELSE IF (pos .GT. 0) THEN
        pos = 2*numbins*binsv
        sigmavx(pos) = sigmavx(pos) + u(1)
      END IF
    END IF
  END SUBROUTINE analysis_calc_cell

  
  SUBROUTINE analysis_calc_end(n,ghost,t,nstep)
    INTEGER :: n, nstep, ghost
    REAL :: t
   
    kinetic = kinetic / (n+2*ghost)**3
    compressional = compressional / (n+2*ghost)**3
    sigma = sigma
  
    if (whichAnalysis(1) ) call outputone("kinetic",-1,kinetic,t,nstep)
    if (whichAnalysis(2) ) call outputone("totalrho",-1,totalrho,t,nstep)
    if (whichAnalysis(3) ) call outputone("compress",-1,compressional,t,nstep)
    if (whichAnalysis(4) ) call outputn(numbins*binsv,"sigma",nstep,sigma,t,nstep)
    if (whichAnalysis(5) ) call outputn(2*numbins*binsv,"sigmavx",nstep,sigmavx,t,nstep)
  END SUBROUTINE analysis_calc_end
  
  SUBROUTINE outputone(fileprefix, filesuffix, value, t, nstep)
    CHARACTER*256 filename
    CHARACTER*(*) fileprefix
    INTEGER :: filesuffix, nstep
    REAL :: value, t
    
    if (node .eq. 0) then
      WRITE(filename,800) TRIM(aodir), TRIM(fileprefix)
      800 format(A,'output-times-',A)
      OPEN(UNIT=2, FILE=TRIM(filename), ACCESS='APPEND')
      WRITE(2,805) t,nstep
      805 FORMAT(E15.6,' ',I10.10)
      CLOSE(2)
    end if
    
    IF (filesuffix .eq. -1) THEN
      WRITE(filename,810) TRIM(aodir), TRIM(fileprefix), node
      810 format(A,'output-',A,'-',I3.3)
      OPEN(UNIT=2,FILE=TRIM(filename),ACCESS='APPEND')
    ELSE
      WRITE(filename,820) TRIM(aodir), TRIM(fileprefix), filesuffix, node
      820 format(A,'output-',A,'-',I8.8,'-',I3.3)
      OPEN(UNIT=2,FILE=TRIM(filename))
    END IF
    
    WRITE(2,"(E15.6)") value
    CLOSE(2)
  
  END SUBROUTINE outputone
  
  SUBROUTINE outputn(size, fileprefix, filesuffix, value, t, nstep)
    CHARACTER*256 filename
    CHARACTER*(*) fileprefix
    INTEGER :: filesuffix, nstep, size, i
    REAL :: t
    REAL, DIMENSION(size) :: value
    
    if (node .eq. 0) then
      WRITE(filename,900) TRIM(aodir), TRIM(fileprefix)
      900 format(A,'output-times-',A)
      OPEN(UNIT=2, FILE=TRIM(filename), ACCESS='APPEND')
      WRITE(2,905) t,nstep
      905 FORMAT(E15.6,' ',I10.10)
      CLOSE(2)
    end if
    
    IF (filesuffix .eq. -1) THEN
      WRITE(filename,910) TRIM(aodir), TRIM(fileprefix), node
      910 format(A,'output-',A,'-',I3.3)
      OPEN(UNIT=2,FILE=TRIM(filename),ACCESS='APPEND')
    ELSE
      WRITE(filename,920) TRIM(aodir), TRIM(fileprefix), filesuffix, node
      920 format(A,'output-',A,'-',I8.8,'-',I3.3)
      OPEN(UNIT=2,FILE=TRIM(filename))
    END IF
    
    DO i=1,size
      WRITE(2,"(E15.6)") value(i)
    END DO
    
    CLOSE(2)
  
  END SUBROUTINE outputn

END MODULE analysis
