MODULE analysis
  IMPLICIT NONE
  
  INTEGER :: binsv, numbins, node
  REAL :: charv, meanrho, totalrho, kinetic, compressional
  REAL, DIMENSION(:), ALLOCATABLE :: sigma, sigmavx
  
  CONTAINS
  
  SUBROUTINE analysis_general(u,n,ghost,t,dt,nstep,S0,on,doanalysis)
    INTEGER :: n, ghost, nstep, i,j,k
    REAL, DIMENSION(n,n,n,4) :: u
    REAL :: dt, S0, on, t
    LOGICAL :: doanalysis
    
    if (doanalysis) call analysis_calc_init(u,n,ghost)
    
    do k=ghost+1,n-ghost
      do j=ghost+1,n-ghost
        do i=ghost+1,n-ghost
          
          if (doanalysis) call analysis_calc_cell(u(i,j,k,:))
          
          IF (rand() .LE. dt*S0*u(i,j,k,1)**on) THEN
          
          END IF
          
        end do
      end do
    end do
    
  if (doanalysis) call analysis_calc_end(n,t,nstep)
  END SUBROUTINE analysis_general
  
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
  
  SUBROUTINE analysis_final()
    DEALLOCATE(sigma)
    DEALLOCATE(sigmavx)
  END SUBROUTINE analysis_final
  
  SUBROUTINE analysis_calc_init(u,n,ghost)
    INTEGER :: k, n, ghost
    REAL, DIMENSION(n,n,n,4) :: u
    
    !Get some values ready for the main calculations
    sigma = 0
    sigmavx = 0
    !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(u,n,ghost) PRIVATE(k) DEFAULT(none) REDUCTION(+:totalrho)
    do k=ghost+1,n-ghost
      totalrho = SUM(u(k,ghost+1:n-ghost,ghost+1:n-ghost,1))
    end do
    meanrho = totalrho / (n*n*n)
  END SUBROUTINE analysis_calc_init
  
  SUBROUTINE analysis_calc_cell(u)
    REAL, DIMENSION(4) :: u
    REAL :: meanmom
    INTEGER :: pos
    
    meanmom = (u(2)**2 + u(3)**2 + u(4)**2)**(.5)
  
    !Calculate kinetic energy
    kinetic = kinetic + .5 * meanmom**2 / u(1)
        
    !Calculate compressional energy
    compressional = compressional + u(1)*log(u(1)/meanrho)
    
    !Update sigma
    pos = INT(meanmom / (charv * u(1)) * binsv)+1
    IF (pos .LT. numbins*binsv .AND. pos .GT. 0) THEN
      sigma(pos) = sigma(pos) + u(1)
    ELSE IF (pos .GT. 0) THEN
      pos = numbins*binsv
      sigma(pos) = sigma(pos) + u(1)
    END IF
    
    !update sigmavx
    pos = INT(u(2) / (charv * u(1)) * binsv + charv*binsv)+1
    IF (pos .LT. 2*numbins*binsv .AND. pos .GT. 0) THEN
      sigmavx(pos) = sigmavx(pos) + u(1)
    ELSE IF (pos .GT. 0) THEN
      pos = 2*numbins*binsv
      sigmavx(pos) = sigmavx(pos) + u(1)
    END IF
  END SUBROUTINE analysis_calc_cell

  
  SUBROUTINE analysis_calc_end(n,t,nstep)
    INTEGER :: n, nstep
    REAL :: t
    
    kinetic = kinetic / (n*n*n)
    compressional = compressional / (n*n*n)
    sigma = sigma
  
    call outputone("totalrho",-1,totalrho,t,nstep)
    call outputone("kinetic",-1,kinetic,t,nstep)
    call outputone("compress",-1,compressional,t,nstep)
    call outputn(numbins*binsv,"sigma",nstep,sigma,t,nstep)
    call outputn(2*numbins*binsv,"sigmavx",nstep,sigmavx,t,nstep)
  END SUBROUTINE analysis_calc_end
  
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
