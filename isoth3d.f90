!     TODO:
!     
!     BY WEDNESDAY (Aug 27) DO:
!         1) Fix up the normalization parameter in the collimated outflow case
!     
!         2) Make injection time random (use the random distribution chris talks about)
!     
!         3) Have the random numbers come from an external file (See ./MPI for details)
!     
!         4) Find out why outflows aren't as predominant as in the cunningham paper
!
!         5) !!!!!!----MPI----!!!! (See ./MPI for details)
!
!     
!     THEN:
!         +START DOCUMENTING!!!
!     
!         +clean up local/global variables!
!          
!         +optimize OMP
!     


PROGRAM isoth3d
!SIMULATION PARAMS
INCLUDE "omp_lib.h"
INTEGER, PARAMETER :: n = 126, ghost=3
REAL, PARAMETER :: CFL = 0.65
INTEGER :: i0, nsteps=1, outputfreq=200, snapshotfreq=200
INTEGER :: MAXSTEPS=0
REAL :: MAXTIME = 0
REAL*4 rand, second, cputime
INTEGER*4 timeArray(3)
CHARACTER (LEN=56) :: returnmsg = ""
REAL, PARAMETER :: PI = 3.1415926535897932384626433832795029

!PHYSICAL PARAMS
REAL :: dt, t
REAL, DIMENSION( n,n,n, 4) :: u

!OUTFLOW PARAMS
REAL :: Imp = 15012.6, S_norm = 2.91e-4*n**3
REAL :: d_inject = 0.05, softangle = 0.0
REAL :: ci,cj,ck
INTEGER :: r=8, numinject=0, maxinject=0, ii, ij, ik

!------INIT--------
PRINT*,"INITIALIZATION"
ii = rand(Time())
CALL setup(u,n)
CALL timestep(dt, u,n,CFL)
IF (snapshotfreq .NE. 0) CALL printout_file(u,n,0.0,0)
cputime = 0.0

!-----MAINLOOP-----
t = 0; 
PRINT*,"ENTERING LOOP"
PRINT*,""
DO while(returnmsg .EQ. "")
  !!!!!!!!!!!!!!!!!!
  !Outflow Injection
  DO WHILE(INT(t*S_norm) .GT. numinject .AND. (maxinject .EQ. 0 .OR. numinject .LT. maxinject))
    numinject = numinject + 1

    !ii = INT(rand(0)*(n-2*(r+ghost))+r+ghost)
    !ij = INT(rand(0)*(n-2*(r+ghost))+r+ghost)
    !ik = INT(rand(0)*(n-2*(r+ghost))+r+ghost)
    ii = INT(rand(0)*n)
    ij = INT(rand(0)*n)
    ik = INT(rand(0)*n)
    
    ci = rand(0)*2 - 1
    cj = rand(0)*2 - 1
    ck = rand(0)*2 - 1
    
    PRINT 10, numinject,ii, ij, ik, Imp, r, d_inject, t
    10 FORMAT("Outflow #", I6.6," at (i,j,k) = (",I3.3,',',I3.3,',',I3.3, &
        ') with (Imp,r,d_inj,t) = (',F10.2,',',I3.3,',',F10.2,',',F10.5,')')
    PRINT 15, ci,cj,ck,softangle
    15 FORMAT("    (ci,cj,ck, softangle) = (",F10.4,',',F10.4,',',F10.4,',',F10.2,')')
        
    call inject_outflow(n,ghost, ii,ij,ik, ci,cj,ck, Imp, r, d_inject, softangle)
  END DO
  
  !!!!!!!!!!!!!!!!!!
  !Boundary
  call boundary_periodic(u,n,ghost)
  
  !!!!!!!!!!!!!!!!!!
  !Calculate timestep
  CALL timestep(dt, u,n,CFL)
  t = t+2*dt;
  
  !!!!!!!!!!!!!!!!!!
  !Strang splitting of operators
  CALL doX(u,n,dt); 
   CALL doY(u,n,dt); 
    CALL doZ(u,n,dt); 
    CALL doZ(u,n,dt); 
   CALL doY(u,n,dt); 
  CALL doX(u,n,dt);
  
  !!!!!!!!!!!!!!!!!!
  !Timer
  cputime = second()
  
  !!!!!!!!!!!!!!!!!!
  !Output
  IF (MOD(nsteps,outputfreq) .eq. 0) THEN
    PRINT*,"nsteps = ", nsteps
    PRINT*,"cputime = ", cputime
    PRINT*,"numinject = ", numinject
    PRINT*,"t = ", t
    PRINT*,"dt = ", dt
    !PRINT*,"dt_vx = ", (CFL/(maxval(1 + abs(u(:,:,:,2))/u(:,:,:,1))))
    !PRINT*,"dt_vy = ", (CFL/(maxval(1 + abs(u(:,:,:,3))/u(:,:,:,1))))
    !PRINT*,"dt_vz = ", (CFL/(maxval(1 + abs(u(:,:,:,4))/u(:,:,:,1))))
    PRINT*,"min(rho) = ", minval(u(:,:,:,1))
    PRINT*,""
  END IF
  IF (snapshotfreq .NE. 0 .AND. MOD(nsteps,snapshotfreq) .EQ. 0) THEN
    CALL printout_file(u,n,t,nsteps)
  END IF
  
  !!!!!!!!!!!!!!!!!!
  !Check if it's time to leave
  IF (MAXTIME .NE. 0 .AND. MAXTIME .LE. cputime) THEN 
    returnmsg = "System Time Limit"
  ELSE IF (MAXSTEPS .NE. 0 .AND. MAXSTEPS .LE. nsteps) THEN 
    returnmsg = "Timestep Limit"
  ELSE IF (minval(u(:,:,:,1)) .LT. 0) THEN
    returnmsg = "Negative Density"
  END IF
  
  nsteps = nsteps + 1
END DO !i0

!------END-------
IF (snapshotfreq .NE. 0) CALL printout_file(u, n, t, nsteps); 
PRINT*,"-----FINAL REPORT-----"
PRINT*,"nsteps = ", nsteps
PRINT*,"cputime = ", cputime
PRINT*,"numinject = ", numinject
PRINT*,"t = ", t
PRINT*,"dt = ", dt
PRINT*,"dt_vx = ", (CFL/(maxval(1 + abs(u(:,:,:,2))/u(:,:,:,1))))
PRINT*,"dt_vy = ", (CFL/(maxval(1 + abs(u(:,:,:,3))/u(:,:,:,1))))
PRINT*,"dt_vz = ", (CFL/(maxval(1 + abs(u(:,:,:,4))/u(:,:,:,1))))
PRINT*,"min(rho) = ", minval(u(:,:,:,1))
PRINT*,"returnmsg = ", returnmsg
PRINT*,""
stop

contains 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE boundary_periodic(u,n,ghost)
    IMPLICIT NONE
    INTEGER :: n, ghost, a,b,c, inside
    REAL ::  u(n,n,n, 4)
    
    !$omp parallel do schedule(static) shared(u,n,ghost) private(a,b,c) default(none)
    DO c=1,n
      DO b=1,n
        DO a=1,ghost
          !periodic X
          u(a,b,c,:) = u(n-2*ghost+a,b,c,:)
          u(n-a+1,b,c,:) = u(2*ghost-a+1,b,c,:)
        END DO
      END DO
    END DO
    
    !$omp parallel do schedule(static) shared(u,n,ghost) private(a,b,c) default(none)
    DO c=1,n
      DO a=1,ghost
        DO b=1,n
          !periodic Y
          u(b,a,c,:) = u(b,n-2*ghost+a,c,:)
          u(b,n-a+1,c,:) = u(b,2*ghost-a+1,c,:)
        END DO
      END DO
    END DO
    
    !$omp parallel do schedule(static) shared(u,n,ghost) private(a,b,c) default(none)
    DO a=1,ghost 
      DO c=1,n
        DO b=1,n
          !periodic Z
          u(b,c,a,:) = u(b,c,n-2*ghost+a,:)
          u(b,c,n-a+1,:) = u(b,c,2*ghost-a+1,:)
        END DO
      END DO
    END DO
  END SUBROUTINE boundary_periodic

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE timestep(dt, u,n,CFL)! result(dt)
    IMPLICIT NONE
    INTEGER :: n
    REAL, PARAMETER :: cs=1
    REAL ::  CFL, dt, u(n,n,n, 4)
    REAL, DIMENSION(n,n,n) :: rho, rhovx, rhovy, rhovz
    !REAL :: rho(n,n,n), rhovx(n,n,n), rhovy(n,n,n), &
    !   rhovz(n,n,n), e(n,n,n), p(n,n,n), cs(n,n,n)
    rho   = u(:,:,:, 1);
    rhovx = u(:,:,:, 2); 
    rhovy = u(:,:,:, 3);            
    rhovz = u(:,:,:, 4); 
    dt = CFL/(maxval(cs + max(abs(rhovx), abs(rhovy), abs(rhovz))/rho) );
  END SUBROUTINE timestep 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE setup(u,n)
    IMPLICIT NONE
    !REAL, PARAMETER :: Mach = 5; 
    REAL :: u(:,:,:, :)
    INTEGER :: i,j,k, n

    ! constant initial density; 
    u(:,:,:,1) = 1;

    ! Symmetric wall-shock test: mach-5 collision
    !u(1:n/2,:,:, 2) = Mach;  u(n/2+1:n,:,:, 2) = -Mach; 

    ! Zero velocity field
    u(:,:,:, 2) = 0
    u(:,:,:, 3) = 0;   
    u(:,:,:, 4) = 0
  END SUBROUTINE setup

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE inject_outflow(n,ghost, ii,ij,ik, ci,cj,ck, Imp, r, inj, softangle)
    IMPLICIT NONE
    INTEGER :: ii,ij,ik, i,j,k, ni,nj,nk, r, V, n,ghost
    REAL :: x, inj, Imp, softangle, mu, collen, P, colnormalization
    REAL :: ci,cj,ck
    
    V = (4.0 * PI / 3.0) * r**3
    collen = (ci**2 + cj**2 + ck**2)**.5
    colnormalization = 1/log(2/softangle)
    
    !$omp parallel do schedule(static) &
    !$omp shared(u,inj,Imp,V,r,ii,ij,ik,ci,cj,ck,softangle,collen,colnormalization,ghost,n)&
    !$omp private(x,mu,P,ni,nj,nk) default(none)
    DO k=ik-r,ik+r
      DO j=ij-r,ij+r
        DO i=ii-r,ii+r
        x = ((i-ii)**2 + (j-ij)**2 + (k-ik)**2)**.5
        IF (x .LE. r .AND. x .NE. 0) THEN
        
          IF (softangle .GT. 0) THEN
            ! First find the angle separation between the point and the
            !  pole with the angle normalized to -1..1
            mu = ACOS( DOT_PRODUCT((/ i-ii,j-ij,k-ik /),(/ ci,cj,ck /)) &
                 / (x*collen) ) * 2.0 / PI - 1
            P = colnormalization / (1 + softangle**2 - mu**2)
            
            !TFC is an ad-hoc fix for sum(P)!=1
            !P = P/560
          ELSE
            P = 1.0
          END IF

          !Now to normalize i,j,k to the grid
          ni = i; nj = j; nk = k
          IF (ni .LE. ghost) THEN
            ni = n - ghost - abs(ghost - ni)
          ELSEIF (ni .GT. n-ghost) THEN
            ni = ni - n + 2*ghost
          END IF
          IF (nj .LE. ghost) THEN
            nj = n - ghost - abs(ghost - nj)
          ELSEIF (nj .GT. n-ghost) THEN
            nj = nj - n + 2*ghost
          END IF
          IF (nk .LE. ghost) THEN
            nk = n - ghost - abs(ghost - nk)
          ELSEIF (nk .GT. n-ghost) THEN
            nk = nk - n + 2*ghost
          END IF
          
          !PRINT*,"(",i,",",j,",",k,") -> (",ni,",",nj,",",nk,")"

          !density injection
          ! v--------these break momentum conservation
          ! u(i,j,k,2) = u(i,j,k,2) * (1 + (r-x)*inj/u(i,j,k,1))
          ! u(i,j,k,3) = u(i,j,k,3) * (1 + (r-x)*inj/u(i,j,k,1))
          ! u(i,j,k,4) = u(i,j,k,4) * (1 + (r-x)*inj/u(i,j,k,1))
          ! v--------this breaks mass conservation
          u(ni,nj,nk,1) = u(ni,nj,nk,1) + (r-x)*inj
          
          !now actually inject the momentum per unit volume
          u(ni,nj,nk,2) = u(ni,nj,nk,2) + P*Imp/V * (i-ii)/x
          u(ni,nj,nk,3) = u(ni,nj,nk,3) + P*Imp/V * (j-ij)/x
          u(ni,nj,nk,4) = u(ni,nj,nk,4) + P*Imp/V * (k-ik)/x
        END IF
        END DO
      END DO
    END DO
  END SUBROUTINE inject_outflow

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE doX(u,n,dt)
    IMPLICIT NONE
    INTEGER :: n, j, k
    INTEGER, DIMENSION(4) :: reorder = (/1, 2, 3, 4/) 
    REAL u(:,:,:, :); 
    REAL :: dt
    REAL, DIMENSION(n,n, 4) :: u2d
    REAL, DIMENSION(n, 4) :: u1d
    ! X-operation -- i of (i,j,k). 
    !$omp parallel do private(u2d, u1d, k) shared(u,n,dt, reorder) default(none) schedule(dynamic)
    DO j = 1, n
    u2d = u(:,j,:, :); ! pick out xz planes 
      DO k = 1, n
       u1d = u2d(:,k, reorder ) ! pick out x lines from xz planes
       CALL tvdeuler(u1d, n, dt)
       u2d(:,k, reorder ) = u1d;     
      END DO !k
      u(:,j,:, :) = u2d;     
    END DO !j 
  END SUBROUTINE doX
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE doY(u,n,dt)
    IMPLICIT NONE
    INTEGER :: n, j, k
    INTEGER, DIMENSION(4) :: reorder = (/1, 3, 2, 4/) !switch vy&vx
    REAL u(:,:,:, :); 
    REAL :: dt
    REAL, DIMENSION(n,n, 4) :: u2d
    REAL, DIMENSION(n, 4) :: u1d
    ! Y-operation -- j of (i,j,k). 
    !$omp parallel do private(u2d, u1d, k) shared(u,n,dt,reorder) default(none) schedule(dynamic)
    DO j = 1, n
      u2d = u(j,:,:, : ); ! pick out yz planes
      DO k = 1, n
       u1d = u2d(:,k, reorder ) ! pick out y lines from yz plane
       CALL tvdeuler(u1d, n, dt)
       u2d(:,k, reorder ) = u1d;     
      END DO !k
      u(j,:,:, :) = u2d;     
    END DO !j 
  END SUBROUTINE doY

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE doZ(u,n,dt)
    IMPLICIT NONE
    INTEGER :: n, j, k
    INTEGER, DIMENSION(4) :: reorder = (/1, 4, 3, 2 /) !switch vz&vx
    REAL u(:,:,:, :); 
    REAL :: dt
    REAL, DIMENSION(n,n, 4) :: u2d
    REAL, DIMENSION(n, 4) :: u1d
    ! Z-operation -- k of (i,j,k). 
    !$omp parallel do private(u2d, u1d, k) shared(u,n,dt,reorder) default(none) schedule(dynamic)
    DO j = 1, n
      u2d = u(j,:,:, : ); ! pick out yz planes
      DO k = 1, n
       u1d = u2d(k,:, reorder ) ! pick out z lines from yz planes
       CALL tvdeuler(u1d, n, dt)
       u2d(k,:, reorder ) = u1d;     
      END DO !k
      u(j,:,:, :) = u2d;     
    END DO !j 
  END SUBROUTINE doZ

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE tvdeuler(u,n,dt) 
    IMPLICIT NONE
    INTEGER :: n, k, i
    REAL, PARAMETER :: csound = 1
    REAL :: dt
    REAL, DIMENSION(n) :: p,rho,v, rhovy,rhovz, e, c, rhov !,csound
    REAL, DIMENSION(n,4) :: u, uhalf, flux, psi, f0, f
    REAL, DIMENSION(n,4) :: fr, fl, r, cc, fluxR, fluxL
    REAL, DIMENSION(n,4) :: wr, wl, fwr, fwl
    REAL, DIMENSION(n, 3) :: state
    flux = getflux(u); 
    uhalf = u  
    DO k = 2,1,-1 !R-K stepper
      rho=u(:,1);  
      p=rho*csound**2; 
      v=u(:,2)/u(:,1);  !rhovx/rho 
      c =  csound + abs(v); 
      !vvv smooth the freezing speed -- a useful conditioning step. 
      DO i=1,3
         c = max(c, cshift(c,1), cshift(c,-1),cshift(c,2),cshift(c,-2)); 
      ENDDO 
      c = (c+cshift(c,1)+cshift(c,2)+cshift(c,-1)+cshift(c,-2))/5.; 
      !^^^ END of smoothing
      cc = spread(c,2,4); 
      flux = getflux(uhalf); 
      wr = u+flux/cc; wl = u-flux/cc 
      ! TVD is only proven for constant c. One can also try tricks like 
      ! dividing fluxes by some function (e.g., c) to smooth them out, 
      ! THEN DOing tvd on them, THEN multiplying them by c to reconstruct new 
      ! versions. (Not implemented now.)
      fwr= wr*cc;     fwl= -wl*cc; 
      fluxR = tvdflux(wr, fwr); 
      fluxL = cshift(reverse(tvdflux(reverse(wl),reverse(fwl))),1,1);  
      flux = (fluxR + fluxL)/2; 
      uhalf = u-(flux-cshift(flux,-1,1))*dt/k;
    END DO !k
    u = uhalf; 
  END SUBROUTINE tvdeuler

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function reverse(a) result(reversed)
    IMPLICIT NONE
    REAL ::  a(:,:)
    REAL, DIMENSION(size(a,1),size(a,2)) :: reversed 
    INTEGER :: N
    N = size(a,1); 
    reversed(:,:) = a(N:1:-1,:)
  END function reverse

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function tvdflux(u,f) result(flux) 
    IMPLICIT NONE
    !rightward flux -> tvd, 2d-order in x
    REAL :: u(:,:), f(:,:)
    REAL, DIMENSION(size(u,1),size(u,2)) :: flux, fr, fl, r, psi
    fr = (cshift(f,1,1)-f)/2      ! deltaflux-right
    fl = (f - cshift(f,-1,1))/2   ! deltaflux-left
    r=0
    where (fr*fl>0) r = fl/fr; 
    !psi = (r+abs(r))/(1+r)        !van Leer
    psi = max(0.,min(abs(r),1.));   !minmod
    ! psi = max(0., min(2*r,1.),min(r,2.)); !superbee
    flux = f + psi*fr; 
  END function tvdflux

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function getflux(u) result(flux)!derive physical flux from state --EULER
    REAL, PARAMETER :: cssquared = 1 
    REAL :: u(:,:)
    REAL, DIMENSION(size(u,1),size(u,2)) :: flux, cc
    REAL, DIMENSION(size(u,1)) :: rho, rhov, v, p, e
    REAL, DIMENSION(size(u,1), 3) :: state

    rhov = u(:,2); rho = u(:,1); v=rhov/rho;
    flux(:,1) = rhov; 
    flux(:,2) = rhov*v + rho*cssquared; 
    flux(:,3) = u(:,3)*v;  ! passive transport of rho v_y 
    flux(:,4) = u(:,4)*v;  ! passive transport of rho v_z
  END function getflux

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE printout_stdout(u,n) 
    IMPLICIT NONE
    REAL, PARAMETER :: eps = 1.E-10
    INTEGER :: n, i,j,k
    REAL ::  u(:,:,:,:)
    DO i = 1,n
      DO j = 1,n
        DO k = 1, n
          write(*,*) i,j,k, u(i,j,k,1),u(i,j,k,2),u(i,j,k,3),u(i,j,k,4); 
               !i,j,k, rho, rho vx, rho vy, rho vz. 
        END DO !k
      END DO !j
    END DO !i
  END SUBROUTINE printout_stdout

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE printout_file(u,n,t,nsteps)
    REAL, PARAMETER :: eps = 1.E-10
    INTEGER :: n,nsteps, i,j,k
    REAL ::  u(:,:,:,:), t
    CHARACTER*128 filename

    PRINT*,"Outputting to file at nsteps=",nsteps, " & t=",t

    !Create the filename
    WRITE(filename,800) nsteps
    800 format('output-',I8.8)
    OPEN(UNIT=1, FILE=TRIM(filename))
    
    !Output time information
    OPEN(UNIT=2, FILE='outputtimes', ACCESS='APPEND')
    WRITE(2,850) t,nsteps
    850 FORMAT(E15.6,' ',I10.10)
    CLOSE(2)

    !Begin with output
    DO i = 1,n
      DO j = 1,n
        DO k = 1, n
          write(1,900) u(i,j,k,1) , u(i,j,k,2),  u(i,j,k,3), u(i,j,k,4); 
          900 format(E15.6,' ',E15.6,' ',E15.6,' ',E15.6)
               !i,j,k, rho, rho vx, rho vy, rho vz. 
        END DO !k
      END DO !j
    END DO !i
    CLOSE(1)
  END SUBROUTINE printout_file
END 

