program gridcontrol

use analysis
implicit none
include "mpif.h"
include "omp_lib.h"

! Variable Definitions:
!     procs - number of processes
!     n - length of local node
!     MAXVEL - maximum v_char to inject
!     globaln - length of global grid
!     ghost - number of ghost cells
!     CFL - courant number
!     outputfreq - frequency of stdout dumps
!     snapshotfreqt - how many physical seconds to output (analysis, slice, snapshot)
!     snapshotfreqnstep - how many steps to output (analysis, slice, snapshot)
!     dims - # of nodes for each dimension
!     u - test grid
!     ghost - number of ghost zones
!     offset - coordinate in node space
!     coords - origin of local grid in the global grid
!     node - node ID in grid
!     ierr - MPI error code
!     randfile - filename with the random numbers
!     or - radius of outflow
!     islocal - is true if a particular event is local to this node
!     isedge - +1 for top edge, -1 for bottom edge
integer, parameter :: procs = 8, n=30, ghost=3
REAL, PARAMETER :: PI = 3.1415926535897932384626433832795029,  CFL = 0.6
integer, parameter :: MAXTIME = 42840.0, MAXSTEPS = 0, MAXINJECT=0  !maxtime = 5.75 hours
real, parameter :: MAXVEL = 75.0
integer, parameter :: outputfreq=100
INTEGER, DIMENSION(3) :: snapshotfreqnstep = (/ 0,0,0 /)
REAL, DIMENSION(3) :: snapshotfreqt = (/ 1.0/20, 1.0/10, 1.0/2 /) * 1, lastsavet = 0
real, parameter :: sqrt2=sqrt(2.0)
integer, DIMENSION(3) :: dims
REAL :: dt, t
real, DIMENSION(n,n,n,4) :: u
integer, DIMENSION(3) :: offset, coords, isedge
integer :: node, globaln,tmp=0
integer :: ierr, COMM_CART
character*64 :: randfile = "random.txt", returnmsg = ""
REAL(8) :: cputimeoffset, cputime, tcputime
integer :: nstep=1, numinject = 0
REAL, PARAMETER :: oImp=9012.6,oSnorm=2.91e-5*n**3*procs,odinj=0.05,osoft=0.0
INTEGER, PARAMETER :: or = 8

!Check that variables make sense
IF (procs**(1.0/3) .ne. int(procs**(1.0/3))) THEN
  returnmsg = "# Procs must be a perfect cube"
  GOTO 666
END IF
globaln = (procs)**(1.0/3)*n
dims = globaln/n

!Initialize MPI
call MPI_INIT(ierr)
call MPI_Comm_size(MPI_COMM_WORLD, tmp, ierr)
IF (tmp .NE. procs) THEN
  print*,tmp
  returnmsg = "Didn't create enough threads for MPI"
  GOTO 666
END IF
cputimeoffset = MPI_Wtime()
cputime = 0

!Create the cartesian grid and find the local ID in it
call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,(/ .TRUE.,.TRUE.,.TRUE. /), &
     .FALSE.,COMM_CART,ierr)
call MPI_CART_MAP(MPI_COMM_WORLD,3,dims,(/ .TRUE.,.TRUE.,.TRUE. /), &
     node,ierr)

!Collects the local coordinates in node space and calculates the offset
call MPI_CART_COORDS(COMM_CART,node,3,offset,ierr)
coords = offset * n
isedge = INT(2*coords/(globaln-n)-1)
print*,"node=",node,"(x,y,z) = ",coords

!Initialize the analysis module
call analysis_init(30,20,oImp**(4.0/7.0)*(oSnorm/n**3.0)**(3.0/7.0),node)

!Now we open the random number file and set it to IO unit 1
! note: read mode is default
OPEN(UNIT=1, FILE=randfile)


!Initialization
CALL setup(u,n)
CALL timestep(dt, u,n,CFL)
t=0.0;

call MPI_BARRIER(MPI_COMM_WORLD,ierr)
do while (TRIM(returnmsg) .eq. "")
  if (node .eq. 0) print*,"Starting nstep=",nstep,"t=",t
  
  !Manage outflows
  CALL outflow_manager(u,n,dt,numinject)
  
  !Now test out the boundry conditions
  call boundary(u,n,ghost)

  !Calculate timestep
  CALL timestep(dt, u,n,CFL)
  t = t+2*dt;
  if (node .eq. 0) print*,"  dt=",dt
  
  !Strang splitting of operators
  CALL doX(u,n,dt); 
    CALL doY(u,n,dt); 
      CALL doZ(u,n,dt); 
      CALL doZ(u,n,dt); 
    CALL doY(u,n,dt); 
  CALL doX(u,n,dt);

  !update time
  cputime = MPI_Wtime() - cputimeoffset

  !Output
  if ((SNAPSHOTFREQNSTEP(3).NE.0 .AND. MOD(nstep, SNAPSHOTFREQNSTEP(3)).EQ.0) &
  .OR. (SNAPSHOTFREQT(3) .NE. 0 .AND. t .GE. SNAPSHOTFREQT(3) + lastsavet(3))) THEN
      call output_file_cube(u,n,t,nstep)
      lastsavet(3) = t
  END IF
  
  if ((SNAPSHOTFREQNSTEP(2).NE.0 .AND. MOD(nstep, SNAPSHOTFREQNSTEP(2)).EQ.0) &
  .OR. (SNAPSHOTFREQT(2) .NE. 0 .AND. t .GE. SNAPSHOTFREQT(2) + lastsavet(2))) THEN
      call output_file_slice(u,n,t,nstep)
      lastsavet(2) = t
  END IF
  
  if ((SNAPSHOTFREQNSTEP(1).NE.0 .AND. MOD(nstep, SNAPSHOTFREQNSTEP(1)).EQ.0) &
  .OR. (SNAPSHOTFREQT(1) .NE. 0 .AND. t .GE. SNAPSHOTFREQT(1) + lastsavet(1))) THEN
      call analysis_calc(u,n,ghost,t,nstep)
      lastsavet(1) = t
  END IF
  
  if (OUTPUTFREQ .NE. 0 .AND. MOD(nstep, OUTPUTFREQ) .EQ. 0 ) &
      CALL output_stdout(nstep,cputime,numinject,t,dt,minval(u(:,:,:,1)))

  !Check if it's time to leave
  IF (MAXTIME .NE. 0 .AND. MAXTIME .LE. cputime) THEN 
    returnmsg = "System Time Limit"
  ELSE IF (MAXSTEPS .NE. 0 .AND. MAXSTEPS .LE. nstep) THEN 
    returnmsg = "Timestep Limit"
  ELSE IF (minval(u(:,:,:,1)) .LT. 0) THEN
    returnmsg = "Negative Density"
  END IF

  if (node .eq. 0) PRINT*,"DONE"
  nstep = nstep + 1
end do

666 PRINT*,"Node: ", node, "Quiting: ",returnmsg
CALL output_stdout(nstep,cputime,numinject,t,dt,minval(u(:,:,:,1)))

if (node .eq. 0) then
  CALL CPU_TIME(tcputime)
  PRINT*,"Average Walltime per CPU:",cputime/procs
  PRINT*,"Total runtime: ", tcputime
  PRINT*,"Efficency: ", 100*cputime/(procs*tcputime), "%"
end if
call MPI_FINALIZE(ierr)
CLOSE(1)
STOP


CONTAINS

  SUBROUTINE outflow_manager(u,n,dt,numinject)
    INTEGER :: n, oi,oj,ok, numinject, events
    REAL, DIMENSION(n,n,n,4) :: u
    REAL :: dt, ci=0.0,cj=0.0,ck=0.0

    events = 0
    DO WHILE( myrand() .LE. exp(-1.0*dt/oSnorm)*(dt/oSnorm)**events/factorial(events) .AND. &
  	(MAXINJECT .EQ. 0 .OR. numinject .LT. MAXINJECT))
    
      !Find coordinates of the outflow
      oi = ANINT(myrand()*(globaln-2*procs**(1/3.0)*ghost-1)+1)
      oj = ANINT(myrand()*(globaln-2*procs**(1/3.0)*ghost-1)+1)
      ok = ANINT(myrand()*(globaln-2*procs**(1/3.0)*ghost-1)+1)
      
      !Create random orientation.  We still poll myrand() even if softangle==0
      ! to keep the random number file sync'd for different runs.
      ci = myrand()*2-1
      cj = myrand()*2-1
      ck = myrand()*2-1
    
      if (node .eq. 0) print*,"  Creating outflow at: ",oi,oj,ok
      !Now we check if the outflow occures in the local grid.  This is done by
      !   finding the components of the distance at the closest approach and 
      !   seeing if any component is larger than the local grid width
      !NOTE: we could do a more direct checking my finding the distance at 
      !   closest approach to avoid false-positive however it would add 
      !   computational time.
      IF (MAXVAL(MIN(ABS(coords-(/oi,oj,ok/)-2*offset*ghost+n/2.0), &
          ABS(coords-((/oi,oj,ok/)-2*offset*ghost+isedge*globaln)+n/2.0)))-or &
          .LT. n/2) THEN
        CALL generate_outflow(u,n,oi,oj,ok,ci,cj,ck)
      END IF
      
      events = events + 1
      numinject = numinject + 1
    END DO
  END SUBROUTINE outflow_manager

  SUBROUTINE generate_outflow(u,n,oi,oj,ok,ci,cj,ck)
    !The following variables are defined:
    !     oi/oj/ok - global coordinates of injection site's origin
    !     ci/cj/ck - position of the orientation pole
    !     ni/nj/nk - local coordinate of arbitrary point of outflow
    !     x - distance between arbitrary point and center of outflow
    INTEGER :: n,oi,oj,ok, i,j,k, ni,nj,nk
    REAL, DIMENSION(n,n,n,4) :: u
    REAL :: r, ci,cj,ck, V, mu, collen, P, colnorm, vmax
    
    V = (4.0 * PI / 3.0) * or**3
    collen = (ci**2 + cj**2 + ck**2)**.5
    colnorm = .05
    vmax = oImp**(4.0/7.0)*(oSnorm/n**3.0)**(3.0/7.0) * MAXVEL
    
    !First we normalize to coordinates to the grid as to wrap any
    !   outflows around the periodic grid
    oi = oi-coords(1)+2*ghost*offset(1)
    oj = oj-coords(2)+2*ghost*offset(2)
    ok = ok-coords(3)+2*ghost*offset(3)
    
    if (oi .lt. -1*or) then
      oi = oi + (globaln-2*procs**(1/3.0)*ghost)
    else if (oi .gt. n+or) then
      oi = oi - (globaln-2*procs**(1/3.0)*ghost)
    end if
    if (oj .lt. -1*or) then
      oj = oj + (globaln-2*procs**(1/3.0)*ghost)
    else if (oj .gt. n+or) then
      oj = oj - (globaln-2*procs**(1/3.0)*ghost)
    end if
    if (ok .lt. -1*or) then
      ok = ok + (globaln-2*procs**(1/3.0)*ghost)
    else if (ok .gt. n+or) then
      ok = ok - (globaln-2*procs**(1/3.0)*ghost)
    end if
    
    !$OMP PARALLEL DO SCHEDULE(STATIC) &
    !$OMP shared(globaln,oi,oj,ok,u,ci,cj,ck,V,collen,colnorm,coords,n,offset,vmax) &
    !$OMP PRIVATE(i,j,k,r,ni,nj,nk,P,mu) DEFAULT(none)
    DO k=ok-or,ok+or
      DO j=oj-or,oj+or
        DO i=oi-or,oi+or
          r = ((i-oi)**2 + (j-oj)**2 + (k-ok)**2)**.5
          IF (r .LE. or .AND. r .NE. 0) THEN
          
            !Now we check if the given coordinate is inside this grid
            IF (MAXVAL((/i,j,k/)) .LE. n-ghost+1 .AND. &
                MINVAL((/i,j,k/)) .GE. ghost-1) THEN
          
              !Take care of collimation
              IF (osoft .GT. 0) THEN
                mu = ACOS( DOT_PRODUCT((/ i-oi,j-oj,k-ok /),(/ ci,cj,ck /)) &
                     / (r*collen) ) * 2.0 / PI - 1
                P = colnorm / (1 + osoft**2 - mu**2)
                PRINT*,"COLLIMATED!",osoft,mu,(/ci,cj,ck/)
              ELSE
                P = 1.0
              END IF
              
              !Mass injection
              u(i,j,k,1) = u(i,j,k,1) + (or-r)*odinj
              
              !Inject correct momentum per unit volume
              if (MIN(u(i,j,k,2) + P*oImp/V,vmax*u(i,j,k,1)) .eq. vmax*u(i,j,k,1)) then
                print*,"CAPPING INJECT VELOCITY"
              END IF
              u(i,j,k,2) = MIN(u(i,j,k,2) + P*oImp/V,vmax*u(i,j,k,1)) * (i-oi)/r
              u(i,j,k,3) = MIN(u(i,j,k,3) + P*oImp/V,vmax*u(i,j,k,1)) * (j-oj)/r
              u(i,j,k,4) = MIN(u(i,j,k,4) + P*oImp/V,vmax*u(i,j,k,1)) * (k-ok)/r
            END IF
          END IF
        END DO
      END DO
    END DO
    
  
  END SUBROUTINE generate_outflow

  SUBROUTINE boundary(u,n,ghost)
    integer :: n, ghost, dim, ndown, nup, count, ierr, tmp, comm
    integer, dimension(3) :: sds, sde, sus, sue, rd, ru, dims
    real, DIMENSION(n,n,n,4) :: u
    
    !Set up the boundries to be sent.  The variables are defined as:
    !     rd - end of boundry to be recieved from downstream
    !     ru - end of boundry to be recieved from upstream
    !     sds - start of boundry to be sent downstream
    !     sde - end of boundry to be sent downstream
    !     sus - start of boundry to be sent upstream
    !     sue - end of boundry to be sent upstream
    rd = (/ ghost, n, n /)
    ru = n-rd+1
    sds = (/ ghost+1, 1, 1 /)
    sde = (/ 2*ghost, n, n/)
    sus = n - sde + 1
    sue = n - sds + 1
    
    count = n*n*ghost*4
    call MPI_Barrier (COMM_CART,ierr)
    call MPI_Comm_dup(COMM_CART,comm,ierr)
    DO dim=0,2
      !Find neighboors
      call MPI_Cart_shift (comm, dim, -1, &
           ndown, tmp, ierr)
      call MPI_Cart_shift (comm, dim, +1, &
           nup, tmp, ierr)
                          
      !First we send to the down
      !print*,"D",node,"->",ndown
      call MPI_Sendrecv(u(sus(1):sue(1),sus(2):sue(2),sus(3):sue(3),:), count,&
                        MPI_REAL, ndown, 1, u(1:rd(1),1:rd(2),1:rd(3),:), &
                        count, MPI_REAL, ndown, 1, comm, &
                        MPI_STATUS_IGNORE, ierr)
      
      !Now we send to the up
      !print*,"U",node,"->",nup
      call MPI_Sendrecv(u(sds(1):sde(1),sds(2):sde(2),sds(3):sde(3),:), count,&
                        MPI_REAL, nup, 2, u(ru(1):n,ru(2):n,ru(3):n,:), &
                        count, MPI_REAL, nup,2, comm, &
                        MPI_STATUS_IGNORE,ierr)
                       
      !Adjust boundries to be sent for next iteration 
      rd = cshift(rd,2)
      ru = cshift(ru,2)
      sds = cshift(sds,2)
      sde = cshift(sde,2)
      sus = cshift(sus,2)
      sue = cshift(sue,2)
    END DO
    call MPI_Comm_free(comm,ierr)
  END SUBROUTINE boundary
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE timestep(dt, u,n,CFL)! result(dt)
    IMPLICIT NONE
    INTEGER :: n
    REAL, PARAMETER :: cs=1
    REAL ::  CFL, ldt,dt, u(n,n,n, 4)
    REAL, DIMENSION(n,n,n) :: rho, rhovx, rhovy, rhovz
    !REAL :: rho(n,n,n), rhovx(n,n,n), rhovy(n,n,n), &
    !   rhovz(n,n,n), e(n,n,n), p(n,n,n), cs(n,n,n)
    rho   = u(:,:,:, 1);
    rhovx = u(:,:,:, 2); 
    rhovy = u(:,:,:, 3);            
    rhovz = u(:,:,:, 4); 
    ldt = CFL/(maxval(cs + max(abs(rhovx), abs(rhovy), abs(rhovz))/rho) )
    call MPI_Allreduce (ldt, dt, 1, MPI_REAL, MPI_MIN, COMM_CART, ierr )
    
  END SUBROUTINE timestep
  
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
  
  
  SUBROUTINE output_stdout(nstep,cputime,numinject,t,dt,min)
    INTEGER :: nstep,numinject
    REAL :: t,dt,min
    REAL(8) :: cputime
    PRINT*,"Report from node ", node
    PRINT*,node,"   nsteps = ", nstep
    PRINT*,node,"   cputime = ", cputime
    PRINT*,node,"   numinject = ", numinject
    PRINT*,node,"   t = ", t
    PRINT*,node,"   dt = ", dt
    PRINT*,node,"   min(rho) = ", minval(u(:,:,:,1))
  END SUBROUTINE output_stdout
  
  SUBROUTINE output_file_slice(u,n,t,nstep)
    integer :: n,nstep,i,j,k
    CHARACTER*128 filename
    real, dimension(n,n,n,4) :: u
    real :: t
    
    !Write timing info
    if (node .eq. 0) then
      OPEN(UNIT=2, FILE='output-times-slice', ACCESS='APPEND')
      WRITE(2,850) t,nstep
      850 FORMAT(E15.6,' ',I10.10)
      CLOSE(2)
    end if
    
    !Create the filename
    WRITE(filename,800) nstep, node
    800 format('output-slice-',I8.8,'-',I3.3)
    OPEN(UNIT=2, FILE=TRIM(filename))
    print*,"Writing to: ",filename,"@ nstep=",nstep
    k = INT(n/2)
    DO k = ghost+1,n-ghost
      DO j = ghost+1,n-ghost
          write(2,870) u(i,j,k,1) , u(i,j,k,2),  u(i,j,k,3), u(i,j,k,4); 
          870 format(E15.6,' ',E15.6,' ',E15.6,' ',E15.6)
               !i,j,k, rho, rho vx, rho vy, rho vz. 
      END DO !j
    END DO !k
    CLOSE(2)
  END SUBROUTINE output_file_slice
  
  SUBROUTINE output_file_cube(u,n,t,nstep)
    integer :: n,nstep,i,j,k
    CHARACTER*128 filename
    real, dimension(n,n,n,4) :: u
    real :: t
    
    !Write timing info
    if (node .eq. 0) then
      OPEN(UNIT=2, FILE='output-times-cube', ACCESS='APPEND')
      WRITE(2,951) t,nstep
      951 FORMAT(E15.6,' ',I10.10)
      CLOSE(2)
    end if
    
    !Create the filename
    WRITE(filename,901) nstep, node
    901 format('output-cube-',I8.8,'-',I3.3)
    print*,"Writing to: ",filename,"@ nstep=",nstep
    OPEN(UNIT=2, FILE=TRIM(filename))
    DO k = ghost+1,n-ghost
      DO j = ghost+1,n-ghost
        DO i = ghost+1,n-ghost
          write(2,971) u(i,j,k,1) , u(i,j,k,2),  u(i,j,k,3), u(i,j,k,4); 
          971 format(E15.6,' ',E15.6,' ',E15.6,' ',E15.6)
               !i,j,k, rho, rho vx, rho vy, rho vz. 
        END DO !i
      END DO !j
    END DO !k
    CLOSE(2)
  END SUBROUTINE output_file_cube
  
  SUBROUTINE setup(u,n)
    INTEGER :: n
    REAL, DIMENSION(n,n,n,4) :: u
    
    !Constant uniform initial density
    u(:,:,:,1) = 1
    !if (node .eq. 0) u(:,:,:,1) = 1.2
    
    !Zero velocity field
    u(:,:,:,2) = 0
    u(:,:,:,3) = 0
    u(:,:,:,4) = 0
  END SUBROUTINE setup
  
  REAL FUNCTION factorial(n)
    INTEGER :: i,n
    REAL :: result
    
    result = 1
    do i=1,n
      result = result * i
    end do
    factorial = result
    return
  END FUNCTION factorial

  REAL FUNCTION myrand()
    INTEGER :: stat
    READ(1,"(F12.10)", IOSTAT=stat) myrand
    IF (stat .ne. 0) THEN
      PRINT*,"An error occured while reading the random file:"
      IF (stat .lt. 0) PRINT*,"  End of file"
      IF (stat .gt. 0) PRINT*,"  Undetermined"
      STOP
    END IF
    return
  END FUNCTION myrand
    
end
