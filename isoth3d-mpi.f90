program gridcontrol

use physicalparams
use analysis
implicit none
include "mpif.h"
include "omp_lib.h"

type olist
    type(olist), pointer :: next
    real, dimension(3) :: orient
    integer, dimension(3) :: coords
    logical :: filled = .False.
end type olist

! Variable Definitions:
!     procs - number of processes
!     n - length of local node
!     MAXVEL - maximum v_char to inject
!     MAXINJECT - soft maximum outflow restriction... only implemented at the BEGINNING of a timestep.
!     MAXTIME - max number of merging times to go
!     globaln - length of global grid
!     ghost - number of ghost cells
!     CFL - courant number
!     outputfreq - frequency of stdout dumps
!     snapshotfreqt - fractions of a merging time to output (analysis, slice, snapshot)
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
!     colnorm - normalization for outflow
REAL, PARAMETER :: PI = 3.1415926535897932384626433832795029
integer, parameter :: MAXTIME = 20, MAXINJECT=0  !maxtime ~6hr 
real, parameter :: MAXVEL = 75.0
INTEGER, DIMENSION(3) :: snapshotfreqnstep = (/ 0,0,0 /)
REAL, DIMENSION(3) :: lastsavet = 0
LOGICAL :: doanalysis = .false.
real, parameter :: sqrt2=sqrt(2.0)
integer, DIMENSION(3) :: dims
REAL :: dt, t
real, DIMENSION(n,n,n,4) :: u
integer, DIMENSION(3) :: offset, coords, isedge
integer :: globaln, tmp=0
integer :: ierr, COMM_CART
character*64 :: returnmsg = ""
REAL(8) :: cputimeoffset, cputime, tcputime
integer :: nstep=1, numinject = 0, newinject=0
REAL :: colnorm = -1.0

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
call srand(seed+node)

!Initialize the analysis module
call analysis_init(30,20,oImp**(4.0/7.0)*(oSnorm0)**(3.0/7.0),node)

!Initialization
CALL setup(u,n)
CALL timestep(dt, u,n,CFL)
snapshotfreqt = snapshotfreqt / ((oImp**(3.0)) * (oSnorm0**(4.0)))**(1.0/7.0)
t=0.0;

call MPI_BARRIER(MPI_COMM_WORLD,ierr)
do while (TRIM(returnmsg) .eq. "")
  if (node .eq. 0) print*,"Starting nstep=",nstep,"t=",t
  
  !Do the boundry conditions
  call boundary(u,n,ghost)
  
  !Analyse grid and do outflows
  CALL gridanalysis(u,n,dt,t,newinject)
  CALL MPI_Allreduce ( newinject, newinject, 1, &
                     MPI_INTEGER, MPI_SUM, COMM_CART, ierr)
  numinject = numinject + newinject
  newinject = 0

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
      !call analysis_calc(u,n,ghost,t,nstep)
      doanalysis = .true.
      lastsavet(1) = t
  END IF
  
  if (OUTPUTFREQ .NE. 0 .AND. MOD(nstep, OUTPUTFREQ) .EQ. 0 ) &
      CALL output_stdout(nstep,cputime,numinject,t,dt,minval(u(:,:,:,1)))

  !Check if it's time to leave
  IF (MAXWALLTIME .NE. 0 .AND. MAXWALLTIME .LE. cputime) THEN 
    returnmsg = "System Time Limit"
  ELSE IF (MAXTIME .NE. 0 .AND. MAXTIME .LE. t*((oImp**(3.0))*(oSnorm0**(4.0)))**(1.0/7.0)) THEN 
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
call analysis_final()
call MPI_FINALIZE(ierr)
call EXIT(0)


CONTAINS

  SUBROUTINE gridanalysis(u,n,dt,t,newinject)
    INTEGER :: n, oi,oj,ok, newinject, i,j,k !,events
    INTEGER :: sendto, lastsend=0, ti,tj,tk, si,sj,sk
    REAL, DIMENSION(6) :: sendrecv
    LOGICAL :: finished = .false.
    INTEGER stat(MPI_STATUS_SIZE)
    REAL, DIMENSION(n,n,n,4) :: u
    REAL :: dt, ci=0.0,cj=0.0,ck=0.0,t
    type(olist), pointer :: ll, cur,lst,cln
    INTEGER :: numdone = 0
    
    allocate(ll)
    cur => ll
    
    if (doanalysis) call analysis_calc_init(u,n,ghost)

    numdone = 0
    !$OMP PARALLEL DO SCHEDULE(STATIC) &
    !$OMP shared(u,numinject,n,dt,doanalysis,MPI_STATUS_IGNORE,offset,lastsend,comm_cart,cur,numdone,newinject) &
    !$OMP PRIVATE(i,j,k,oi,oj,ok,ci,cj,ck,ierr,sendto,sendrecv,stat,si,sj,sk) DEFAULT(none)
    DO k=ghost+1,n-ghost
      DO j=ghost+1,n-ghost
        DO i=ghost+1,n-ghost
        
          if (doanalysis) call analysis_calc_cell(u(i,j,k,:))
          
          !events = 0
!exp(-1.0*dt/(oSnorm0*p**on))*(dt/(oSnorm0*p**on))**events/factorial(events)
          IF ((MAXINJECT .EQ. 0 .OR. numinject .LT. MAXINJECT) .AND. &
          rand() .LE. dt*oSnorm0*(u(i,j,k,1)**on)) THEN

            !Find coordinates of the outflow
            oi = ANINT(rand()*(n-2*ghost-1)+ghost+1)
            oj = ANINT(rand()*(n-2*ghost-1)+ghost+1)
            ok = ANINT(rand()*(n-2*ghost-1)+ghost+1)
            
            !Create random orientation.  We still poll rand() even if softangle==0
            ! to keep the random number file sync'd for different runs.
            ci = rand()*2-1
            cj = rand()*2-1
            ck = rand()*2-1

            !$OMP CRITICAL
            !create linked list of outflows
            cur%orient = (/ci,cj,ck/)
            cur%coords = (/oi,oj,ok/)
            cur%filled = .True.
            allocate(cur%next)
            cur => cur%next
            cur%next => NULL() 
            !$OMP END CRITICAL
            newinject = newinject + 1
            
            ! CALL generate_outflow(u,n,oi,oj,ok,ci,cj,ck)
            !             !events = events + 1
            !             newinject = newinject + 1
            
            !Now we check if the outflow is near an edge and 
            ! send it to the correct node:
            do tk=-1,1
              do tj=-1,1
                do ti=-1,1
                  if (( ti .eq. 0 .OR. oi+ti*or .GT. n-ghost .OR. oi+ti*or .LE. ghost) .AND. &
                      ( tj .eq. 0 .OR. oj+tj*or .GT. n-ghost .OR. oj+tj*or .LE. ghost) .AND. &
                      ( tk .eq. 0 .OR. ok+tk*or .GT. n-ghost .OR. ok+tk*or .LE. ghost) .AND. &
                      (any((/ti,tj,tk/) .ne. 0 )) ) THEN
                    !$OMP CRITICAL
                    !find out who needs the information
                    call MPI_Cart_rank (COMM_CART, offset+(/ti,tj,tk/), sendto, ierr)
                    !Wait for the send buffer to clear
                    if (lastsend .ne. 0) call MPI_WAIT(lastsend,MPI_STATUS_IGNORE,ierr)
                    !normalize the coordinates for the reciving grid
                    si = abs(n - 2*ghost - ti* oi)
                    sj = abs(n - 2*ghost - tj* oj) 
                    sk = abs(n - 2*ghost - tk* ok) 
                    if (ti .eq. 0) si = oi
                    if (tj .eq. 0) sj = oj
                    if (tk .eq. 0) sk = ok
                    sendrecv = (/ REAL(si),REAL(sj),REAL(sk),ci,cj,ck/)
                    !Send
                    call MPI_ISEND(sendrecv,6,MPI_REAL,sendto,6,COMM_CART,lastsend,ierr)
                    
                    !$OMP END CRITICAL
                  end if
                end do
              end do
            end do
          end if
          
          !check the message buffer and take care of any outstanding outflows
          if (MODULO(k,500) .eq. 0) call blind_outflow_recv(cur,n,numdone)

        END DO
      END DO
    END DO
    
    !Inform neighboors we are done
    do tk=-1,1
      do tj=-1,1
        do ti=-1,1
          if ( any((/ti,tj,tk/) .ne. 0 ) ) THEN
            call blind_outflow_recv(cur,n,numdone)
            if (lastsend .ne. 0) call MPI_WAIT(lastsend,MPI_STATUS_IGNORE,ierr)
            call MPI_Cart_rank (COMM_CART, offset+(/ti,tj,tk/), sendto, ierr)
            call MPI_ISEND((/-1.0,-1.0,-1.0,-1.0,-1.0,-1.0/),6,MPI_REAL,sendto,6,COMM_CART,lastsend,ierr)
          END IF
        END DO
      END DO
    END DO
    
    lst => cur
    ! cur => ll
    ! do while (associated(cur%next) .and. cur%filled .eqv. .True.)
    !   print*,node,"==",cur%coords,"=="
    !   cur=>cur%next
    ! end do
    cur => ll
    do while (numdone .ne. 30) !26 + extra time
      if (associated(cur%next) .and. cur%filled .eqv. .True.) then
        !inject current outflow
        PRINT*,node,"Injecting at ", cur%coords, cur%orient
        CALL generate_outflow(u,n,cur%coords(1),cur%coords(2),cur%coords(3),&
            cur%orient(1),cur%orient(2),cur%orient(3))
        !move on to next outflow
        cur => cur%next
        deallocate(ll)
        ll => cur
      ELSE
        continue
      END IF
      if (numdone .GE. 26) numdone = numdone + 1
      call blind_outflow_recv(lst,n,numdone)
    end do
    
  if (doanalysis) call analysis_calc_end(n,t,nstep)
  doanalysis = .false.
  finished = .false.
  END SUBROUTINE gridanalysis
  
  SUBROUTINE blind_outflow_recv(cur,n, ndone)
    INTEGER :: oi,oj,ok, n, ndone
    REAL :: ci,cj,ck
    LOGICAL :: ismessage = .false.
    REAL, DIMENSION(6) :: recv
    INTEGER stat(MPI_STATUS_SIZE)
    type(olist), pointer :: cur
    
    !$OMP CRITICAL
      1337 continue !hacked up do-while loop
      ismessage = .false.
      call MPI_IPROBE(MPI_ANY_SOURCE,6,COMM_CART,ismessage,MPI_STATUS_IGNORE,ierr)
      IF (ismessage .eqv. .true.) THEN
        call MPI_RECV(recv,6,MPI_REAL,MPI_ANY_SOURCE,6,COMM_CART,stat,ierr)!MPI_STATUS_IGNORE
        
        IF ( ALL(recv .eq. -1.0) ) THEN
          ndone = ndone + 1
          GOTO 1337
        END IF
        
        oi = INT(recv(1))
        oj = INT(recv(2))
        ok = INT(recv(3))
        ci = recv(4)
        cj = recv(5)
        ck = recv(6)
        
        cur%orient = (/ci,cj,ck/)
        cur%coords = (/oi,oj,ok/)
        cur%filled = .True.
        allocate(cur%next)
        cur => cur%next
        cur%next => NULL()
        
        GOTO 1337
      END IF
    !$OMP END CRITICAL
  end subroutine blind_outflow_recv
    

  SUBROUTINE generate_outflow(u,n,oi,oj,ok,ci,cj,ck)
    !The following variables are defined:
    !     oi/oj/ok - global coordinates of injection site's origin
    !     ci/cj/ck - position of the orientation pole
    !     ni/nj/nk - local coordinate of arbitrary point of outflow
    !     x - distance between arbitrary point and center of outflow
    INTEGER :: n,oi,oj,ok, i,j,k, ni,nj,nk
    REAL, DIMENSION(n,n,n,4) :: u
    REAL :: r, ci,cj,ck, V, mu, collen, P, vmax, cosphi
    
    V = (4.0 * PI / 3.0) * or**3
    collen = (ci**2 + cj**2 + ck**2)**.5
    vmax = oImp**(4.0/7.0)*(oSnorm0*u(oi,oj,ok,1)**on)**(3.0/7.0) * MAXVEL
    
    IF (osoft .GT. 0 .AND. colnorm .eq. -1.0) THEN
      colnorm = 0.0
      DO k=ok-or,ok+or
        DO j=oj-or,oj+or
          DO i=oi-or,oi+or
            r = ((i-oi)**2 + (j-oj)**2 + (k-ok)**2)**.5
            IF (r .LE. or .AND. r .NE. 0) THEN
              cosphi = DOT_PRODUCT((/ i-oi,j-oj,k-ok /),(/ ci,cj,ck /)) / (r*collen)
              !if we are at the asymptote, move over the smallest discrete angle
              IF (cosphi .ge. 1) THEN
                cosphi = 1-atan(1.0/or)
              ELSE IF (cosphi .le. -1) THEN 
                cosphi = atan(1.0/or) - 1
              END IF
              mu = ACOS(cosphi) * 2.0 / PI - 1
              colnorm = colnorm + 1/ (1 + osoft**2 - mu**2)
            END IF
          END DO
        END DO
      END DO
      colnorm = V/colnorm
      print*,"colnorm=",colnorm
    END IF

    !$OMP PARALLEL DO SCHEDULE(STATIC) &
    !$OMP shared(globaln,oi,oj,ok,u,ci,cj,ck,V,collen,colnorm,coords,n,offset,vmax,cosphi) &
    !$OMP PRIVATE(i,j,k,r,ni,nj,nk,P,mu) DEFAULT(none)
    DO k=ok-or,ok+or
      DO j=oj-or,oj+or
        DO i=oi-or,oi+or
          r = ((i-oi)**2 + (j-oj)**2 + (k-ok)**2)**.5
          IF (r .LE. or .AND. r .NE. 0 .AND. minval((/i,j,k/)) .GT. 0 &
              .AND. maxval((/i,j,k/)) .LE. n) THEN
            !PRINT*,"DOING INJECTION @ ",i,j,k
            !Take care of collimation
            IF (osoft .GT. 0) THEN
              cosphi = DOT_PRODUCT((/ i-oi,j-oj,k-ok /),(/ ci,cj,ck /)) / (r*collen)
    			    !if we are at the asymptote, move over the smallest discrete angle
    			    IF (cosphi .ge. 1) THEN
    			      cosphi = 1-atan(1.0/r)
    			    ELSE IF (cosphi .le. -1) THEN 
      		      cosphi = atan(1.0/r) - 1
    			    END IF
              mu = ACOS(cosphi) * 2.0 / PI - 1
              P = colnorm / (1 + osoft**2 - mu**2)
            ELSE
              P = 1.0
            END IF
            
            !Mass injection
            u(i,j,k,1) = u(i,j,k,1) + (or-r)*odinj
            
            !Inject correct momentum per unit volume
            ! if (MIN(u(i,j,k,2) + P*oImp/V,vmax*u(i,j,k,1)) .eq. vmax*u(i,j,k,1)) then
            !   print*,"CAPPING INJECT VELOCITY"
            ! END IF
            u(i,j,k,2) = MIN(u(i,j,k,2) + P*oImp/V,vmax*u(i,j,k,1)) * (i-oi)/r
            u(i,j,k,3) = MIN(u(i,j,k,3) + P*oImp/V,vmax*u(i,j,k,1)) * (j-oj)/r
            u(i,j,k,4) = MIN(u(i,j,k,4) + P*oImp/V,vmax*u(i,j,k,1)) * (k-ok)/r
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
    i = INT(n/2)
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

  ! DEPRICATED
  ! REAL FUNCTION rand()
  !   INTEGER :: stat
  !   READ(1,"(F12.10)", IOSTAT=stat) rand
  !   IF (stat .ne. 0) THEN
  !     PRINT*,"An error occured while reading the random file:"
  !     IF (stat .lt. 0) PRINT*,"  End of file"
  !     IF (stat .gt. 0) PRINT*,"  Undetermined"
  !     STOP
  !   END IF
  !   return
  ! END FUNCTION rand
    
end
