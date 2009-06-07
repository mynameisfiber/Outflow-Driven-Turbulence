program gridcontrol

use params
use analysis
implicit none
include "mpif.h"
include "omp_lib.h"

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
REAL, PARAMETER :: PI = 3.1415926535897932384626433832795029
integer, parameter :: MAXINJECT=0
REAL, DIMENSION(3) :: lastsavet = 0
LOGICAL :: doanalysis = .false.
real, parameter :: sqrt2=sqrt(2.0)
integer, DIMENSION(3) :: dims
REAL :: dt, dteff, dtdone, t, tmerge
real, DIMENSION(n,n,n,4) :: u, ubak
integer, DIMENSION(3) :: offset, coords, isedge
integer :: globaln, tmp=0
integer :: ierr, COMM_CART
character*64 :: returnmsg = ""
CHARACTER*256 :: logfile = ""
REAL(8) :: cputimeoffset, cputime, tcputime
integer :: nstep=1, numinject = 0, isinstable
REAL :: maxv


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
call srand(seed+node)

!Create log file
PRINT*,node,"starting logger @ ",TRIM(odir)
WRITE(logfile,101) odir,node
101 format(A,'nodelog-',I2.2)
OPEN(UNIT=6, FILE=TRIM(logfile), ACCESS='APPEND')
IF (resume) WRITE(6,*),"RESUMING FROM SAVED FILES"
WRITE(6,*),"node=",node,"(x,y,z) = ",coords
FLUSH(6)

!Initialize the analysis module
call analysis_init(30,20,oImp**(4.0/7.0)*(oSnorm0)**(3.0/7.0),node,odir)

!Initialization
tmerge = (op**(3.0/7.0)) / ((oImp**(3.0/7.0)) * (oSnorm0**(4.0/7.0)))
snapshotfreqt = snapshotfreqt * tmerge
IF (resume) THEN
  CALL read_resume_file()
  CALL read_resume_cube(u,n,ghost,nstep)
  call boundary(u,n,ghost)
  lastsavet = t
  WRITE(6,*),"Resumed at t=",t," nstep=",nstep
ELSE
  CALL setup(u,n)
  t=0.0;
END IF

call MPI_BARRIER(MPI_COMM_WORLD,ierr)
do while (TRIM(returnmsg) .eq. "")
  WRITE(6,*),"Starting nstep=",nstep,"t=",t,"(",t/tmerge,"t_merge)"; flush(6)
  
!  !Do the boundry conditions
!  call boundary(u,n,ghost) 
!
  !Analyse grid and do outflows
  CALL gridanalysis(u,n,dt,t,numinject,maxv)

  !Do the boundry conditions
  if (MOD(nstep,bndryfreq) .eq. 0) &
    call boundary(u,n,ghost) 

  !Calculate timestep
  call timestep(dt, maxv, CFL)
  t = t+2*dt;
  !,"REALmaxvcomp",max(maxval(abs(u(:,:,:,2)/u(:,:,:,1))),maxval(abs(u(:,:,:,3)/u(:,:,:,1))),maxval(abs(u(:,:,:,4)/u(:,:,:,1))))

 !Strang splitting of operators
  ubak = u
  dtdone = 0.0
  dteff = dt
  do while (dtdone .lt. dt) 
    
    isinstable = 0
    if (dtdone+dteff .gt. dt) dteff = dt-dtdone

    IF (doX(u,n,dteff,maxv) .lt. 0) THEN
      WRITE(6,*),"Negative at X1"
      isinstable = 1
      ELSE IF (doY(u,n,dteff,maxv) .lt. 0) THEN
        WRITE(6,*),"Negative at Y1"
        isinstable = 1
        ELSE IF (doZ(u,n,dteff,maxv) .lt. 0) THEN
          WRITE(6,*),"Negative at Z1"
          isinstable = 1
        ELSE IF (dteff*(1+maxv) .ge. 1) THEN
          WRITE(6,*),"Courant Condition fail with maxv=",maxv," dt*(1+maxv)=",dteff*(1+maxv)
          isinstable = 2
        ELSE IF (doZ(u,n,dteff,maxv) .lt. 0) THEN
          WRITE(6,*),"Negative at Z2"
          isinstable = 1
        ELSE IF (doY(u,n,dteff,maxv) .lt. 0) THEN
          WRITE(6,*),"Negative at Y2"
          isinstable = 1
      ELSE IF (doX(u,n,dteff,maxv) .lt. 0) THEN
        WRITE(6,*),"Negative at X2"
        isinstable = 1
    END IF 

    IF (isinstable .ne. 0) THEN
      FLUSH(6)
      u = ubak
      dtdone = 0
      !IF (isinstable .eq. 1) THEN
        dteff = dteff / 2
      !ELSE
      !  dteff = CFL/(1+maxv)
      !END IF
      isinstable = 0
      WRITE(6,*),"Instable condition reached. maxv=",maxv,"dt'=",dteff,"numiter=",dt/dteff
      FLUSH(6)
    ELSE
      dtdone = dtdone + dteff
    END IF
  end do
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
      CALL output_stdout(nstep,cputime,numinject,t,dt)

  !Check if it's time to leave
  IF (MAXWALLTIME .NE. 0 .AND. MAXWALLTIME .LE. cputime) THEN 
    returnmsg = "System Time Limit"
  ELSE IF (MAXTIME .NE. 0 .AND. MAXTIME .LE. t*((oImp**(3.0))*(oSnorm0**(4.0)))**(1.0/7.0)) THEN 
    returnmsg = "Timestep Limit"
  ELSE IF (MAXSTEP .NE. 0 .AND. MAXSTEP .LE. nstep) THEN
    returnmsg = "Integration Step Limit"
  END IF

  nstep = nstep + 1
  FLUSH(6)
end do

666 WRITE(6,*),"Node: ", node, "Quiting: ",returnmsg

WRITE(6,*),"Writing resume data"
IF (node .eq. 0) call create_resume_file()
CALL output_file_cube(u,n,t,nstep)

CALL output_stdout(nstep,cputime,numinject,t,dt)

if (node .eq. 0) then
  CALL CPU_TIME(tcputime)
  WRITE(6,*),"Average Walltime per CPU:",cputime/procs
  !WRITE(6,*),"Total runtime: ", tcputime
  WRITE(6,*),"Efficency: ", 100*cputime/(procs*tcputime), "%"
end if
call analysis_final()
call MPI_FINALIZE(ierr)
CLOSE(6)
call EXIT(0)


CONTAINS

  SUBROUTINE gridanalysis(u,n,dt,t,numinject,maxv)
    INTEGER, PARAMETER :: bufsize = 40
    INTEGER :: n, newinject, numinject, i,j,k, idx, d, chk
    LOGICAL :: isinject 
    REAL, DIMENSION(n,n,n,4) :: u
    REAL :: dt, ci=0.0,cj=0.0,ck=0.0,t, oi, oj, ok, rhoaddtot, maxv
    REAL, DIMENSION(7,bufsize) :: lflows
    REAL, DIMENSION(7,bufsize*procs) :: gflows
    INTEGER, DIMENSION(3) :: inj, inloc

    if (doanalysis) call analysis_calc_init(u,n,ghost)

    gflows(:,:) = -1
    lflows(:,:) = -1
    idx = 1 
    maxv = 0
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) &
    !$OMP shared(u,numinject,n,dt,doanalysis,lflows,idx) &
    !$OMP PRIVATE(i,j,k,oi,oj,ok,ci,cj,ck) DEFAULT(SHARED) &
    !$OMP REDUCTION(MAX:maxv)
    DO k=1+ghost,n-ghost
      DO j=1+ghost,n-ghost
        DO i=1+ghost,n-ghost

      !    !VVVVVVVVV------  THE NEXT SECTION FOR ACTIVE GRID ONLY  -----VVVVVVVVVVVV
      !    IF (min(i,j,k) .gt. ghost .AND. max(i,j,k) .LE. n-ghost) THEN
          if (u(i,j,k,1) .lt. 0) THEN
            WRITE(6,*),"Negative: ",i,j,k,"|",u(i,j,k,:)
          END IF
          maxv = MAX(maxv,MAXVAL(abs(u(i,j,k,2:4))/u(i,j,k,1)))

!          if (maxval(u(i,j,k,2:4))/u(i,j,k,1) .gt. MAXVEL) THEN
!            rhotemp = u(i,j,k,1)
!            u(i,j,k,1) = maxval(u(i,j,k,2:4))/MAXVEL 
!            rhoadd = rhoadd + u(i,j,k,1) - rhotemp
!          END IF
          
          if (doanalysis) call analysis_calc_cell(u(i,j,k,:))
          IF ((MAXINJECT .EQ. 0 .OR. numinject .LT. MAXINJECT) .AND. &
          rand() .LE. 2*dt*oSnorm0*(u(i,j,k,1)**on) ) THEN
            !         ^-- multiply by two because one step represents a change
            !             of 2*dt

            oi = REAL(i)
            oj = REAL(j)
            ok = REAL(k)
            
            !Create random orientation.  We still poll rand() even if softangle==0
            ! to keep the random number file sync'd for different runs.
            ci = rand()*2-1
            cj = rand()*2-1
            ck = rand()*2-1
            
            !$OMP CRITICAL
            lflows(:,idx) = (/ oi, oj, ok, ci, cj, ck, REAL(node) /)
            idx = idx + 1
            !$OMP END CRITICAL
           
            !WRITE(6,*),idx,"Queuing ", lflows(:,idx-1); flush(6)

          end if
        END DO
      END DO
    END DO

    CALL MPI_ALLGATHER(lflows,bufsize*7,MPI_REAL,gflows,bufsize*7,MPI_REAL, &
                                COMM_CART,ierr)
   
    ! SHARED(u,node,comm_cart,offset,bufsize,procs,gflows,n,ghost,or) DEFAULT(NONE) &
    newinject = 0
    ! PARALLEL DO DEFAULT(SHARED) &
    ! PRIVATE(idx,d,isinject,inj,ierr,chk,inloc) &
    ! REDUCTION(+:newinject) REDUCTION(MAX:maxv) SCHEDULE(DYNAMIC)
    do idx=1,bufsize*(procs)
      !WRITE(6,*),idx,"|",gflows(:,idx)
      IF (any(gflows(1:3,idx) .ne. -1)) THEN
        newinject = newinject + 1
        d = INT(gflows(7,idx))
        if (d .eq. node) THEN
          inj = gflows(1:3,idx)
          isinject = .true.
        else
          CALL MPI_CART_COORDS(COMM_CART,d,3,inloc,ierr)
          isinject = .true.
          do chk=1,3
            !First check standard grid overlap
            inj(chk) = gflows(chk,idx) - (offset(chk)-inloc(chk))*(n-2*ghost)
            if (inj(chk) .le. n-ghost+or .and. inj(chk) .gt. ghost-or) THEN
              isinject = isinject .and. .True.
            else
              !Now check periodicity
              inj(chk) = gflows(chk,idx) + (offset(chk)-inloc(chk))*(n-2*ghost)
              if (inj(chk) .le. n-ghost+or .and. inj(chk) .gt. ghost-or) THEN
                isinject = isinject .and. .true.
              else
                !Both checks failed... this outflow is not in the grid
                isinject = isinject .and. .false.
              end if
            end if
          end do
        end if
        IF (isinject) THEN
          WRITE(6,*),"Injecting at ",inj,gflows(4:6,idx)
          CALL generate_outflow(u,n,inj(1),inj(2),inj(3),gflows(4,idx),gflows(5,idx),gflows(6,idx),maxv)
        end if
      END IF
    end do
    numinject = numinject + newinject

    rhoaddtot = oImp / MAXVEL * newinject
    u(:,:,:,:) = u(:,:,:,:) * ((n-2*ghost)**3*procs*op / ((n-2*ghost)**3*procs*op + rhoaddtot))
    WRITE(6,*),"injfactor=",((n-2*ghost)**3*procs*op / ((n-2*ghost)**3*procs*op+rhoaddtot))


    if (doanalysis) call analysis_calc_end(n,ghost,t,nstep)
    doanalysis = .false.
  END SUBROUTINE gridanalysis
  
  SUBROUTINE generate_outflow(u,n,oi,oj,ok,ci,cj,ck,maxv)
    !The following variables are defined:
    !     oi/oj/ok - global coordinates of injection site's origin
    !     ci/cj/ck - position of the orientation pole
    !     ni/nj/nk - local coordinate of arbitrary point of outflow
    !     x - distance between arbitrary point and center of outflow
    INTEGER :: n,oi,oj,ok, i,j,k, ni,nj,nk, avg=2, tmp
    REAL, DIMENSION(n,n,n,4) :: u
    REAL :: r, ci,cj,ck, V=0, collen, P, cosphi
    REAL :: colnorm,maxv

    collen = (ci**2 + cj**2 + ck**2)**.5
    IF (osoft .GT. 0)  colnorm = 0.0
    IF (osoft .GT. 0 .OR. V .eq. 0) THEN
      !$OMP PARALLEL DO SCHEDULE(DYNAMIC) &
      !$OMP PRIVATE(i,j,k,cosphi,r) &
      !$OMP SHARED(oi,oj,ok,collen) &
      !$OMP REDUCTION(+:colnorm) REDUCTION(+:V)
      DO k=ok-or,ok+or
        DO j=oj-or,oj+or
          DO i=oi-or,oi+or
            r = ((i-oi+.5)**2 + (j-oj+.5)**2 + (k-ok+.5)**2)**.5
            IF (r .LE. or .AND. r .GT. 0) THEN
              IF (osoft .GT. 0) THEN
                cosphi = DOT_PRODUCT((/ i-oi+.5,j-oj+.5,k-ok+.5 /),(/ ci,cj,ck /)) / (r*collen)
                colnorm = colnorm + 1/ (1 + osoft**2 - cosphi**2)
             END IF
             V = V + 1
            END IF
          END DO
        END DO
      END DO
      if (osoft .gt. 0) colnorm = colnorm/V
    END IF

    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) &
    !$OMP shared(oi,oj,ok,u,ci,cj,ck,V,collen,colnorm,n) &
    !$OMP PRIVATE(i,j,k,r,P,cosphi) DEFAULT(SHARED) &
    !$OMP REDUCTION(MAX:maxv)
    DO k=ok-or-1,ok+or+1
      DO j=oj-or-1,oj+or+1
        DO i=oi-or-1,oi+or+1
          r = ((i-oi+.5)**2 + (j-oj+.5)**2 + (k-ok+.5)**2)**.5
          IF (r .LE. or .AND. r .GT. 0 .AND. min(i,j,k) .GE. 1 &
              .AND. max(i,j,k) .LE. n) THEN
            !WRITE(6,*),"DOING INJECTION @ ",i,j,k
            !Take care of collimation
            IF (osoft .GT. 0) THEN
              cosphi = DOT_PRODUCT((/ i-oi+.5,j-oj+.5,k-ok+.5 /),(/ ci,cj,ck /)) / (r*collen)
              P = 1/(colnorm * (1 + osoft**2 - cosphi**2))
            ELSE
              P = 1
            END IF
            
            !u(i,j,k,1) = MAX(u(i,j,k,1), (P*oImp/V)/MAXVEL)
            u(i,j,k,1) = u(i,j,k,1) + (P*oImp/V)/MAXVEL
            u(i,j,k,2) = u(i,j,k,2) + P * oImp/V * (i-oi+.5)/r
            u(i,j,k,3) = u(i,j,k,3) + P * oImp/V * (j-oj+.5)/r
            u(i,j,k,4) = u(i,j,k,4) + P * oImp/V * (k-ok+.5)/r
          
            maxv = MAX(maxv,MAXVAL(abs(u(i,j,k,2:4))/u(i,j,k,1)))

          END IF
        END DO
      END DO
    END DO

    DO tmp=1,1
      !randomize the direction we average in to reduce any directional bias
      ni = INT(SIGN(1.0, rand() - .5))
      nj = INT(SIGN(1.0, rand() - .5))
      nk = INT(SIGN(1.0, rand() - .5))
      !$OMP PARALLEL DO PRIVATE(r,i,j,k) DEFAULT(SHARED) SCHEDULE(DYNAMIC)
      DO k=ok-nk*(or+1),ok+nk*(or+1),nk
        DO j=oj-nj*(or+1),oj+nj*(or+1),nj
          DO i=oi-ni*(or+1),oi+ni*(or+1),ni
            r = ((i-oi+.5)**2 + (j-oj+.5)**2 + (k-ok+.5)**2)
            IF (r .LT. or+1 .AND. r .GT. 0 .AND. min(i,j,k) .GE. 1 &
                .AND. max(i,j,k) .LE. n) THEN
              !write(6,*),"Before Ave: ",i,j,k,"|",u(i,j,k,:); flush(6);
              !u(i,j,k,2:4) = MAX(MIN(u(i,j,k,2:4),(/1,1,1/)*MAXVEL*velparam(u(i,j,k,1))*u(i,j,k,1)),&
              !                                 -1*(/1,1,1/)*MAXVEL*velparam(u(i,j,k,1))*u(i,j,k,1))
              u(i,j,k,2:4) = SUM(SUM(SUM( &
                                u( maxmin(i-avg,1,n):maxmin(i+avg,1,n), &
                                   maxmin(j-avg,1,n):maxmin(j+avg,1,n), &
                                   maxmin(k-avg,1,n):maxmin(k+avg,1,n), &
                                   2:4) &
                              ,3),2),1)&
                               / ((maxmin(i+avg,1,n)-maxmin(i-avg,1,n)+1) * &
                                  (maxmin(j+avg,1,n)-maxmin(j-avg,1,n)+1) * &
                                  (maxmin(k+avg,1,n)-maxmin(k-avg,1,n)+1) )
 !             IF (tmp .eq. 1 .and. u(i,j,k,1) .lt. 5e-2) WRITE(6,*),"injinfo",i,j,k,"|",u(i,j,k,:)
              !write(6,*),"After Ave: ",i,j,k,"|",u(i,j,k,:); flush(6);
            END IF
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE generate_outflow

  INTEGER FUNCTION maxmin(num,low,high) result(i)
    INTEGER :: high, low, num
    i = min(max(num,low),high)
  END FUNCTION maxmin

  SUBROUTINE boundary(u,n,ghost)
    integer :: n, ghost, dim, ndown, nup, count, ierr, tmp
    integer, dimension(3) :: sds, sde, sus, sue, rd, ru
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
!    IF (ierr .ne. MPI_SUCCESS) THEN
!      WRITE(6,*),"Error at boundary barrier"; flush(6)
!    end if
!    call MPI_Comm_dup(COMM_CART,comm,ierr)
!    IF (ierr .ne. MPI_SUCCESS) THEN
!      WRITE(6,*),"Error at boundary dup"; flush(6)
!    end if
    DO dim=0,2
      !Find neighboors
      call MPI_Cart_shift (COMM_CART, dim, -1, &
           ndown, tmp, ierr)
      call MPI_Cart_shift (COMM_CART, dim, +1, &
           nup, tmp, ierr)
                          
      !First we send to the down
      !WRITE(6,*),"D",node,"->",ndown; FLUSH(6)
      !WRITE(6,*),sus(1),sue(1),",",sus(2),sue(2),",",sus(3),sue(3)
      !WRITE(6,*),1,rd(1),",",1,rd(2),",",1,rd(3)
      call MPI_Sendrecv(u(sus(1):sue(1),sus(2):sue(2),sus(3):sue(3),:), count,&
                        MPI_REAL, ndown, dim*node*ndown, u(1:rd(1),1:rd(2),1:rd(3),:), &
                        count, MPI_REAL, ndown, dim*node*ndown, COMM_CART, &
                        MPI_STATUS_IGNORE, ierr)

      !Now we send to the up
      !WRITE(6,*),"U",node,"->",nup; FLUSH(6)
      !WRITE(6,*),sds(1),sde(1),",",sds(2),sde(2),",",sds(3),sde(3)
      !WRITE(6,*),ru(1),n,",",ru(2),n,",",ru(3),n
      call MPI_Sendrecv(u(sds(1):sde(1),sds(2):sde(2),sds(3):sde(3),:), count,&
                        MPI_REAL, nup, dim*node*nup*2, u(ru(1):n,ru(2):n,ru(3):n,:), &
                        count, MPI_REAL, nup, dim*node*nup*2, COMM_CART, &
                        MPI_STATUS_IGNORE,ierr)


      !Adjust boundries to be sent for next iteration 
      rd = cshift(rd,-1)
      ru = cshift(ru,-1)
      sds = cshift(sds,-1)
      sde = cshift(sde,-1)
      sus = cshift(sus,-1)
      sue = cshift(sue,-1)
    END DO
!    call MPI_Comm_free(comm,ierr)

  END SUBROUTINE boundary
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE timestep(dt, maxv, CFL)!u,n,CFL)! result(dt)
    IMPLICIT NONE
    REAL ::  CFL, ldt,dt, maxv!, u(n,n,n, 4), maxv=0, maxp=0
    !maxvcomp = max(maxval(abs(u(:,:,:,2))/u(:,:,:,1)),maxval(abs(u(:,:,:,3))/u(:,:,:,1)),maxval(abs(u(:,:,:,4))/u(:,:,:,1)))
    
    if (maxv .eq. 0) maxv=MAXVEL !we set maxvcomp to 10 so that the first
                                     !timesteps are not overly large
    ldt = CFL/(1+maxv)
    call MPI_Allreduce (ldt, dt, 1, MPI_REAL, MPI_MIN, COMM_CART, ierr )
    WRITE(6,*),"dt=",dt,"maxv=",maxv
  END SUBROUTINE timestep

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL FUNCTION doX(u,n,dt,maxv) result(minrho)
    IMPLICIT NONE
    INTEGER :: n, j, k
    INTEGER, DIMENSION(4) :: reorder = (/1, 2, 3, 4/) 
    REAL u(:,:,:, :); 
    REAL :: dt, minrho, maxv
    REAL, DIMENSION(n,n, 4) :: u2d
    REAL, DIMENSION(n, 4) :: u1d
    minrho = u(1,1,1,1)
    maxv = 0
    ! X-operation -- i of (i,j,k). 
    !$omp parallel do private(u2d, u1d, k) shared(u,n,dt, reorder) DEFAULT(SHARED) schedule(dynamic) &
    !$omp reduction(MIN:minrho) reduction(MAX:maxv)
    DO j = 1, n
    u2d = u(:,j,:, :); ! pick out xz planes 
      DO k = 1, n
       u1d = u2d(:,k, reorder ) ! pick out x lines from xz planes
       CALL tvdeuler(u1d, n, dt)
       minrho = MIN(minrho, MINVAL(u1d(:,1)))
       maxv = MAX(maxv, MAXVAL(u1d(:,2:4)/SPREAD(u1d(:,1),2,3)))
       u2d(:,k, reorder ) = u1d;     
      END DO !k
      u(:,j,:, :) = u2d;     
    END DO !j 
 END FUNCTION doX

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL FUNCTION doY(u,n,dt,maxv) result(minrho)
    IMPLICIT NONE
    INTEGER :: n, j, k
    INTEGER, DIMENSION(4) :: reorder = (/1, 3, 2, 4/) !switch vy&vx
    REAL u(:,:,:, :); 
    REAL :: dt, minrho, maxv
    REAL, DIMENSION(n,n, 4) :: u2d
    REAL, DIMENSION(n, 4) :: u1d
    minrho = u(1,1,1,1)
    maxv = 0
    ! Y-operation -- j of (i,j,k). 
    !$omp parallel do private(u2d, u1d, k) shared(u,n,dt,reorder) DEFAULT(SHARED) schedule(dynamic) &
    !$omp reduction(MIN:minrho) reduction(MAX:maxv)
    DO j = 1, n
      u2d = u(j,:,:, : ); ! pick out yz planes
      DO k = 1, n
       u1d = u2d(:,k, reorder ) ! pick out y lines from yz plane
       CALL tvdeuler(u1d, n, dt)
       minrho = MIN(minrho, MINVAL(u1d(:,1)))
       maxv = MAX(maxv, MAXVAL(u1d(:,2:4)/SPREAD(u1d(:,1),2,3)))
       u2d(:,k, reorder ) = u1d;     
      END DO !k
      u(j,:,:, :) = u2d;     
    END DO !j 
  END FUNCTION doY

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL FUNCTION doZ(u,n,dt,maxv) result(minrho)
    IMPLICIT NONE
    INTEGER :: n, j, k
    INTEGER, DIMENSION(4) :: reorder = (/1, 4, 3, 2 /) !switch vz&vx
    REAL u(:,:,:, :); 
    REAL :: dt, minrho, maxv
    REAL, DIMENSION(n,n, 4) :: u2d
    REAL, DIMENSION(n, 4) :: u1d
    minrho = u(1,1,1,1)
    maxv = 0
    ! Z-operation -- k of (i,j,k). 
    !$omp parallel do private(u2d, u1d, k) shared(u,n,dt,reorder) DEFAULT(SHARED) schedule(dynamic) &
    !$omp reduction(MIN:minrho) reduction(MAX:maxv)
    DO j = 1, n
      u2d = u(j,:,:, : ); ! pick out yz planes
      DO k = 1, n
       u1d = u2d(k,:, reorder ) ! pick out z lines from yz planes
       CALL tvdeuler(u1d, n, dt)
       minrho = MIN(minrho, MINVAL(u1d(:,1)))
       maxv = MAX(maxv, MAXVAL(u1d(:,2:4)/SPREAD(u1d(:,1),2,3)))
       u2d(k,:, reorder ) = u1d;     
      END DO !k
      u(j,:,:, :) = u2d;     
    END DO !j 
  END FUNCTION doZ

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE tvdeuler(u,n,dt) 
    IMPLICIT NONE
    INTEGER :: n, k
    REAL, PARAMETER :: csound = 1
    REAL :: dt
    REAL, DIMENSION(n) :: v, c
    REAL, DIMENSION(n,4) :: u, uhalf, flux 
    REAL, DIMENSION(n,4) :: cc, fluxR, fluxL
    REAL, DIMENSION(n,4) :: wr, wl, fwr, fwl
    uhalf = u  
    DO k = 2,1,-1 !R-K stepper
      v=u(:,2)/u(:,1);  !rhovx/rho 
      c =  MAXVAL(abs(v) + csound); 
      !vvv smooth the freezing speed -- a useful conditioning step. 
      !DO i=1,3
      !   c = min(c, cshift(c,1), cshift(c,-1),cshift(c,2),cshift(c,-2)); 
      !ENDDO 
      !c = (c+cshift(c,1)+cshift(c,2)+cshift(c,-1)+cshift(c,-2))/5.; 
      !^^^ END of smoothing
      cc = spread(c,2,4); 
      flux = getflux(uhalf); 
      wr = u+flux/cc; wl = u-flux/cc 
      ! TVD is only proven for constant c. One can also try tricks like 
      ! dividing fluxes by some function (e.g., c) to smooth them out, 
      ! THEN DOing tvd on them, THEN multiplying them by c to reconstruct new 
      ! versions. (Not implemented now.)
      !fwr= wr*cc;     fwl= -wl*cc; 
      fwr= u*cc+flux;     fwl= flux-u*cc; 
      fluxR = tvdflux(wr, fwr); 
      fluxL = cshift(reverse(tvdflux(reverse(wl),reverse(fwl))),1,1);  
      flux = (fluxR + fluxL)/2; !flux(n,:) = 0
      !uhalf = u-(flux-cshift(flux,-1,1))*dt/k;
      uhalf(2:n-1,:) = u(2:n-1,:)-(flux(2:n-1,:)-flux(1:n-2,:))*dt/k;
    END DO !k
!      if (minval(uhalf(ghost+1:n-ghost,1)) .lt. 0) then
!        WRITE(6,*),"k=",k
!        call debugoutput(uhalf,u,flux,fluxR,fluxL,n)
!      end if
    u = REAL(uhalf); 
  END SUBROUTINE tvdeuler

  subroutine debugoutput(u,du,v,fr,fl,n)
    integer :: n, k
    REAL, dimension(n,4) :: du,u, v,fl,fr

    write(6,*),"        n       u_old        u        flux          |v|         u+flux          fluxR             fluxL"
    do k=2,n-1
      !if (MOD(k,25) .eq. 0) write(6,*),"           n       u_old        u        flux          |v|         fluxR             fluxL"
      if (PRODUCT(u(k-1:k+1,1)) .lt. 0) THEN
        write(6,*),k,du(k,1),u(k,1),v(k,1),sqrt(du(k,2)**2+du(k,3)**2+du(k,4)**2)/du(k,1),fr(k,1),fl(k,1)
      end if
    end do
  end subroutine debugoutput

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
    REAL, DIMENSION(size(u,1),size(u,2)) :: flux
    REAL, DIMENSION(size(u,1)) :: rho, rhov, v 

    rhov = u(:,2); rho = u(:,1); v=rhov/rho;
    flux(:,1) = rhov; 
    flux(:,2) = rhov*v! + rho*cssquared; 
    flux(:,3) = u(:,3)*v;  ! passive transport of rho v_y 
    flux(:,4) = u(:,4)*v;  ! passive transport of rho v_z
  END function getflux
  
  
  SUBROUTINE output_stdout(nstep,cputime,numinject,t,dt)
    INTEGER :: nstep,numinject
    REAL :: t,dt
    REAL(8) :: cputime
    WRITE(6,*),"Report from node ", node
    WRITE(6,*),node,"   nsteps = ", nstep
    WRITE(6,*),node,"   cputime = ", cputime
    WRITE(6,*),node,"   numinject = ", numinject
    WRITE(6,*),node,"   t = ", t
    WRITE(6,*),node,"   dt = ", dt
    FLUSH(6)
  END SUBROUTINE output_stdout
  
  SUBROUTINE output_file_slice(u,n,t,nstep)
    integer :: n,nstep,i
    CHARACTER*256 filename
    real, dimension(n,n,n,4) :: u
    real :: t
    
    !Write timing info
    if (node .eq. 0) then
      WRITE(filename,849) odir
      849 format(A,'output-times-slice')
      OPEN(UNIT=2, FILE=trim(filename), ACCESS='APPEND')
      WRITE(2,850) t,nstep
      850 FORMAT(E15.6,' ',I10.10)
      CLOSE(2)
    end if
    
    !Create the filename
    WRITE(filename,800) odir, nstep, node
    800 format(A,'output-slice-',I8.8,'-',I3.3)
    i = INT(n/2)
    OPEN(UNIT=2, FILE=TRIM(filename),form='unformatted',access='direct', &
         recl=sizeof(u(i,ghost+1:n-ghost,ghost+1:n-ghost,:)))
    WRITE(6,*),"Writing to: ",filename,"@ nstep=",nstep
    WRITE(2,rec=1),u(i,ghost+1:n-ghost,ghost+1:n-ghost,:)
    CLOSE(2)
  END SUBROUTINE output_file_slice
  
  SUBROUTINE output_file_cube(u,n,t,nstep)
    integer :: n,nstep
    CHARACTER*256 filename
    real, dimension(n,n,n,4) :: u
    real :: t
    
    !Write timing info
    if (node .eq. 0) then
      WRITE(filename,952) odir
      952 format(A,'output-times-cube')
      OPEN(UNIT=2, FILE=trim(filename), ACCESS='APPEND')
      WRITE(2,951) t,nstep
      951 FORMAT(E15.6,' ',I10.10)
      CLOSE(2)
    end if
    
    !Create the filename
    WRITE(filename,901) odir,nstep, node
    901 format(A,'output-cube-',I8.8,'-',I3.3)
    WRITE(6,*),"Writing to: ",filename,"@ nstep=",nstep
    OPEN(UNIT=2,FILE=TRIM(filename),form='unformatted',access='direct', &
         recl=sizeof(u(ghost+1:n-ghost,ghost+1:n-ghost,ghost+1:n-ghost,:)))
    WRITE(2,rec=1),u(ghost+1:n-ghost,ghost+1:n-ghost,ghost+1:n-ghost,:)
    CLOSE(2)
  END SUBROUTINE output_file_cube

  SUBROUTINE create_resume_file()
    CHARACTER*256 filename

    WRITE(filename,1001) odir
    1001 format(A,"resume_data")
    OPEN(UNIT=99, FILE=trim(filename), ACTION="WRITE")
    
    WRITE(99,1002) t, nstep, numinject
    1002 FORMAT(E15.6, I10.10, I10.10)

    CLOSE(99)
  END SUBROUTINE create_resume_file

  SUBROUTINE read_resume_file()
    CHARACTER*256 filename

    WRITE(filename,1001) odir
    1001 format(A,"resume_data")
    OPEN(UNIT=99, FILE=trim(filename), ACTION="READ")
    
    READ(99,1012) t, nstep, numinject
    1012 FORMAT(E15.6, I10.10, I10.10)

    CLOSE(99)
  END SUBROUTINE read_resume_file

  SUBROUTINE read_resume_cube(u,n,ghost,nstep)
    INTEGER :: n, nstep, ghost
    REAL, DIMENSION(n,n,n,4) :: u
    CHARACTER*256 filename

    WRITE(filename,1021) odir,nstep, node
    1021 format(A,'output-cube-',I8.8,'-',I3.3)
    OPEN(UNIT=99, FILE=trim(filename), FORM='unformatted', ACCESS='direct', &
         recl=sizeof(u(ghost+1:n-ghost,ghost+1:n-ghost,ghost+1:n-ghost,:)))
    READ(99,rec=1),u(ghost+1:n-ghost,ghost+1:n-ghost,ghost+1:n-ghost,:)

    CLOSE(2)
  END SUBROUTINE read_resume_cube
  
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
  
end
