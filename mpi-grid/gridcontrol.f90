program gridcontrol
  
implicit none
include "mpif.h"

! Variable Definitions:
!     procs - number of processes
!     n - length of local node
!     globaln - length of global grid
!     ghost - number of ghost cells
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
integer, parameter :: procs = 8, n=16, ghost=3
real, parameter :: sqrt2=sqrt(2.0)
integer, DIMENSION(3) :: dims
real, DIMENSION(n,n,n,4) :: u
integer, DIMENSION(3) :: offset, coords, isedge
integer :: node, globaln
integer :: ierr, COMM_CART
character*64 :: randfile = "random.txt"
integer :: nstep=0, or = 4


!Check that all variables make sense
IF (procs**(1.0/3) .ne. int(procs**(1.0/3))) THEN
  PRINT *,"# Procs must be a perfect cube"
END IF
globaln = (procs)**(1.0/3)*n
dims = globaln/n


!Initialize MPI
call MPI_INIT(ierr)

!Create the cartesian grid and find the local ID in it
call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,(/ .TRUE.,.TRUE.,.TRUE. /), &
     .FALSE.,COMM_CART,ierr)
call MPI_CART_MAP(MPI_COMM_WORLD,3,dims,(/ .TRUE.,.TRUE.,.TRUE. /), &
     node,ierr)

!Collects the local coordinates in node space and calculates the offset
call MPI_CART_COORDS(COMM_CART,node,3,offset,ierr)
coords = offset * n
isedge = INT(2*coords/(globaln-n)-1)
print*,"node=",node,"(x,y,z) = ",coords,"isedge = ",isedge


!Now we open the random number file and set it to IO unit 1
! note: read mode is default
OPEN(UNIT=1, FILE=randfile)


!Initialize the grid
CALL setup(u,n)

do
  if (node .eq. 0) print*,"nstep=",nstep
  
  !Manage outflows
  CALL outflow_manager(u,n)
  
  !Now test out the boundry conditions
  call boundary(u,n,ghost)

  !Output
  if (node .eq. 0) call output(u,n,nstep)
  nstep = nstep + 1
  GOTO 666
end do

666 call MPI_FINALIZE(ierr)
CLOSE(1)
return


CONTAINS

  SUBROUTINE outflow_manager(u,n)
    INTEGER :: n, oi,oj,ok
    REAL, DIMENSION(n,n,n,4) :: u
    
    !Find coordinates of the outflow
    oi = ANINT(myrand()*(globaln-1)+1)
    oj = ANINT(myrand()*(globaln-1)+1)
    ok = ANINT(myrand()*(globaln-1)+1)
    
    if (node .eq. 0) print*,"Creating outflow at: ",oi,oj,ok
    !Now we check if the outflow occures in the local grid.  This is done by
    !   finding the components of the distance at the closest approach and seeing
    !   if any component is larger than the local grid width
    !NOTE: we could do a more direct checking my finding the distance at closest
    !   approach to avoid false-positive however it would add computational time.
    IF (MAXVAL(MIN(ABS(coords-(/oi,oj,ok/)+n/2.0), &
        ABS(coords-((/oi,oj,ok/)+isedge*globaln)+n/2.0)))-or &
        .LT. n/2) THEN
      CALL generate_outflow(u,n,oi,oj,ok)
      PRINT*,"Outflow is in node:",node
    END IF
    
  END SUBROUTINE outflow_manager

  SUBROUTINE generate_outflow(u,n,oi,oj,ok)
    !The following variables are defined:
    !     oi/oj/ok - global coordinates of injection site's origin
    !     ni/nj/nk - local coordinate of arbitrary point of outflow
    !     x - distance between arbitrary point and center of outflow
    INTEGER :: n,oi,oj,ok, i,j,k, ni,nj,nk
    REAL, DIMENSION(n,n,n,4) :: u
    REAL :: r
    
    !$OMP PARALLEL DO SCHEDULE(STATIC) &
    !$OMP shared(globaln,oi,oj,ok,r,u) PRIVATE(i,j,k) DEFAULT(none)
    DO k=ok-r,ok+r
      DO j=oj-r,oj+r
        DO i=oi-r,oi+r
          r = ((i-oi)**2 + (j-oj)**2 + (k-ok)**2)**.5
          IF (r .LT. or) THEN
          
            !First we normalize to coordinates to the grid as to wrap any
            !   outflows around the periodic grid
            ni = MOD(i,globaln); nj = MOD(j,globaln); nk = MOD(k,globaln)
            IF( MAXVAL(ABS(coords+n/2-(/ni,nj,nk/))) .GT. n/2 ) THEN
              
              
              !NOW CHECK IF ni/nj/nk IS IN THE LOCAL GRID!
            
              !Now we inject the outflow into the local grid
             END IF
          END IF
        END DO
      END DO
    END DO
    
  
  END SUBROUTINE generate_outflow

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
    ! if (node .eq. 0)then
    !   print*,"rd",rd
    !   print*,"ru",ru
    !   print*,"sds",sds
    !   print*,"sde",sde
    !   print*,"sus",sus
    !   print*,"sue",sue
    !   print*,"size(ru)",size(u(ru(1):n,ru(2):n,ru(3):n,:))
    !   print*,"size(rd)",size(u(1:rd(1),1:rd(2),1:rd(3),:))
    !   print*,"size(sd)",size(u(sds(1):sde(1),sds(2):sde(2),sds(3):sde(3),:))
    !   print*,"size(su)",size(u(sus(1):sue(1),sus(2):sue(2),sus(3):sue(3),:))
    !   print*,"count",count
    ! end if
    DO dim=0,2
      !Find neighboors
      call MPI_Cart_shift (COMM_CART, dim, 1, &
           ndown, tmp, ierr)
      call MPI_Cart_shift (COMM_CART, dim, -1, &
           nup, tmp, ierr)
                          
      !First we send to the down
      !print*,"D",node,"->",ndown
      call MPI_Sendrecv(u(sds(1):sde(1),sds(2):sde(2),sds(3):sde(3),:), count,&
                        MPI_REAL, ndown, node, u(ru(1):n,ru(2):n,ru(3):n,:), &
                        count, MPI_REAL, nup, nup, COMM_CART, &
                        MPI_STATUS_IGNORE, ierr)
      
      !Now we send to the up
      !print*,"U",node,"->",nup
      call MPI_Sendrecv(u(sus(1):sue(1),sus(2):sue(2),sus(3):sue(3),:), count,&
                        MPI_REAL, nup, node, u(1:rd(1),1:rd(2),1:rd(3),:), &
                        count, MPI_REAL, ndown,ndown,COMM_CART, &
                        MPI_STATUS_IGNORE,ierr)
                       
      !Adjust boundries to be sent for next iteration 
      rd = cshift(rd,1)
      ru = cshift(ru,1)
      sds = cshift(sds,1)
      sde = cshift(sde,1)
      sus = cshift(sus,1)
      sue = cshift(sue,1)
    END DO
  END SUBROUTINE boundary
  
  SUBROUTINE output(u,n,nstep)
    integer :: n,nstep,i,j,k
    CHARACTER*128 filename
    real, dimension(n,n,n,4) :: u
    
    !Create the filename
    WRITE(filename,800) nstep, node
    800 format('output-',I8.8,'-',I3.3)
    OPEN(UNIT=2, FILE=TRIM(filename))
    print*,"Writing to: ",filename,"@ nstep=",nstep
    DO i = 1,n
      DO j = 1,n
        DO k = 1, n
          write(2,900) u(i,j,k,1) , u(i,j,k,2),  u(i,j,k,3), u(i,j,k,4); 
          900 format(E15.6,' ',E15.6,' ',E15.6,' ',E15.6)
               !i,j,k, rho, rho vx, rho vy, rho vz. 
        END DO !k
      END DO !j
    END DO !i
    CLOSE(2)
  END SUBROUTINE output
  
  SUBROUTINE setup(u,n)
    INTEGER :: n
    REAL, DIMENSION(n,n,n,4) :: u
    
    !Constant uniform initial density
    u(:,:,:,1) = 1
    
    !Zero velocity field
    u(:,:,:,2) = 0
    u(:,:,:,3) = 0
    u(:,:,:,4) = 0
  END SUBROUTINE setup

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
