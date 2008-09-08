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
!     r - radius of outflow
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
integer :: r = 4, oi,oj,ok, nstep=0
logical :: islocal = .FALSE.

!Check that all variables make sense
IF (procs**(1.0/3) .ne. int(procs**(1.0/3))) THEN
  PRINT *,"# Procs must be a perfect cube"
END IF
globaln = (procs)**(1.0/3)*n
dims = globaln/n

!Now we open the random number file and set it to IO unit 1
! note: read mode is default
OPEN(UNIT=1, FILE=randfile)

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

!Initialize the grid
u = node

!MAIN LOOP
do
  if (node .eq. 0) print*,"nstep=",nstep
  
  !Find coordinates of the outflow
  oi = ANINT(myrand()*(globaln-1)+1)
  oj = ANINT(myrand()*(globaln-1)+1)
  ok = ANINT(myrand()*(globaln-1)+1)
  
  if (node .eq. 0) print*,"Creating outflow at: ",oi,oj,ok
  !Now we check if the outflow occures in the local grid.  This is done by
  !   finding the minimum distance between the center of the grid and the
  !   outflow (using isedge to relocate the outflow based on periodicity)
  islocal = .FALSE.
  IF ( SQRT( SUM( MIN( ABS(coords-(/oi,oj,ok/)+n/2.0), &
      ABS(coords-((/oi,oj,ok/)+isedge*globaln)+n/2.0) )**2 ) ) -r .lt. &
      sqrt2*(n-.5)/2) THEN
    islocal = .TRUE.
  END IF
  IF (islocal) PRINT*,"Outflow is in node:",node
  
  
  !Now test out the boundry conditions
  call boundary(u,n,ghost,node)

  !Output
  if (node .eq. 0) call output(u,n,nstep,node)
  nstep = nstep + 1
end do

666 call MPI_FINALIZE(ierr)
CLOSE(1)
return


CONTAINS

  SUBROUTINE boundary(u,n,ghost,node)
    integer :: n, node, ghost, dim, ndown, nup, count, ierr, tmp
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
  
  SUBROUTINE output(u,n,nstep,node)
    integer :: n,nstep,node,i,j,k
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
