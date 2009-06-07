MODULE params
  REAL, PARAMETER :: CFL = 0.65
  INTEGER, PARAMETER :: procs = 8
  INTEGER, PARAMETER :: n=108
  INTEGER, PARAMETER :: ghost=4
  INTEGER, PARAMETER :: seed=4212 !4224
  REAL, PARAMETER :: MAXVEL = 50.0
  INTEGER, PARAMETER :: bndryfreq = 5
  CHARACTER(LEN=*), PARAMETER :: odir="/home/fiber/Programming/outflow/output/"
  !CHARACTER(LEN=*), PARAMETER :: odir="/mnt/node_scratch/mgorelick/"

  LOGICAL, PARAMETER :: resume = .False.

  ! I/S s.t. v=1.3M and l=20 and t=15.3
  !REAL, PARAMETER :: oImp=10400
  !REAL, PARAMETER :: oSnorm0=8.124999e-6
  
  ! I/S s.t. v=2.5M and l=20 and t=8
  !REAL, PARAMETER :: oImp=20000
  !REAL, PARAMETER :: oSnorm0=1.5625e-5
  
  ! testing
  !REAL, PARAMETER :: oImp=20000
  !REAL, PARAMETER :: oSnorm0=1.5625e-2

  !run A in Proposed Runs
  REAL, PARAMETER :: oImp = 160000.0
  REAL, PARAMETER :: oSnorm0 = 9.765e-07
  INTEGER, PARAMETER :: or = 15 
  
  REAL, PARAMETER :: odinj=0.0
  REAL, PARAMETER :: osoft=0.05
  REAL, PARAMETER :: on = 0.0
  REAL, PARAMETER :: op = 1.0
 
  INTEGER, PARAMETER :: MAXWALLTIME = 300
  INTEGER, PARAMETER :: MAXSTEP =0
  INTEGER, PARAMETER :: outputfreq=100
  INTEGER, PARAMETER :: MAXTIME = 100
  
  !REAL, DIMENSION(3) :: snapshotfreqt = (/ 0.0, 0.0, 0.0 /)
  !INTEGER, DIMENSION(3) :: snapshotfreqnstep = (/ 0,0,2 /)
  REAL, DIMENSION(3) :: snapshotfreqt = (/ 1.0/20.0, 1.0/20.0, 1.0 /)
  INTEGER, DIMENSION(3) :: snapshotfreqnstep = (/ 0,0,0 /)
END MODULE params
