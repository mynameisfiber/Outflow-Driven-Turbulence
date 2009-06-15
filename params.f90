MODULE params
  REAL, PARAMETER :: CFL = 0.65
  INTEGER, PARAMETER :: procs = 1
  INTEGER, PARAMETER :: n=108
  INTEGER, PARAMETER :: ghost=4
  INTEGER, PARAMETER :: seed=4212 !4224
  REAL, PARAMETER :: MAXVEL = 50.0
  INTEGER, PARAMETER :: bndryfreq = 5
  CHARACTER(LEN=*), PARAMETER :: odir="/home/fiber/Programming/outflow/output/"
  !CHARACTER(LEN=*), PARAMETER :: odir="/mnt/node_scratch/mgorelick/"

  CHARACTER(LEN=*), PARAMETER :: resume = ""

  !run A in Proposed Runs
  REAL, PARAMETER :: oImp = 160000.0
  REAL, PARAMETER :: oSnorm0 = 9.765e-07
  INTEGER, PARAMETER :: or = 15
  
  !run I in Proposed Runs
  !REAL, PARAMETER :: oImp = 1280000.0
  !REAL, PARAMETER :: oSnorm0 = 6.1035e-08
  !INTEGER, PARAMETER :: or = 30
   
  REAL, PARAMETER :: odinj=0.0
  REAL, PARAMETER :: osoft=0.05
  REAL, PARAMETER :: on = 0.0
  REAL, PARAMETER :: op = 1.0
 
  INTEGER, PARAMETER :: MAXWALLTIME = 171900 !~48hr
  INTEGER, PARAMETER :: MAXSTEP =0
  INTEGER, PARAMETER :: outputfreq=500
  INTEGER, PARAMETER :: MAXTIME = 100
  
  REAL, DIMENSION(3) :: snapshotfreqt = (/ 1.0/20.0, 1.0/20.0, 1.0 /)
  INTEGER, DIMENSION(3) :: snapshotfreqnstep = (/ 0,0,0 /)
END MODULE params
