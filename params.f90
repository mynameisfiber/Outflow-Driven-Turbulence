MODULE physicalparams
  INTEGER, PARAMETER :: procs = 8
  INTEGER, PARAMETER :: n=86
  INTEGER, PARAMETER :: ghost=3
  INTEGER, PARAMETER :: seed=4352
  
  ! I/S s.t. v=1.3M and l=20 and t=15.3
  REAL, PARAMETER :: oImp=10400
  REAL, PARAMETER :: oSnorm0=8.124999e-6
  
  REAL, PARAMETER :: odinj=0.05
  REAL, PARAMETER :: osoft=0.1
  REAL, PARAMETER :: colnorm=4.4918960E-04
  REAL, PARAMETER :: on = 0.0
  REAL, PARAMETER :: op = 1.0
  INTEGER, PARAMETER :: or = 6
END MODULE physicalparams