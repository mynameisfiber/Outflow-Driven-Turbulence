PROGRAM test

INTEGER,DIMENSION(3) :: coords, a, b, c
INTEGER :: oi,oj,ok
n = 10

coords = (/ 0,0,0 /)
oi = 1
oj = 1
ok = 1

PRINT*,(coords - (/oi,oj,ok/) + n/2.0)
PRINT*,(coords - (/oi,oj,ok/) + n/2.0)**2
PRINT*,SUM((coords - (/oi,oj,ok/) + n/2.0)**2)
PRINT*,SQRT(SUM((coords - (/oi,oj,ok/) + n/2.0)**2))


a = (/ 1,1,1 /)
b = (/ 10,10,10 /)
c = (/ 5,7,17 /)

PRINT*,"min(a,c)",min(a,c)
PRINT*,"min(b,c)",min(b,c)

END PROGRAM test
