PROGRAM outflowlist
type node
    type(node), pointer :: next
    real :: num
end type node
INTEGER :: i
type(node), pointer :: ll, cur, first
call srand(234)

PRINT*,"45%40=",MODULO(45,40)
PRINT*,"-45%40=",MODULO(-45,40)

allocate(ll)
cur => ll
first => cur

do i=1,10
  cur%num = i
  allocate(cur%next)
  cur => cur%next
end do
cur%num = 10
allocate(cur%next)
cur => cur%next
cur%num = 10
allocate(cur%next)
cur => cur%next



cur => first
do while (associated(cur%next))
  print*,cur%num
  do while (associated(cur%next) .and. cur%num .eq. cur%next%num)
    ll => cur%next
    cur => cur%next%next
    deallocate(ll)
    ll => cur
  end do
  cur=>cur%next
  deallocate(ll)
  ll => cur
end do

END PROGRAM outflowlist