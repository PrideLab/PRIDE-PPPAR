!!   INPUT  : q
!!
!!   OUTPUT : u

subroutine pnormal(q, u)
implicit none

real*8 q,u
real*8 p,y,b0,bn(10)

b0 = 1.570796288d0
bn(1) = 0.3706987906d-1
bn(2) = -0.8364353598d-3
bn(3) = -0.2250947176d-3
bn(4) = 0.6841218299d-5
bn(5) = 0.5824238515d-5
bn(6) = -0.1045274970d-5
bn(7) = 0.8360937017d-7
bn(8) = -0.3231081277d-8
bn(9) = 0.3657763036d-10
bn(10) = 0.6936233982d-12

if (q .ge. 1.d0 .or. q .le. 0.d0) then
  write (*, '(a)') '***ERROR(pnormal): input q can not >=1 or <=0 '
  call exit(1)
else if (q .eq. 0.5d0) then
  u = 0.d0
  return
endif

if (q .lt. 0.5d0 .and. q .gt. 0.d0) then
  p = q
else if (q .gt. 0.5d0) then
  p = 1.d0 - q
endif

y = -1.d0*dlog(4.d0*p*(1.d0 - p))

y = y*(b0 + y*(bn(1) + y*(bn(2) + y*(bn(3) + y*(bn(4) + y*(bn(5) + y*(bn(6) + y*(bn(7) + y*(bn(8) + y*(bn(9) + y*bn(10)))))))))))

u = dsqrt(y)
if (q .gt. 0.5d0) u = -1.d0*u

return
end
