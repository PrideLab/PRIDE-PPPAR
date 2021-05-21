!!   INPUT  : n,q,typ
!!
!!   OUTPUT : valu

subroutine pchi2(n, q, typ, valu)
implicit none

integer*4 n,typ
real*8    q,valu
!
real*8    p,w,x,x0,pp,d
integer*4 k

if (q .le. 0.d0) then
  write (*, '(a,f5.1,a)') ' ****ERROR(pchi2): input q can not <= 0'
  call exit(1)
endif

if (n .eq. 1) then
  call pnormal(q/2.d0, x)
  valu = x*x
  return
endif
if (n .eq. 2) then
  valu = -2.d0*dlog(q)
  return
endif

p = 1.d0 - q
call pnormal(q, x)
w = 2.d0/(9.d0*dble(n))
if (w .lt. 0.d0) w = 0.d0
x0 = dble(n)*(1.d0 - w + x*dsqrt(w))**3

if (typ .eq. 0) then
  valu = x0
  return
endif

k = 0

100 call chi2(n, x0, pp, d)

if (abs((d - 0.d0)) .lt. 1.d-20) then
  valu = x0
  return
endif

valu = x0 - (pp - p)/d
if (abs(x0 - valu) .le. ((1.d-4)*abs(valu))) return

k = k + 1
if (k .ge. 30) return

x0 = valu
goto 100


return
end
