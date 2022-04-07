
! Mathematical function -----------------------------------------------------
real*8 function SQR(x)
implicit none
include 'file_para.h'
real*8, intent(in) :: x
SQR=((x)*(x))
end function

real*8 function SQRT2(x)
implicit none
include 'file_para.h'
real*8, intent(in) :: x
real*8, external :: rcond
SQRT2=rcond(x<=0.d0.or.x/=x,0.d0,dsqrt(x))
end function

! average of array (except 0)
real*8 function raver(x,n)
implicit none
include 'file_para.h'
real*8, intent(in) :: x(n)  !x(:)
integer*4, intent(in) :: n
integer*4 :: i, numZero
numZero=0
do i=1,n  !size(x)
    if(dabs(x(i))<=1d-20) numZero=numZero+1
enddo
raver=sum(x)/(n-numZero)
end function

! rms of array (except 0)
real*8 function rrms(x,n)
implicit none
include 'file_para.h'
real*8, intent(in) :: x(n)  !x(:)
integer*4, intent(in) :: n
real*8, allocatable :: x1(:)
real*8 aver,raver
integer*4 i, numZero
external :: raver

numZero=0
allocate(x1(n))  !(x1(size(x,dim=1)))
x1=x
do i=1,n  !size(x1,dim=1)
    if(dabs(x1(i))<=1d-20) numZero=numZero+1
enddo
aver=raver(x1,n)
do i=1,n  !size(x1,dim=1)
    if(dabs(x1(i))<=1d-20) x1(i)=aver
enddo
x1=(x1-aver)*(x1-aver)
rrms=dsqrt(sum(x1)/(n-numZero))
deallocate(x1)
end function

! delete values that exceed 3σ
subroutine errorfilter(a, b, n)  ! (except 0)
implicit none
include 'file_para.h'
integer*4, intent(in) :: n
real*8, intent(in) :: a(n)  !a(:)
real*8, intent(out) :: b(n) !b(:)
integer*4 i
real*8 aver, rms, raver, rrms
external :: raver, rrms
aver=raver(a,n)
rms=rrms(a,n)
b=a
do i=1,n  !size(a,dim=1)
    if(dabs(b(i)-aver)>=rms*3d0) b(i)=0
enddo
end subroutine

integer*4 function ROUND(x)
implicit none
include 'file_para.h'
real*8, intent(in) :: x
ROUND=int(floor(x+0.5d0))
end function

real*8 function ROUNDF(x, n)
implicit none
include 'file_para.h'
real*8, intent(in) :: x
integer*4, intent(in) :: n
integer*4 n1, tmp
integer*4, external :: icond, ROUND
n1=icond(n<=0,0,n)
ROUNDF=ROUND(x*10**n)*1.d0/(10**n)
end function

! least common multiple (int) -----------------------------------------------
integer*4 function iminmul(a, b)
implicit none
include 'file_para.h'
integer*4, intent(in) :: a, b
integer*4 a1, b1, product, tmp
product=a*b
a1=max(a,b)
b1=min(a,b)
do while (mod(a1,b1)/=0)
    tmp=b1
    b1=mod(a1,b1)
    a1=tmp
enddo
iminmul=product/b1
end function
