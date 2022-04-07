
! Matrix function -----------------------------------------------------------
! inner product -------------------------------------------------------------
! inner product of vectors
! args   : real*8 *a,*b     I   vector a,b (n x 1)
!          integer*4 n      I   size of vector a,b
! return : a'*b
!----------------------------------------------------------------------------
real*8 function dot(a, b, n)
implicit none
include 'file_para.h'
integer*4, intent(in) :: n
real*8, intent(in) :: a(n),b(n)
real*8 c(n)
c=a*b
dot=sum(c)
end function

! euclid norm ---------------------------------------------------------------
! euclid norm of vector
! args   : real*8 *a        I   vector a (n x 1)
!          integer*4 n      I   size of vector a
! return :  .or.  a  .or. 
!----------------------------------------------------------------------------
real*8 function norm(a, n)
implicit none
include 'file_para.h'
integer*4, intent(in) :: n
real*8, intent(in) :: a(n)
real*8, external :: dot
norm=dsqrt(dot(a,a,n))
end function

subroutine mattran(n, m, A, B)  ! A(n,m) B(m,n)
implicit none
include 'file_para.h'
integer*4, intent(in) :: n,m
real*8, intent(in) :: A(n,m)
real*8, intent(out) :: B(m,n)
B=reshape(A,(/m,n/),order=(/2,1/))
end subroutine

subroutine matmul2(tr, n, k, m, A, B, C)
implicit none
include 'file_para.h'
character(2), intent(in) :: tr
integer*4, intent(in) :: n, k, m
real*8, intent(in)  :: A(n,m),B(m,k)
real*8, intent(out) :: C(n,k)
real*8 A2(n,m),B2(m,k)
A2=A; B2=B; C=0.d0

if(tr(1:1)=='N')then
    if(tr(2:2)=='T') B2=reshape(B,(/m,k/),order=(/2,1/))
elseif(tr(1:1)=='T')then
    A2=reshape(A,(/n,m/),order=(/2,1/))
    if(tr(2:2)=='T') B2=reshape(B,(/m,k/),order=(/2,1/))
endif
C=matmul(A2,B2)
end subroutine

! gauss-jordan method to compute the inverse of matrix a
! a -- is the n*n matrix ,it's input and output(inverse of a)
! n--- is the n ,input
! l--- if the matrix a is dsingular ,l=0 ,else l /=0
! is,js-----is a matrix is(n),and js(n)
subroutine matinv(a, n, l)
implicit none
include 'file_para.h'
integer*4, intent(in) :: n
integer*4, intent(out) :: l
real*8    a(n,n),t,d
integer*4 i,j,k,is(n),js(n)
l=0
do k=1,n
    d=0.d0
    do i=k,n
        do j=k,n
            if (abs(a(i,j)).gt.d) then
                d=abs(a(i,j))
                is(k)=i
                js(k)=j
            endif
        enddo
    enddo
    if (d<=1d-20) then
        l=-1; !write(*,"(1x,'err**not inv')")
        return
    endif

    do j=1,n
        t=a(k,j)
        a(k,j)=a(is(k),j)
        a(is(k),j)=t
    enddo
    do i=1,n
        t=a(i,k)
        a(i,k)=a(i,js(k))
        a(i,js(k))=t
    enddo
    a(k,k)=1.d0/a(k,k)
    do j=1,n
        if (j.ne.k) then
            a(k,j)=a(k,j)*a(k,k)
        endif
    enddo
    do i=1,n
        if (i.ne.k) then
            do j=1,n
                if (j.ne.k) then
                    a(i,j)=a(i,j)-a(i,k)*a(k,j)
                endif
            enddo
        endif
    enddo
    do i=1,n
        if (i.ne.k) then
            a(i,k)=-a(i,k)*a(k,k)
        endif
    enddo
enddo
do k=n,1,-1
    do j=1,n
        t=a(k,j)
        a(k,j)=a(js(k),j)
        a(js(k),j)=t
    enddo
    do i=1,n
        t=a(i,k)
        a(i,k)=a(i,is(k))
        a(i,is(k))=t
    enddo
enddo
end subroutine

!! solve linear equation ----------------------------------------------------
subroutine solve(tr, A, Y, n, m, X,stat)
implicit none
include 'file_para.h'
character(*), intent(in) :: tr
integer*4, intent(in) :: n,m
real*8, intent(in) :: A(n,n),Y(n,m)
real*8, intent(out) :: X(n,m)
integer*4, intent(out) :: stat
real*8 B(n,n)
integer*4 info
B=A
call matinv(B,n,info)
if(info==0)then
    if(tr(1:1)=='N') call matmul2('NN',n,m,n,B,Y,X)
    if(tr(1:1)/='N') call matmul2('TN',n,m,n,B,Y,X)
endif
stat=info
end subroutine

! least square estimation ---------------------------------------------------
! least square estimation by solving normal equation (x=(A*A')^-1*A*y)
! args   : real*8 *A        I   transpose of (weighted) design matrix (n x m)
!          real*8 *y        I   (weighted) measurements (m x 1)
!          integer*4 n,m    I   number of parameters and measurements (n<=m)
!          real*8 *x        O   estmated parameters (n x 1)
!          real*8 *Q        O   esimated parameters covariance matrix (n x n)
! return : status (0:ok,0>:error)
! notes  : for weighted least square, replace A and y by A*w and w*y (w=W^(1/2))
!          matirix stored by column-major order (fortran convention)
!----------------------------------------------------------------------------
subroutine lsq(A, y, n, m, x, Q, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: n,m
real*8, intent(in) :: A(n,m),y(m,1)
real*8, intent(out) :: Q(n,n),x(n,1)
integer*4, intent(out) :: stat
real*8 Ay(n,1)
integer*4 info
Q=0.d0
x=0.d0
if(m<n)then
    stat=-1; return
endif
Ay=0.d0
call matmul2('NN',n,1,m,A,y,Ay)
call matmul2('NT',n,n,m,A,A,Q)
call matinv(Q,n,info)
if(info==0) call matmul2('NN',n,1,n,Q,Ay,x)
stat=info
end subroutine
