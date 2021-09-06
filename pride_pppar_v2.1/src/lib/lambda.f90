!! reference to RTKLIB
!! rewritten by Jihang Lin 2020

! lambda/mlambda integer least-square estimation -------------------------------
! integer least-square estimation. reduction is performed by lambda (ref.[1]),
! and search by mlambda (ref.[2]).
! args   : integer    n      I  number of float parameters
!          integer    m      I  number of fixed solutions!
!          real*8     a      I  float parameters (n x 1)
!          real*8     Q      I  covariance matrix of float parameters (n x n)
!          real*8     F      O  fixed solutions (n x m)
!          real*8     s      O  sum of squared residulas of fixed solutions (1 x m)
!          integer    info   O  status (0:ok,other:error)
! notes  : matrix stored by column-major order (fortran convension)
! ------------------------------------------------------------------------------
subroutine Lambda(n, m, a, Q, F, s, info)
    IMPLICIT NONE
    dimension a(n), za(n), s(m), E(n*m), F(n*m)
    dimension D(n), Q(n*n), L(n*n), Z(n*n), Zt(n*n), Zi(n*n)
    integer n, m
    real*8  a, za, s, E, F
    real*8  D, Q, L, Z, Zt, Zi
    integer info
    
    info = 0
    if (n <= 0.or.m <= 0) then
        info = -1
        return
    end if
    
    call lam_Eye(Z, n)
    
    ! LD factorization
    call LD(n, Q, L, D, info)
    if (info == 0) then
        ! lambda reduction
        call Reduction(n, L, D, Z)
        call lam_MatTrans(Z, n, n, Zt)
        call lam_MatMult(Zt, a, n, 1, n, za)
        ! mlambda search
        call Search(n, m, L, D, za, E, s, info)
        if (info == 0) then
            Zi = Zt
            call lam_MatInv(Zi, n, info)
            call lam_MatMult(Zi, E, n, m, n, F)
        end if
    end if
end subroutine Lambda

! LD factorize, i.e. Q = L' * diag(D) * L --------------------------------------
subroutine LD(n, Q, L, D, info)
    IMPLICIT NONE
    integer n, i, j, k
    real*8  t
    real*8, dimension(n) :: D
    real*8, dimension(n*n) :: Q, L, A
    integer info
    
    info = 0
    A = Q

    do i = n, 1, -1
        D(i) = A(i+(i-1)*n)
        if (D(i) < 0.0D0) then
            info = -1
            exit
        end if
        t = SQRT(D(i))
        do j = 1, i
            L(i+(j-1)*n) = A(i+(j-1)*n) / t
        end do
        do j = 1, i-1, +1
            do k = 1, j
                A(j+(k-1)*n) = A(j+(k-1)*n) - L(i+(k-1)*n) * L(i+(j-1)*n)
            end do
        end do
        do j = 1, i
            L(i+(j-1)*n) = L(i+(j-1)*n) / L(i+(i-1)*n)
        end do
    end do
    
    if (info /= 0) then
        Print *, "LD factorization error\n"
    else 
        info = 0
    end if
end subroutine LD

! lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ------------------
subroutine Reduction(n, L, D, Z)
    IMPLICIT NONE
    integer n, i, j, k
    real*8  del
    real*8, dimension(n) :: D
    real*8, dimension(n*n) :: L, Z
    j = n - 1
    k = n - 1
    do while(j > 0)
        if (j <= k) then
            do i = j+1, n, +1
                call lam_Gauss(n, L, Z, i, j)
            end do
        end if
        del = D(j) + L(j+1+(j-1)*n)**2 * D(j+1)
        if (del + 1E-6 < D(j+1)) then
        ! compared considering numerical error
            call lam_Perm(n, L, D, j, del, Z)
            k = j
            j = n - 1
        else 
            j = j - 1
        end if
    end do
end subroutine Reduction

! modified lambda (mlambda) search (ref. [2]) ----------------------------------
subroutine Search(n, m, L, D, zs, zn, s, info)
    IMPLICIT NONE
    integer, parameter :: LOOPMAX = 10000 ! maximum count of search loop
    integer n, m, i, j, k, c, nn, imax
    real*8  Sgn, RoundUp
    real*8  newdist, maxdist, y
    real*8, dimension(n) :: D, dist, z, zb, zs, step
    real*8, dimension(m) :: s
    real*8, dimension(n*m) :: zn
    real*8, dimension(n*n) :: L, T
    integer info
    
    T = 0.0D0
    dist = 0.0D0
    k = n
    nn = 1
    imax = 1
    zb(k) = zs(k)
    z(k)  = RoundUp(zb(k))
    y = zb(k) - z(k)
    step(k) = Sgn(y)
    
    maxdist = 1E38
    do c = 0, LOOPMAX
        newdist = dist(k) + y**2 / D(k)
        if (newdist < maxdist) then
            if (k /= 1) then
                k = k - 1
                dist(k) = newdist
                do i = 1, k
                    T(k+(i-1)*n) = T(k+1+(i-1)*n) &
                        + (z(k+1) - zb(k+1)) * L(k+1+(i-1)*n)
                end do
                zb(k) = zs(k) + T(k+(k-1)*n)
                z(k)  = RoundUp(zb(k))
                y = zb(k) - z(k)
                step(k) = SGN(y)
            else
                if (nn <= m) then
                    if (nn == 1.or.newdist > s(imax)) then
                        imax = nn
                    end if
                    do i = 1, n
                        zn(i+(nn-1)*n) = z(i)
                    end do
                    s(nn) = newdist
                    nn = nn + 1
                else
                    if (newdist < s(imax)) then
                        do i = 1, n
                            zn(i+(imax-1)*n) = z(i)
                        end do
                        s(imax) = newdist
                        imax = 1
                        do i = 1, m
                            if(s(imax) < s(i)) imax = i
                        end do
                    end if
                    maxdist = s(imax)
                end if
                z(1) = z(1) + step(1)
                y = zb(1) - z(1)
                step(1) = -step(1) - Sgn(step(1))
            end if
        else
            if (k == n) then
                exit
            else
                k = k + 1
                z(k) = z(k) + step(k)
                y = zb(k) - z(k)
                step(k) = -step(k) - Sgn(step(k))
            end if
        end if
    end do
    
    do i = 1, m
    ! sort by s
        do j = i+1, m-1, +1
            if (s(i) < s(j)) then
                cycle
            end if
            call lam_Swap(s(i), s(j))
            do k = 1, n
                call lam_Swap(zn(k+(i-1)*n), zn(k+(j-1)*n))
            end do
        end do
    end do
    
    if (c >= LOOPMAX) then
        Print *, "search loop count overflow\n"
        info = -1
    else
        info = 0
    end if
end subroutine Search

real*8 function Sgn(x)
    IMPLICIT NONE
    real*8 x
    if (x <= 0) then
        Sgn = -1.0D0
    else
        Sgn = +1.0D0
    endif
end function Sgn

real*8 function RoundUp(x)
    IMPLICIT NONE
    real*8 x
    RoundUp = Floor(x + 0.5D0)
end function RoundUp

subroutine lam_Swap(x, y)
    IMPLICIT NONE
    real*8 x, y, t
    t = x
    x = y
    y = t
end subroutine lam_Swap


! Integer lam_Gauss transformation ------------------------------------------------- 
subroutine lam_Gauss(n, L, Z, i, j)
    IMPLICIT NONE
    integer n, i, j, k, mu
    real*8  RoundUp
    real*8, dimension(n*n) :: L, Z
    mu = NINT(RoundUp(L(i+(j-1)*n)))
    if (mu /= 0) then
        do k = i, n
            L(k+(j-1)*n) = L(k+(j-1)*n) - mu * L(k+(i-1)*n)
        end do
        do k = 1, n
            Z(k+(j-1)*n) = Z(k+(j-1)*n) - mu * Z(k+(i-1)*n)
        end do
    end if
end subroutine lam_Gauss

! lam_Permutation ------------------------------------------------------------------
subroutine lam_Perm(n, L, D, j, del, Z)
    IMPLICIT NONE
    integer n, j, k
    real*8  del, eta, lam, a0, a1
    real*8, dimension(n) :: D
    real*8, dimension(n*n) :: L, Z
    eta = D(j) / del
    lam = D(j+1) * L(j+1+(j-1)*n) / del
    D(j) = eta * D(j+1)
    D(j+1) = del
    do k = 1, j-1, +1
        a0 = L(j+(k-1)*n)
        a1 = L(j+1+(k-1)*n)
        L(j+(k-1)*n) = -L(j+1+(j-1)*n)*a0 + a1
        L(j+1+(k-1)*n) = eta * a0 + lam * a1
    end do 
    L(j+1+(j-1)*n) = lam
    do k = j+2, n, +1
        call lam_Swap(L(k+(j-1)*n), L(k+j*n))
    end do
    do k = 1, n
        call lam_Swap(Z(k+(j-1)*n), Z(k+j*n))
    end do
end subroutine lam_Perm



! lambda reduction -------------------------------------------------------------
! reduction by lambda (ref [1]) for integer least square
! args   : integer    n      I  number of float parameters
!          real*8     Q      I  covariance matrix of float parameters (n x n)
!          real*8     Z      O  lambda reduction matrix (n x n)
!          integer    info   O  status (0:ok,other:error)
! ------------------------------------------------------------------------------
subroutine Lambda_Reduction(n, Q, Z, info)
    IMPLICIT NONE
    dimension D(n), Q(n*n), L(n*n), Z(n*n)
    integer n, i, j, info
    real*8  D, Q, L, Z
    
    info = 0
    if (n <= 0) then
        info = -1
    end if
    
    call lam_Eye(Z, n)
    call LD(n, Q, L, D, info)
    if (info /= 0) return
    call Reduction(n, L, D, Z)
    info = 0
end subroutine Lambda_Reduction

! mlambda search ---------------------------------------------------------------
! integer least-square estimation. reduction is performed by lambda (ref.[1]),
! and search by mlambda (ref.[2]).
! args   : integer    n         I  number of float parameters
!          integer    m         I  number of fixed solutions!
!          real*8     a         I  float parameters (n x 1)
!          real*8     Q         I  covariance matrix of float parameters (n x n)
!          real*8     F         O  fixed solutions (n x m)
!          real*8     s         O  sum of squared residulas of fixed solutions (1 x m)
!          integer    info      O  status (0:ok,other:error)
! notes  : matrix stored by column-major order (fortran convension)
! ------------------------------------------------------------------------------
subroutine Lambda_Search(n, m, a, Q, F, s, info)
    IMPLICIT NONE
    dimension a(n), za(n), s(m), F(n*m)
    dimension D(n), Q(n*n), L(n*n)
    integer n, m
    real*8  a, za, s, F
    real*8  D, Q, L
    integer info
    
    info = 0
    if (n <= 0.or.m <= 0) then
        info = -1
    end if
    
    L = 0
    call LD(n, Q, L, D, info)
    if (info /= 0) return
    call Search(n, m, L, D, a, F, s, info)
end subroutine Lambda_Search

! identity matrix --------------------------------------------------------------
! turns a given matrix into an identity matrix
! args   : integer    n         I  number of rows and columns of matrix
!          real*8     Mat       O  matrix Mat(n*n) 
! (if n<=0, does nothing)
! ------------------------------------------------------------------------------
subroutine lam_Eye(Mat, n)
    IMPLICIT NONE
    integer n, i
    real*8  Mat(n*n)
    Mat = 0.D0
    forall (i = 1:n) Mat(i+(i-1)*n) = 1.0D0
end subroutine lam_Eye

! multiply matrix --------------------------------------------------------------
! multiply matrix by matrix (Mat3=Mat1*Mat2)
! args   : integer    n,k,m     I  size of matrix Mat1, Mat2
!          real*8     Mat1,Mat2 I  matrix Mat1 (n x m), Mat2 (m x k)
!          real*8     Mat3      O  matrix Mat3 (n x k)
! ------------------------------------------------------------------------------
subroutine lam_MatMult(Mat1, Mat2, n, k, m, Mat3)
    IMPLICIT NONE
    integer n, m, k
    integer i, j, t
    real*8  Mat1(n*m), Mat2(m*k), Mat3(n*k)
    do i = 1, n
        do j = 1, k
            Mat3(i+(j-1)*n) = 0.0D0
            do t = 1, m
                Mat3(i+(j-1)*n) = Mat3(i+(j-1)*n) &
                    + Mat1(i+(t-1)*n) * Mat2(t+(j-1)*m)
            end do
        end do
    end do
end subroutine lam_MatMult

! transform matrix -------------------------------------------------------------
! transform matrix (Mat2=Mat1')
! args   : integer    n,m       I  size of matrix Mat1
!          real*8     Mat1      I  matrix Mat1 (n x m)
!          real*8     Mat2      O  matrix Mat2 (m x n)
! ------------------------------------------------------------------------------
subroutine lam_MatTrans(Mat1, n, m, Mat2)
    IMPLICIT NONE
    integer n, m
    integer i, j
    real*8  Mat1(n*m), Mat2(m*n)
    do i = 1, n
        dO j = 1, m
            Mat2(j+(i-1)*n) = Mat1(i+(j-1)*m)
        end do      
    end do
end subroutine lam_MatTrans

! inverse matrix -------------------------------------------------------------
! inverse matrix (Mat2=inv(Mat1))
! args   : integer    n         I  number of rows and columns of matrix
!          real*8     Mat1      I  matrix Mat1 (n x n)
!          real*8     Mat2      O  matrix Mat2 (n x n)
!          integer    info      O  status (0:ok,other:error)
! ------------------------------------------------------------------------------
subroutine lam_MatInv(a, n, info)
    IMPLICIT NONE
    dimension a(n*n), is(n), js(n)
    integer n, i, j, k, info
    integer is, js
    real*8  a, t, d
    
    info = 1
    
    do 100 k = 1, n
        d = 0.0D0
        do 10 i = k, n
            do 10 j = k, n
                if (abs(a(i+(j-1)*n)) > d) then
                    d = abs(a(i+(j-1)*n))
                    is(k) = i
                    js(k) = j
                end if
   10   continue
        if (d + 1.0D0 == 1.0D0) then
            info = 0
            write(*, 20)
            return
        end if
   20   format(1x,'err**not inv')
        do 30 j = 1, n
            t = a(k+(j-1)*n)
            a(k+(j-1)*n) = a(is(k)+(j-1)*n)
            a(is(k)+(j-1)*n) = t
   30   continue
        do 40 i = 1, n
            t = a(i+(k-1)*n)
            a(i+(k-1)*n) = a(i+(js(k)-1)*n)
            a(i+(js(k)-1)*n) = t
   40   continue
        a(k+(k-1)*n) = 1 / a(k+(k-1)*n)
        do 50 j = 1, n
            if (j /= k) then
                a(k+(j-1)*n) = a(k+(j-1)*n) * a(k+(k-1)*n)
            end if
   50   continue
        do 70 i = 1, n
            if (i /= k) then
                do 60 j = 1, n
                    if (j /= k) then
                        a(i+(j-1)*n) = a(i+(j-1)*n) &
                            - a(i+(k-1)*n) * a(k+(j-1)*n)
                    end if
   60           continue
            end if
   70   continue
        do 80 i = 1, n
            if (i /= k) then
                a(i+(k-1)*n) = -a(i+(k-1)*n) * a(k+(k-1)*n)
            end if
   80   continue
100 continue
    do 130 k = n, 1, -1
        do 110 j = 1, n
            t = a(k+(j-1)*n)
            a(k+(j-1)*n) = a(js(k)+(j-1)*n)
            a(js(k)+(j-1)*n) = t
110     continue
        do 120 i = 1, n
            t = a(i+(k-1)*n)
            a(i+(k-1)*n) = a(i+(is(k)-1)*n)
            a(i+(is(k)-1)*n) = t
120     continue
130 continue
    return
end subroutine lam_MatInv
