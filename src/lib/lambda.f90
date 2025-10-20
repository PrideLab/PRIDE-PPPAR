!! reference to RTKLIB (revision 1.1)
!!     [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
!!         a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82,
!!         1995
!!     [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
!!         integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
!!
!! rewritten by Jihang Lin in 2020, modified in 2021

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
subroutine lambda(n, m, a, Q, F, s, add, info)
    implicit none

    dimension a(n), za(n), s(m), E(n*m), F(n*m)
    dimension D(n), Q(n*n), L(n*n), Z(n*n), Zt(n*n), Zi(n*n)

    integer i, j, n, m
    real*8  a, za, s, E, F, add(2)
    real*8  D, Q, L, Z, Zt, Zi
    integer info
    
    info = 0
    if (n .le. 0 .or. m .le. 0) then
        info = -1
        return
    end if
    
    call lam_eye(Z, n)
    
    call LD(n, Q, L, D, info)                      ! LD factorization


    if (info .eq. 0) then
        call reduction(n, L, D, Z)                 ! lambda reduction
        call lam_mattrans(Z, n, n, Zt)
        call lam_matmult(Zt, a, n, 1, n, za)
        call search(n, m, L, D, za, E, s, add, info)    ! mlambda search
        if (info .eq. 0) then
            Zi = Zt
            call lam_matinv(Zi, n, info)
            call lam_matmult(Zi, E, n, m, n, F)
        end if
    end if
end subroutine lambda

! LD factorize, i.e. Q = L' * diag(D) * L --------------------------------------
subroutine LD(n, Q, L, D, info)
    implicit none

    integer n, i, j, k
    real*8  t
    real*8, dimension(n) :: D
    real*8, dimension(n*n) :: Q, L, A
    integer info

    info = 0
    A = Q

    do i = n, 1, -1
        D(i) = A(i+(i-1)*n)
        if (D(i) .lt. 0.0D0) then
            info = -1
            exit
        end if
        t = SQRT(D(i))
        do j = 1, i
            L(i+(j-1)*n) = A(i+(j-1)*n) / t
        end do
        do j = 1, i-1
            do k = 1, j
                A(j+(k-1)*n) = A(j+(k-1)*n) - L(i+(k-1)*n) * L(i+(j-1)*n)
            end do
        end do
        do j = 1, i
            L(i+(j-1)*n) = L(i+(j-1)*n) / L(i+(i-1)*n)
        end do
    end do
    
    if (info .ne. 0) then
        write (*, '(a)') '***ERROR(LD): LD factorization error '
    else 
        info = 0
    end if
end subroutine LD

! lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ------------------
subroutine reduction(n, L, D, Z)
    implicit none

    integer n, i, j, k
    real*8  del
    real*8, dimension(n) :: D
    real*8, dimension(n*n) :: L, Z

    j = n - 1
    k = n - 1
    do while(j .gt. 0)
        if (j .le. k) then
            do i = j+1, n
                call lam_gauss(n, L, Z, i, j)
            end do
        end if
        del = D(j) + L(j+1+(j-1)*n)**2 * D(j+1)
        if (del + 1E-6 .lt. D(j+1)) then           ! compared considering numerical error
            call lam_perm(n, L, D, j, del, Z)
            k = j
            j = n - 1
        else 
            j = j - 1
        end if
    end do
end subroutine reduction

! modified lambda (mlambda) search (ref. [2]) ----------------------------------
subroutine search(n, m, L, D, zs, zn, s, add, info)
    implicit none

    integer, parameter :: LOOPMAX = 1000000          ! maximum count of search loop
    integer n, m, i, j, k, c, nn, imax
    real*8  sgn, roundup, add(2)
    real*8  newdist, maxdist, y
    real*8, dimension(n) :: D, dist, z, zb, zs, step, z2
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
    z(k) = roundup(zb(k))
    y = zb(k) - z(k)
    step(k) = sgn(y)
    
    maxdist = 1E38
    do c = 0, LOOPMAX
        newdist = dist(k) + y**2 / D(k)
        if (newdist .lt. maxdist) then
            if (k .ne. 1) then
                k = k - 1
                dist(k) = newdist
                do i = 1, k
                    T(k+(i-1)*n) = T(k+1+(i-1)*n) &
                        + (z(k+1) - zb(k+1)) * L(k+1+(i-1)*n)
                end do
                zb(k) = zs(k) + T(k+(k-1)*n)
                z(k)  = roundup(zb(k))
                y = zb(k) - z(k)
                step(k) = SGN(y)
            else
                if (nn .le. m) then
                    if (nn .eq. 1 .or. newdist .gt. s(imax)) then
                        imax = nn
                    end if
                    do i = 1, n
                        zn(i+(nn-1)*n) = z(i)
                    end do
                    s(nn) = newdist
                    nn = nn + 1
                else
                    if (newdist .lt. s(imax)) then
                        do i = 1, n
                            zn(i+(imax-1)*n) = z(i)
                        end do
                        s(imax) = newdist
                        imax = 1
                        do i = 1, m
                            if(s(imax) .lt. s(i)) imax = i
                        end do
                    end if
                    maxdist = s(imax)
                end if
                z(1) = z(1) + step(1)
                y = zb(1) - z(1)
                step(1) = -step(1) - sgn(step(1))
            end if
        else
            if (k .eq. n) then
                exit
            else
                k = k + 1
                z(k) = z(k) + step(k)
                y = zb(k) - z(k)
                step(k) = -step(k) - sgn(step(k))
            end if
        end if
    end do
    ! sort by s
    do i = 1, m
        do j = i, m
            if (s(i) .lt. s(j)) then
                cycle
            end if
            call lam_swap(s(i), s(j))
            do k = 1, n
                call lam_swap(zn(k+(i-1)*n), zn(k+(j-1)*n))
            end do
        end do
    end do

    add = 0.d0

    T = 0.d0
    k = n
    z(k) = zn(k) - zn(k + n)
    zb(k) = zn(k + n)
    do k = n - 1, 1, -1
        do i = 1, k
            T(k+(i-1)*n) = T(k+1+(i-1)*n) &
            + (zn(k+1) - zb(k+1)) * L(k+1+(i-1)*n)
         enddo

         zb(k) = zn(k + n) + T(k + (k - 1) * n)
         z(k) = zb(k) - zn(k)
    enddo

    do i = 1, n
        add(1) = add(1) +z(i)**2 / D(i)
    enddo
 

    T = 0.d0
    k = n
    z2(k) = zn(k) - zs(k)
    zb(k) = zs(k)
    do k = n - 1, 1, -1
        do i = 1, k
            T(k+(i-1)*n) = T(k+1+(i-1)*n) &
            + (zn(k + 1) - zb(k+1)) * L(k+1+(i-1)*n)
         enddo

         zb(k) = zs(k) + T(k + (k - 1) * n)
         z2(k) = zb(k) - zn(k)
    enddo
    do i = 1, n
        add(2) = add(2) + dabs(z(i)*z2(i)) / D(i)
    enddo

    if (c .ge. LOOPMAX) then
        write (*, '(a)') '***ERROR(search): search loop count overflow '
        info = -1
    else
        info =  0
    end if
end subroutine search

real*8 function sgn(x)
    implicit none
    real*8 x
    if (x .le. 0) then
        sgn = -1.0D0
    else
        sgn = +1.0D0
    endif
end function sgn

real*8 function roundup(x)
    implicit none
    real*8 x
    roundup = floor(x + 0.5D0)
end function roundup

subroutine lam_swap(x, y)
    implicit none
    real*8 x, y, t
    t = x
    x = y
    y = t
end subroutine lam_swap

! Integer lam_gauss transformation ------------------------------------------------- 
subroutine lam_gauss(n, L, Z, i, j)
    implicit none

    integer n, i, j, k, mu
    real*8, dimension(n*n) :: L, Z
    real*8  roundup

    mu = NINT(roundup(L(i+(j-1)*n)))
    if (mu .ne. 0) then
        do k = i, n
            L(k+(j-1)*n) = L(k+(j-1)*n) - mu * L(k+(i-1)*n)
        end do
        do k = 1, n
            Z(k+(j-1)*n) = Z(k+(j-1)*n) - mu * Z(k+(i-1)*n)
        end do
    end if
end subroutine lam_gauss

! lam_permutation ------------------------------------------------------------------
subroutine lam_perm(n, L, D, j, del, Z)
    implicit none

    integer n, j, k
    real*8  del, eta, lam, a0, a1
    real*8, dimension(n) :: D
    real*8, dimension(n*n) :: L, Z

    eta = D(j) / del
    lam = D(j+1) * L(j+1+(j-1)*n) / del
    D(j) = eta * D(j+1)
    D(j+1) = del
    do k = 1, j-1
        a0 = L(j+(k-1)*n)
        a1 = L(j+1+(k-1)*n)
        L(j+(k-1)*n) = -L(j+1+(j-1)*n)*a0 + a1
        L(j+1+(k-1)*n) = eta * a0 + lam * a1
    end do 
    L(j+1+(j-1)*n) = lam
    do k = j+2, n
        call lam_swap(L(k+(j-1)*n), L(k+j*n))
    end do
    do k = 1, n
        call lam_swap(Z(k+(j-1)*n), Z(k+j*n))
    end do
end subroutine lam_perm

! lambda reduction -------------------------------------------------------------
! reduction by lambda (ref [1]) for integer least square
! args   : integer    n      I  number of float parameters
!          real*8     Q      I  covariance matrix of float parameters (n x n)
!          real*8     Z      O  lambda reduction matrix (n x n)
!          integer    info   O  status (0:ok,other:error)
! ------------------------------------------------------------------------------
subroutine lambda_reduction(n, Q, Z, info)
    implicit none

    dimension D(n), Q(n*n), L(n*n), Z(n*n)
    integer n, i, j, info
    real*8  D, Q, L, Z
    
    info = 0
    if (n .le. 0) then
        info = -1
    end if
    
    call lam_eye(Z, n)
    call LD(n, Q, L, D, info)
    if (info .ne. 0) return
    call reduction(n, L, D, Z)
    info = 0
end subroutine lambda_reduction

! mlambda search ---------------------------------------------------------------
! integer least-square estimation. reduction is performed by lambda (ref.[1]),
! and search by mlambda (ref.[2]).
! args   : integer    n         I  number of float parameters
!          integer    m         I  number of fixed solutions
!          real*8     a         I  float parameters (n x 1)
!          real*8     Q         I  covariance matrix of float parameters (n x n)
!          real*8     F         O  fixed solutions (n x m)
!          real*8     s         O  sum of squared residulas of fixed solutions (1 x m)
!          integer    info      O  status (0:ok,other:error)
! notes  : matrix stored by column-major order (fortran convension)
! ------------------------------------------------------------------------------
subroutine lambda_search(n, m, a, Q, F, s, info)
    implicit none

    dimension a(n), za(n), s(m), F(n*m)
    dimension D(n), Q(n*n), L(n*n)

    integer n, m
    real*8  a, za, s, F
    real*8  D, Q, L, add(2)
    integer info
    
    info = 0
    if (n .le. 0 .or. m .le. 0) then
        info = -1
    end if
    
    L = 0
    call LD(n, Q, L, D, info)
    if (info .ne. 0) return
    call search(n, m, L, D, a, F, s, add, info)
end subroutine lambda_search

! identity matrix --------------------------------------------------------------
! turns a given matrix into an identity matrix
! args   : integer    n         I  number of rows and columns of matrix
!          real*8     Mat       O  matrix Mat (n x n)
! (if n.le.0, does nothing)
! ------------------------------------------------------------------------------
subroutine lam_eye(Mat, n)
    implicit none

    integer n, i
    real*8  Mat(n*n)

    Mat = 0.D0
    forall (i = 1:n) Mat(i+(i-1)*n) = 1.0D0
end subroutine lam_eye

! multiply matrix --------------------------------------------------------------
! multiply matrix by matrix (Mat3=Mat1*Mat2)
! args   : integer    n,k,m     I  size of matrix Mat1, Mat2
!          real*8     Mat1,Mat2 I  matrix Mat1 (n x m), Mat2 (m x k)
!          real*8     Mat3      O  matrix Mat3 (n x k)
! ------------------------------------------------------------------------------
subroutine lam_matmult(Mat1, Mat2, n, k, m, Mat3)
    implicit none

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
end subroutine lam_matmult

! transform matrix -------------------------------------------------------------
! transform matrix (Mat2=Mat1')
! args   : integer    n,m       I  size of matrix Mat1
!          real*8     Mat1      I  matrix Mat1 (n x m)
!          real*8     Mat2      O  matrix Mat2 (m x n)
! ------------------------------------------------------------------------------
subroutine lam_mattrans(Mat1, n, m, Mat2)
    implicit none

    integer n, m
    integer i, j
    real*8  Mat1(n*m), Mat2(m*n)

    do i = 1, n
        dO j = 1, m
            Mat2(j+(i-1)*n) = Mat1(i+(j-1)*m)
        end do      
    end do
end subroutine lam_mattrans

! inverse matrix -------------------------------------------------------------
! inverse matrix (Mat2=inv(Mat1))
! args   : integer    n         I  number of rows and columns of matrix
!          real*8     Mat1      I  matrix Mat1 (n x n)
!          real*8     Mat2      O  matrix Mat2 (n x n)
!          integer    info      O  status (0:ok,other:error)
! ------------------------------------------------------------------------------
subroutine lam_matinv(a, n, info)
    implicit none

    dimension a(n*n), is(n), js(n)

    integer n, i, j, k, info
    integer is, js
    real*8  a, t, d
    
    info = 1
    
    do k = 1, n
        d = 0.0D0
        do i = k, n
            do j = k, n
                if (abs(a(i+(j-1)*n)) .gt. d) then
                    d = abs(a(i+(j-1)*n))
                    is(k) = i
                    js(k) = j
                end if
            end do
        end do
        if (d + 1.0D0 .eq. 1.0D0) then
            info = 0
            write (*, '(a)') '***ERROR(lam_matinv): matrix singularity '
            return
        end if
        do j = 1, n
            t = a(k+(j-1)*n)
            a(k+(j-1)*n) = a(is(k)+(j-1)*n)
            a(is(k)+(j-1)*n) = t
        end do
        do i = 1, n
            t = a(i+(k-1)*n)
            a(i+(k-1)*n) = a(i+(js(k)-1)*n)
            a(i+(js(k)-1)*n) = t
        end do
        a(k+(k-1)*n) = 1 / a(k+(k-1)*n)
        do j = 1, n
            if (j .ne. k) then
                a(k+(j-1)*n) = a(k+(j-1)*n) * a(k+(k-1)*n)
            end if
        end do
        do i = 1, n
            if (i .ne. k) then
                do j = 1, n
                    if (j .ne. k) then
                        a(i+(j-1)*n) = a(i+(j-1)*n) - a(i+(k-1)*n) * a(k+(j-1)*n)
                    end if
                end do
            end if
        end do
        do i = 1, n
            if (i .ne. k) then
                a(i+(k-1)*n) = -a(i+(k-1)*n) * a(k+(k-1)*n)
            end if
        end do
    end do
    do k = n, 1, -1
        do j = 1, n
            t = a(k+(j-1)*n)
            a(k+(j-1)*n) = a(js(k)+(j-1)*n)
            a(js(k)+(j-1)*n) = t
        end do
        do i = 1, n
            t = a(i+(k-1)*n)
            a(i+(k-1)*n) = a(i+(is(k)-1)*n)
            a(i+(is(k)-1)*n) = t
        end do
    end do
    return
end subroutine lam_matinv
