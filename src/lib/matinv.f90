!!
!! matinv.f90
!!
!!    Copyright (C) 2023 by Wuhan University
!!
!!    This program belongs to PRIDE PPP-AR which is an open source software:
!!    you can redistribute it and/or modify it under the terms of the GNU
!!    General Public License (version 3) as published by the Free Software Foundation.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!    GNU General Public License (version 3) for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program. If not, see <https://www.gnu.org/licenses/>.
!!
!! Contributor: Maorong Ge
!!
!!
!! purpose    :  inversion of a matrix.
!!               For small matrices (n <= 8) the legacy full-pivot
!!               Gauss-Jordan implementation is used (LAPACK overhead
!!               dominates). For larger matrices the LAPACK
!!               symmetric-indefinite path is used: Bunch-Kaufman
!!               factorization (dsytrf) with diagonal pivoting followed by
!!               dsytri for the inverse. This handles both SPD and
!!               indefinite normal matrices that arise in kinematic LSQ.
!!               A general LU (dgetrf/dgetri) is kept as a safety net, and
!!               the legacy code is used as a last-resort fallback.
!!
!! parameters : a : matrix to be inverted
!!           ndim : dimension of a matrix
!!              n : dimension of a to be invert
!!              d : value of determinant a (0 on singularity, 1 on success)
!!

subroutine matinv(a, ndim, n, d)

  implicit none
  integer*4 ndim, n
  real*8 a(ndim, ndim), d
!
!! local
  integer*4 i, j, info, lwork
  integer*4, allocatable :: ipiv(:)
  real*8,    allocatable :: work(:)
  real*8,    allocatable :: asave(:, :)
  real*8 work_query(1)
  logical bk_ok

  if (n .lt. 0) then
    d = 0.d0
    return
  endif

!
!! Trivial / tiny matrices: use the legacy hand-rolled routine (it handles
!! the n==0 no-op case the same way the original implementation did, and
!! LAPACK call overhead would dominate for tiny n anyway). Callers like
!! polydf(idgr=0)-> matinv(n=0) rely on d=1 (success) for an empty matrix.
  if (n .le. 8) then
    call matinv_legacy(a, ndim, n, d)
    return
  endif

!
!! Save a copy in case Bunch-Kaufman fails and we need to fall back.
  allocate (asave(n, n))
  do j = 1, n
    do i = 1, n
      asave(i, j) = a(i, j)
    enddo
  enddo

!
!! Defensively symmetrize: mirror upper triangle into lower (the LSQ caller
!! only fills one triangle).
  do j = 1, n
    do i = j + 1, n
      a(i, j) = a(j, i)
    enddo
  enddo

!
!! Try Bunch-Kaufman symmetric-indefinite factorization.
  bk_ok = .false.
  allocate (ipiv(n))
  call dsytrf('U', n, a, ndim, ipiv, work_query, -1, info)
  if (info .eq. 0) then
    lwork = int(work_query(1))
    if (lwork .lt. n) lwork = n
    allocate (work(lwork))
    call dsytrf('U', n, a, ndim, ipiv, work, lwork, info)
    if (info .eq. 0) then
      call dsytri('U', n, a, ndim, ipiv, work, info)
      if (info .eq. 0) bk_ok = .true.
    endif
    deallocate (work)
  endif
  deallocate (ipiv)

  if (bk_ok) then
!   mirror upper triangle into lower so the full inverse is returned
    do j = 1, n
      do i = j + 1, n
        a(i, j) = a(j, i)
      enddo
    enddo
    d = 1.d0
    deallocate (asave)
    return
  endif

!
!! Bunch-Kaufman failed -- restore matrix and try general LU inverse.
  do j = 1, n
    do i = 1, n
      a(i, j) = asave(i, j)
    enddo
  enddo

  allocate (ipiv(n))
  call dgetrf(n, n, a, ndim, ipiv, info)
  if (info .eq. 0) then
    call dgetri(n, a, ndim, ipiv, work_query, -1, info)
    if (info .eq. 0) then
      lwork = int(work_query(1))
      if (lwork .lt. n) lwork = n
      allocate (work(lwork))
      call dgetri(n, a, ndim, ipiv, work, lwork, info)
      deallocate (work)
      if (info .eq. 0) then
        deallocate (ipiv)
        deallocate (asave)
        d = 1.d0
        return
      endif
    endif
  endif
  deallocate (ipiv)

!
!! Both LAPACK paths failed -- restore and fall through to legacy.
  do j = 1, n
    do i = 1, n
      a(i, j) = asave(i, j)
    enddo
  enddo
  deallocate (asave)

  call matinv_legacy(a, ndim, n, d)
  return
end subroutine

!!
!! Legacy full-pivot Gauss-Jordan implementation (kept for small matrices
!! and as a guaranteed-correct fallback). Same interface as matinv.
!!
subroutine matinv_legacy(a, ndim, n, d)

  implicit none
  integer*4 ndim, n
  real*8 a(ndim, ndim), d
!
!! local
  integer*4 i, j, k, l(ndim), m(ndim)
  real*8 max_a, swap

  d = 1.d0
!
!! loop over all row
  do k = 1, n
    l(k) = k
    m(k) = k
    max_a = a(k, k)
!
!! find out the diag. element with max. abs. value.
    do i = k, n
      do j = k, n
        if (dabs(max_a) - dabs(a(i, j)) .lt. 0.d0) then
          max_a = a(i, j)
          m(k) = i
          l(k) = j
        endif
      enddo
    enddo
!
!! rank defection or bad condition
    if (dabs(max_a) .lt. 1d-13) then
      d = 0.d0
      return
    endif
!
!! swap k and l(k) column
    if (l(k) .gt. k) then
      do i = 1, n
        swap = -a(i, k)
        a(i, k) = a(i, l(k))
        a(i, l(k)) = swap
      enddo
    endif
!
!! swap k and m(k) raw
    if (m(k) .gt. k) then
      do j = 1, n
        swap = -a(k, j)
        a(k, j) = a(m(k), j)
        a(m(k), j) = swap
      enddo
    endif
!
!! elemination
    do i = 1, n
      if (i .ne. k) a(k, i) = -a(k, i)/max_a
    enddo

    do i = 1, n
      do j = 1, n
        if (i .ne. k .and. j .ne. k) a(j, i) = a(j, i) + a(k, i)*a(j, k)
      enddo
    enddo
!
!! save invert part
    do j = 1, n
      if (j .ne. k) a(j, k) = a(j, k)/max_a
    enddo
    a(k, k) = 1/max_a
  enddo
!
!! re-range raw and column as input
  do k = n, 1, -1
    if (l(k) .gt. k) then
      do j = 1, n
        swap = a(k, j)
        a(k, j) = -a(l(k), j)
        a(l(k), j) = swap
      enddo
    endif

    if (m(k) .gt. k) then
      do i = 1, n
        swap = a(i, k)
        a(i, k) = -a(i, m(k))
        a(i, m(k)) = swap
      enddo
    endif
  enddo
  return
end subroutine
