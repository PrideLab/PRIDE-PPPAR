!
!! map_invnormal.f90
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
!! Contributor: Jianghui Geng, Jihang Lin
!!
!!
!!
!! purpose  : compute features for ambiguity resolution
!! parameter:
!!    input  : QN      -- inversion matrix info (type invm)
!!             q22     -- normal matrix (lower-triangular packed)
!!             bias_p  -- bias vector used in penalty calculation
!!    output : features   -- 8-element array storing
!
subroutine cal_features(ncad, q22, bias_p, chisq, ratio, nobs, features)
    implicit none
    integer, intent(in) :: ncad
    real*8, intent(inout) :: q22(*)
    real*8, intent(in) :: bias_p(*)
    real*8, intent(in) :: chisq, ratio, nobs
    real*8, intent(out) :: features(8)

    real*8 :: std_dev, trace, penalty, dismin, wratio, pratio, sum_p, ratio_o
    real*8 :: dist(100), add(2)
    integer :: i, j, k

    if (ncad .gt. 1) then
        trace = 0.d0
        std_dev = q22(1)
        do i = 1, ncad
            k = (i - 1) * ncad - (i - 2) * (i - 1) / 2 + 1
            if (std_dev .lt. q22(k)) then
                std_dev = q22(k)
            endif
            trace = trace + q22(k)
        enddo
        trace = dsqrt(trace / std_dev)
        do i = 1, ncad
            do j = 1, ncad - i + 1
                k = (i - 1) * ncad - (i - 2) * (i - 1) / 2 + j
                q22(k) = q22(k) / std_dev
            enddo
        enddo
        call ambpenalty(ncad, q22, bias_p, dist, add)
        sum_p = 0.d0
        do i = 2, 100
            sum_p = sum_p + exp(-0.5d0 * dist(i))
        end do
        penalty = sum_p / exp(-0.5d0 * dist(1))
        dismin = dist(2) - dist(1)
        ratio_o = dist(2)/dist(1)
        wratio = dismin / add(1)
        pratio = add(2) / add(1)
    else
        penalty = -1.d0
        dismin = -1.d0
        wratio = -1.d0
        ratio_o  = -1.d0
        pratio = -1.d0
        trace = -1.d0
    endif

 !   write (*, '(a,i4,4f10.3)') 'Narrow-lane AR(fin): ', ncad, chisq, ratio

    features(1) = 1.d0 * ncad
    features(2) = dismin
    features(3) = ratio_o
    features(4) = penalty
    features(5) = trace
    features(6) = wratio
    features(7) = pratio
    features(8) = nobs
    return
end subroutine cal_features
