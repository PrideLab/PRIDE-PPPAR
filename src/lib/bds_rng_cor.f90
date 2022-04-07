!
!! bds_rng_cor.f90
!!
!!    Copyright (C) 2021 by Wuhan University
!!
!!    This program belongs to PRIDE PPP-AR which is an open source software:
!!    you can redistribute it and/or modify it under the terms of the GNU
!!    General Public License (version 3) as published by the Free Software Foundation.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License (version 3) for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!!
!! Contributor: Jianghui Geng, Songfeng Yang, Yuanxin Pan
!! 
!!
!!
!! prn     : C[0-9][0-9]
!! nchannel: 1-3 -> B1-B3
!! elev    : rad,[0-pi/2]
!! function bds_rng_cor : in meters
!
real*8 function bds_rng_cor(prn,nchannel,elev)
    character(len=3)::prn
    integer*4         nchannel
    real*8            elev
    ! local
    logical   lfirst/.true./
    real*8    cor_coefs(10,9),pi,elev_deg,cor(2)
    save      cor_coefs,pi
    integer*4 col,row,i,j,k
    real*8    coefs(90)
    data coefs/-0.55,  -0.71,  -0.27,   -0.55,  -0.71,  -0.27,   -0.47,  -0.40,  -0.22, &
               -0.40,  -0.36,  -0.23,   -0.40,  -0.36,  -0.23,   -0.38,  -0.31,  -0.15, &
               -0.34,  -0.33,  -0.21,   -0.34,  -0.33,  -0.21,   -0.32,  -0.26,  -0.13, &
               -0.23,  -0.19,  -0.15,   -0.23,  -0.19,  -0.15,   -0.23,  -0.18,  -0.10, &
               -0.15,  -0.14,  -0.11,   -0.15,  -0.14,  -0.11,   -0.11,  -0.06,  -0.04, &
               -0.04,  -0.03,  -0.04,   -0.04,  -0.03,  -0.04,    0.06,   0.09,   0.05, &
                0.09,   0.08,   0.05,    0.09,   0.08,   0.05,    0.34,   0.28,   0.14, &
                0.19,   0.17,   0.14,    0.19,   0.17,   0.14,    0.69,   0.48,   0.27, &
                0.27,   0.24,   0.19,    0.27,   0.24,   0.19,    0.97,   0.64,   0.36, &
                0.35,   0.33,   0.32,    0.35,   0.33,   0.32,    1.05,   0.69,   0.47/
    if(lfirst) then
      do j = 1, 10
        do k = 1, 9
          cor_coefs(j,k) = coefs((j-1)*9+k)
        enddo
      enddo
      lfirst=.false.
      pi=4*datan(1.d0)
    endif
    ! calc row & col
    read(prn(2:3),'(i2)') col
    if(col.lt.6) then
        col=0
    else if(col.lt.11.or.col.eq.13) then
        col=1
    else
        col=2
    endif
    col=col*3+nchannel
    elev_deg=elev/pi*180
    row=idint(elev_deg/10)+1
    ! interpolation
    cor(1)=cor_coefs(row,  col)
    cor(2)=cor_coefs(row+1,col)
    bds_rng_cor=cor(1)+dmod(elev_deg,10.d0)*(cor(2)-cor(1))/10.d0
    return

end
