!
!! codspp.f90
!!
!!    Copyright (C) 2021 by Wuhan University
!!
!!    This program is an open source software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License (version 3) as
!!    published by the Free Software Foundation.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License (version 3) for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!!
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : single point positioning using range observations to check
!!            and/or improving station position
!! parameter:
!!    input : lfncid,lfnrem -- tmp files
!!            jd,sod        -- checkint time
!!            OB       -- site information
!!            ite           -- # of iterations
!!            flag          -- flag for position initializtion
!!    output: deltax        -- position correction
!
subroutine codspp(lfncid, lfnrem, jd, sod, OB, ite, flag, deltax)
  implicit none
  include '../header/const.h'
  include '../header/rnxobs.h'

  integer*4 lfncid, lfnrem, jd, ite, flag
  real*8 sod, deltax(1:*)
  type(rnxobr) OB
!
!! local
  integer*4 nobs, i, j, k, iobs, iptx, ipt(MAXSAT)
  real*8 g, g1, det, sig0, sig1
  real*8 g_G,g1_G,g_R,g1_R,g_E,g1_E,g_C,g1_C,g_J,g1_J
  real*8 coef_G,coef_R,coef_E,coef_C,coef_J
  real*8 amat(MAXSAT, 9), omc(MAXSAT)
  real*8 c(9, 9), w(9), v(MAXSAT), solu(9), c_tmp(81)
  character*1 tmpsys_G,tmpsys_R,tmpsys_E,tmpsys_C,tmpsys_3,tmpsys_J
  integer*4 tmpsys_G_num,tmpsys_R_num,tmpsys_E_num,tmpsys_C_num,tmpsys_3_num,tmpsys_J_num,tmpsys_num
  integer*4 sysnum,min_sat_num
  character*1 sys(6)
  integer*4 prn_int
!
!! function called
  logical*1 chitst
  real*8 dot
 
  ! G
  g_G = freq1_G/freq2_G
  g1_G = g_G/(g_G*g_G - 1.d0)
  coef_G=g_G**2
  ! R
  g_R = 9.0/7.0
  g1_R = g_R/(g_R*g_R - 1.d0)
  coef_R=g_R**2
  ! E
  g_E = freq1_E/freq2_E
  g1_E = g_E/(g_E*g_E - 1.d0)
  coef_E=g_E**2
  ! C
  g_C = freq1_C/freq2_C
  g1_C = g_C/(g_C*g_C - 1.d0)
  coef_C=g_C**2
  ! J
  g_J = freq1_J/freq2_J
  g1_J = g_J/(g_J*g_J - 1.d0)
  coef_J=g_J**2

  tmpsys_G=''
  tmpsys_R=''
  tmpsys_E=''
  tmpsys_C=''
  tmpsys_3=''
  tmpsys_J=''
  tmpsys_G_num=0
  tmpsys_R_num=0
  tmpsys_E_num=0
  tmpsys_C_num=0
  tmpsys_3_num=0
  tmpsys_J_num=0
  do i = 1, OB%nprn
    if (OB%omc(i, 3) .ne. 0.d0) then
      if(OB%prn(i)(1:1) .eq. 'G') then
        tmpsys_G='G'
        tmpsys_G_num=tmpsys_G_num+1
      elseif(OB%prn(i)(1:1) .eq. 'R') then
        tmpsys_R='R'
        tmpsys_R_num=tmpsys_R_num+1
      elseif(OB%prn(i)(1:1) .eq. 'E') then
        tmpsys_E='E'
        tmpsys_E_num=tmpsys_E_num+1
      elseif(OB%prn(i)(1:1) .eq. 'C') then
        read(OB%prn(i),'(1x,i2)') prn_int
        if(prn_int .le. 17)then
          tmpsys_C='C'
          tmpsys_C_num=tmpsys_C_num+1
        else
          tmpsys_3='3'
          tmpsys_3_num=tmpsys_3_num+1
        endif
      elseif(OB%prn(i)(1:1) .eq. 'J') then
        tmpsys_J='J'
        tmpsys_J_num=tmpsys_J_num+1
      endif
    endif
  enddo

  sysnum=0
  sys=''
  if(tmpsys_G .eq. 'G')then
    sysnum=sysnum+1
    sys(sysnum)='G'
  endif
  if(tmpsys_R .eq. 'R')then
    sysnum=sysnum+1
    sys(sysnum)='R'
  endif
  if(tmpsys_E .eq. 'E')then
    sysnum=sysnum+1
    sys(sysnum)='E'
  endif
  if(tmpsys_C .eq. 'C')then
    sysnum=sysnum+1
    sys(sysnum)='C'
  endif
  if(tmpsys_3 .eq. '3')then
    sysnum=sysnum+1
    sys(sysnum)='3'
  endif
  if(tmpsys_J .eq. 'J')then
    sysnum=sysnum+1
    sys(sysnum)='J'
  endif
!
!! range omc and partial derivatives
  nobs = 0
  amat=0.d0
  do i = 1, OB%nprn
    if (OB%omc(i, 3) .ne. 0.d0) then
      nobs = nobs + 1
      do j = 1, 3
        amat(nobs, j) = OB%amat(j, i)
      enddo
!! receiver clock unit (m)
      do j=1,sysnum
        if(OB%prn(i)(1:1) .eq. 'G' .and. sys(j) .eq. 'G') then
          amat(nobs, 3+j) = 1.d0
          omc(nobs) = (OB%omc(i, 3)*coef_G-OB%omc(i, 4))/(coef_G-1.d0)
        elseif(OB%prn(i)(1:1) .eq. 'R' .and. sys(j) .eq. 'R') then
          amat(nobs, 3+j) = 1.d0
          omc(nobs) = (OB%omc(i, 3)*coef_R-OB%omc(i, 4))/(coef_R-1.d0)
        elseif(OB%prn(i)(1:1) .eq. 'E' .and. sys(j) .eq. 'E') then
          amat(nobs, 3+j) = 1.d0
          omc(nobs) = (OB%omc(i, 3)*coef_E-OB%omc(i, 4))/(coef_E-1.d0)
        elseif(OB%prn(i)(1:1) .eq. 'C') then
          read(OB%prn(i),'(1x,i2)') prn_int
          if(prn_int .le. 17)then
            if(sys(j) .eq. 'C')then
              amat(nobs, 3+j) = 1.d0
              omc(nobs) = (OB%omc(i, 3)*coef_C-OB%omc(i, 4))/(coef_C-1.d0)
            endif
          else
            if(sys(j) .eq. '3')then
              amat(nobs, 3+j) = 1.d0
              omc(nobs) = (OB%omc(i, 3)*coef_C-OB%omc(i, 4))/(coef_C-1.d0)
            endif
          endif
        elseif(OB%prn(i)(1:1) .eq. 'J' .and. sys(j) .eq. 'J') then
          amat(nobs, 3+j) = 1.d0
          omc(nobs) = (OB%omc(i, 3)*coef_J-OB%omc(i, 4))/(coef_J-1.d0)
        endif
      enddo
      ipt(nobs) = i
    endif
  enddo
!
!! at least min_sat_num satellites
  min_sat_num=3+sysnum
  iptx = 0
  deltax(1:10) = 0.d0
5 iobs = nobs
  if (iptx .ne. 0) iobs = nobs - 1
  if (iobs .gt. min_sat_num) then
    c = 0.d0
    w = 0.d0
    solu = 0.d0
    sig0 = 0.d0
    do j = 1, nobs
      sig0 = sig0 + omc(j)*omc(j)
    enddo
!
!! add to normal equation
    do i = 1, min_sat_num
      do j = 1, nobs
        w(i) = w(i) + amat(j, i)*omc(j)
      enddo
      do j = 1, i
        do k = 1, nobs
          c(i, j) = c(i, j) + amat(k, i)*amat(k, j)
        enddo
        c(j, i) = c(i, j)
      enddo
    enddo
!
!! solve normal equation
    do i = 1, min_sat_num
      do j= 1, min_sat_num
        c_tmp((i-1)*min_sat_num+j)=c(j,i)
      enddo
    enddo
    call matinv(c_tmp, min_sat_num, min_sat_num, det)
    do i = 1, min_sat_num
      do j= 1, min_sat_num
        c(j,i)=c_tmp((i-1)*min_sat_num+j)
      enddo
    enddo
    if (det .eq. 0.d0) then
      write (*, '(a)') '***ERROR(codspp): matrix singularity '
      call exit(1)
    endif
!
!! solution
    do i = 1, min_sat_num
      do j = 1, min_sat_num
        solu(i) = solu(i) + c(i, j)*w(j)
      enddo
    enddo
!
!! sigma0
    do i = 1, min_sat_num
      sig0 = sig0 - solu(i)*w(i)
    enddo
    sig0 = dsqrt(sig0/(iobs - min_sat_num))
    if (iptx .ne. 0 .and. chitst(-1, iobs - 1, sig1, sig0, 0.99d0)) then
!! in this condition, save the epoch for it may be related to pseudorange quality
      iptx = 0
      deltax(1:10) = 0.d0
      return
    endif
    sig1 = sig0
    deltax(1:min_sat_num) = solu(1:min_sat_num)
    deltax(10) = 0.d0
    do i = 1, 3
      deltax(10) = deltax(10) + deltax(i)*deltax(i)
    enddo
    deltax(10) = dsqrt(deltax(10))
!
!! set OB
    if (iptx .ne. 0) then
      ite = 0
      deltax(10) = 10.d0
      OB%obs(ipt(iptx), 1:4) = 0.d0
      write (*, '(a,i5,f8.1,a,a3)') '###WARNING(codspp): bad range in SIT at ', jd, sod, ' for SAT', OB%prn(ipt(iptx))
      if (lfncid .ne. 0 .and. lfnrem .ne. 0) then
        write (lfncid) 'de'
        write (lfnrem) 1, jd, sod, 1, ipt(iptx)
      endif
      return
    endif
!
!! reliability
    if (ite .eq. 10 .or. sig0 .gt. 0.5 .and. deltax(10) .lt. 1.d0 .and. dabs(deltax(4)/VLIGHT) .lt. 1.d-6 &
       .and. dabs(deltax(5)/VLIGHT) .lt. 1.d-6 .and. dabs(deltax(6)/VLIGHT) .lt. 1.d-6 .and. dabs(deltax(7)/VLIGHT) .lt. 1.d-6 &
       .and. dabs(deltax(8)/VLIGHT) .lt. 1.d-6 .and. dabs(deltax(9)/VLIGHT) .lt. 1.d-6) then
!
!! residual vector
      do i = 1, nobs
        v(i) = -omc(i)
        do j = 1, min_sat_num
          v(i) = v(i) + amat(i, j)*deltax(j)
        enddo
      enddo
!
!! correlation analysis
      call cltasy(nobs, min_sat_num, amat, c, v, iptx)
      if (iptx .ne. 0) then
        omc(iptx) = 0.d0
        amat(iptx, 1:min_sat_num) = 0.d0
        if(OB%prn(ipt(iptx))(1:1) .eq. 'G')then
          tmpsys_G_num=tmpsys_G_num-1
          tmpsys_num=tmpsys_G_num
        elseif(OB%prn(ipt(iptx))(1:1) .eq. 'R')then
          tmpsys_R_num=tmpsys_R_num-1
          tmpsys_num=tmpsys_R_num
        elseif(OB%prn(ipt(iptx))(1:1) .eq. 'E')then
          tmpsys_E_num=tmpsys_E_num-1
          tmpsys_num=tmpsys_E_num
        elseif(OB%prn(ipt(iptx))(1:1) .eq. 'C')then
          tmpsys_C_num=tmpsys_C_num-1
          tmpsys_num=tmpsys_C_num
        elseif(OB%prn(ipt(iptx))(1:1) .eq. '3')then
          tmpsys_3_num=tmpsys_3_num-1
          tmpsys_num=tmpsys_3_num
        elseif(OB%prn(ipt(iptx))(1:1) .eq. 'J')then
          tmpsys_J_num=tmpsys_J_num-1
          tmpsys_num=tmpsys_J_num
        endif
        if(tmpsys_num .eq. 0)min_sat_num=min_sat_num-1
        goto 5
      else
        deltax(1:10) = 0.d0
      endif
    endif
  else
    deltax(1:min_sat_num) = 0.d0
    if (flag .eq. 2) then
      deltax(10) = 10.d0
      write (*, '(a,i5,f8.1)') '###WARNING(codspp): initialization fails in SIT at ', jd, sod
      do i = 1, nobs
        OB%obs(ipt(i), 1:4) = 0.d0
      enddo
      if (lfncid .ne. 0 .and. lfnrem .ne. 0) then
        write (lfncid) 'de'
        write (lfnrem) 1, jd, sod, nobs, (ipt(i), i=1, nobs)
      endif
    else
      deltax(10) = 0.d0
    endif
  endif

  return
end
