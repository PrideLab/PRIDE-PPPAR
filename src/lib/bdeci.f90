!!
!!  INPUT  : estval,sigm,uni,cutdev,cutsig
!!
!!  OUTPUT : prob,decis
!!    Calculate the decision-function value as per Appendix A
!!    of Dong and Bock [1989].
!!       Da-nan Dong  880403
!!    Add damping factor. DND  880602
!!    Change definition of decision function based on hypothesis
!!       test theory.      DND 880929
!!    Remove arbitrary factor of 3., scale wide-lane but not narrow-lane
!!       sigmas.    King 930322
!!
!!    Input:
!!      est:     estimated (real) bias value
!!      sigma:   estimated uncertainty of bias value, scaled by nrms for
!!                 widelane, unscaled for narrow lane
!!      ih :     control for receiver ambiguity
!!                 = 1  unit=one cycle
!!                 = 2  unit=half cycle
!!      cutdev:  threshold deviation for taper function (= 0.4 in Dong and Bock;
!!                 default still 0.4 here and in FIXDRV)
!!      cutsig:  threshold sigma for taper function (=0.33 in Dong and Bock;
!!                 default still 0.4 here and in FIXDRV)
!!
!!    Output:
!!      prob : area of decision-function region, set = 1.0 for
!!             each bias on input to BDECI, reduced by the probability
!!             of an error in rounding each bias. (Since we now search
!!             only one bias at a time, the cumulative probability is
!!             not calculated by NBIASR or NBIASP.)
!!      deci : decision function d(x,sigma) = FT/Q, inverse of
!!             1. - allowable rate for a type 1 error (alpha),
!!             compared with input wlcut or nlcut in NBIASR.
!!             (But F is no longer used because we search one bias at
!!             a time.)

subroutine bdeci(estval, sigm, uni, cutdev, cutsig, prob, decis)
implicit none
include '../header/const.h'

real*8 estval, sigm, cutdev, cutsig, prob, dev, decis
real*8 tm1, tm2, c, d1, a1, bint, taper, add, cdev, csig, trun, s1, s2
real*8 b1, b2, gammpq,  erfcb1, erfcb2
integer*4 uni, j
external gammpq

s2=1.414213562373095d0

prob = 1.d0

add = estval*dble(uni)
bint = dint(add + 0.5d0*dsign(1.d0, add))/dble(uni)
dev = dabs(bint - estval)

cdev = cutdev
if (cutdev .lt. 1.d-3) cdev = 0.15d0
if (uni .eq. 2) cdev = 0.5d0*cdev

if (dev .ge. cdev) then
  prob = 0.d0
  cdev = 0.5d0
endif

s1 = 1.d0/(sigm*s2)
trun = 1.d-9

tm1 = 1.d0 - dev/cdev
csig = cutsig
if (cutsig .lt. 1.d-3) csig = 0.15d0
if (sigm .gt. csig) then
  prob = 0.d0
  csig = 0.5d0
endif
tm2 = (csig - sigm)*3.d0
if (tm2 .lt. 0.d0) tm2 = 0.d0
taper = tm1**2*tm2

c = 0.d0
do j = 1, 50
  a1 = dble(j)
  b1 = (a1 - dev)*s1
  b2 = (a1 + dev)*s1
  if (b1 .lt. 0.d0 .or. b1 .gt. 15.d0) then
    erfcb1 = 0.d0
  else
    if (b1 .lt. 0.d0) then
      erfcb1 = 1.d0 + gammpq(0, 0.5d0, b1**2)
    else
      erfcb1 = gammpq(1, 0.5d0, b1**2)
    endif
  endif
  if (b2 .lt. 0.d0 .or. b2 .gt. 15.d0) then
    erfcb2 = 0.d0
  else
    if (b2 .lt. 0.d0) then
      erfcb2 = 1.d0 + gammpq(0, 0.5d0, b2**2)
    else
      erfcb2 = gammpq(1, 0.5d0, b2**2)
    endif
  endif
  d1 = erfcb1 - erfcb2
  c = c + d1
  if (d1 .lt. trun) goto 100
enddo

100 continue
if (prob .gt. 0.d0) prob = 1.d0 - c
if (c .lt. 1.d-9) c = 1.d-9
decis = taper/c
continue
return
end

!!!
real*8 function gammpq(typ, a, x)
implicit none
integer*4 typ ! 0:p 1:q
real*8  a, x, gln, mcf_gam,ser_gam
if (x .lt. 0.d0 .or. a .le. 0.d0) read(*,*) !pause
if(typ .eq. 0) then
  if (x .lt. a + 1.d0) then
    call ser_g(gammpq, a, x, gln)
  else
    call cf_g(mcf_gam, a, x, gln)
    gammpq = 1.d0 - mcf_gam
  endif
else if(typ .eq. 1) then
  if (x .lt. a + 1.d0) then
    call ser_g(ser_gam, a, x, gln)
    gammpq = 1.d0 - ser_gam
  else
    call cf_g(gammpq, a, x, gln)
  endif
endif
return
end

real*8 function mln_gam(xx)
implicit none
real*8 cof(6), stp, one, x, xx, tmp, ser
integer*4 j
data cof, stp/76.18009173d0, -86.50532033d0, 24.01409822d0, -1.231739516d0, .120858003d-2, -.536382d-5, 2.50662827465d0/
one=1.0d0
x = xx - one
tmp = x + 5.5d0
tmp = (x + 0.5d0)*log(tmp) - tmp
ser = one
do j = 1, 6
  x = x + one
  ser = ser + cof(j)/x
enddo
mln_gam = tmp + log(stp*ser)
return
end

subroutine cf_g(mcf_gam, a, x, gln)
implicit none
real*8 mcf_gam, mln_gam, a, an_f, x, gln, gold, a0, a1, b0, b1, fac, an, an_a, g, eps
integer*4 n, itmax
parameter(itmax=100, eps=3.d-7)
gln = mln_gam(a)
g = 0.d0
gold = 0.d0
a0 = 1.d0
a1 = x
b0 = 0.d0
b1 = 1.d0
fac = 1.d0
do n = 1, itmax
  an = float(n)
  an_a = an - a
  a0 = (a1 + a0*an_a)*fac
  b0 = (b1 + b0*an_a)*fac
  an_f = an*fac
  a1 = x*a0 + an_f*a1
  b1 = x*b0 + an_f*b1
  if (a1 .ne. 0.) then
    fac = 1.d0/a1
    g = b1*fac
    if (abs((g - gold)/g) .lt. eps) then
      mcf_gam = dexp(-x + a*dlog(x) - gln)*g
      return
    endif
    gold = g
  endif
enddo
write (*, '(a)') '***ERROR(bdeci/gcf): a too large, itmax to small '
call exit(1)
end

subroutine ser_g(ser_gam, a, x, gln)
implicit none
real*8 ser_gam, a, x, gln, mln_gam, ap, sum, del, eps
integer*4 n, itmax
parameter(itmax=100, eps=3.d-7)
gln = mln_gam(a)
if (x .le. 0.d0) then
  if (x .lt. 0.d0) then
    write (*, '(a)') '***ERROR(bdeci/gser): x < 0 '
    call exit(1)
  endif
  ser_gam = 0.d0
  return
endif
ap = a
sum = 1.d0/a
del = sum
do n = 1, itmax
  ap = ap + 1.d0
  del = del*x/ap
  sum = sum + del
  if (abs(del) .lt. abs(sum)*eps) then
    ser_gam = sum*exp(-x + a*log(x) - gln)
    return
  endif
enddo
write (*, '(a)') '***ERROR(bdeci/gser): a too large, itmax too small'
call exit(1)
end
