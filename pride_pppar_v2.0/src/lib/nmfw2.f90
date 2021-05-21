!!   INPUT  : lat,elev
!!
!!   OUTPUT : wetmf

subroutine nmfw2(lat, elev, wetmf)
implicit none

real*8 lat, elev,wetmf(2)
integer*4 i
real*8 a, b, c, beta, cose, gamma, sine, contop
real*8 wetmf_lats(5), w2p0_abc(5, 3)
real*8 dl, da, db, dc, la, deg2rad

data wetmf_lats/15.d0, 30.d0, 45.d0, 60.d0, 75.d0/

data w2p0_abc/ &
  5.8021897d-4, 5.6794847d-4, 5.8118019d-4, 5.9727542d-4, 6.1641693d-4, &
  1.4275268d-3, 1.5138625d-3, 1.4572752d-3, 1.5007428d-3, 1.7599082d-3, &
  4.3472961d-2, 4.6729510d-2, 4.3908931d-2, 4.4626982d-2, 5.4736038d-2/

a=0.d0
b=0.d0
c=0.d0
deg2rad = 3.14159265d0/180.d0
la = abs(lat)

if (la .le. wetmf_lats(1)) then
  a = w2p0_abc(1, 1)
  b = w2p0_abc(1, 2)
  c = w2p0_abc(1, 3)
endif

do i = 1, 4
  if (la .gt. wetmf_lats(i) .and. la .le. wetmf_lats(i + 1)) then
    dl = (la - wetmf_lats(i))/(wetmf_lats(i + 1) - wetmf_lats(i))
    da = w2p0_abc(i + 1, 1) - w2p0_abc(i, 1)
    a = w2p0_abc(i, 1) + dl*da
    db = w2p0_abc(i + 1, 2) - w2p0_abc(i, 2)
    b = w2p0_abc(i, 2) + dl*db
    dc = w2p0_abc(i + 1, 3) - w2p0_abc(i, 3)
    c = w2p0_abc(i, 3) + dl*dc
  endif
end do

if (la .ge. wetmf_lats(5)) then
  a = w2p0_abc(5, 1)
  b = w2p0_abc(5, 2)
  c = w2p0_abc(5, 3)
endif

sine = sin(elev*deg2rad)
cose = cos(elev*deg2rad)
beta = b/(sine + c)
gamma = a/(sine + beta)
contop = (1.d0 + a/(1.d0 + b/(1.d0 + c)))

wetmf(1) = contop/(sine + gamma)
wetmf(2) = -contop/(sine + gamma)**2*(cose - a/(sine + beta)**2*cose*(1.d0 - b/(sine + c)**2))

return
end

