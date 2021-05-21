!!   INPUT  : doy,lat,height,elev
!!
!!   OUTPUT : hmf

subroutine nmfh2p1(doy, lat, height, elev, hmf)
implicit none

real*8 doy,lat,height,elev,hmf(2)
!
integer*4 i
real*8 a, b, c, beta, cose, gamma, sine, contop
real*8 hs_km, la, dl, doy_atm, doyr_atm, cost, doy2rad, deg2rad
real*8 hmf_lats(5),avg_abc(5, 3), amp_abc(5, 3)
real*8 daavg, daamp, dbavg, dbamp, dcavg, dcamp
real*8 aavg, aamp, bavg, bamp, cavg, camp
real*8 a_ht, b_ht, c_ht, ht_corr_coef, ht_corr, dhcc_del, dht_corr_del

data hmf_lats/15.d0, 30.d0, 45.d0, 60.d0, 75.d0/

data avg_abc/ &
  1.2769934d-3, 1.2683230d-3, 1.2465397d-3, 1.2196049d-3, 1.2045996d-3, &
  2.9153695d-3, 2.9152299d-3, 2.9288445d-3, 2.9022565d-3, 2.9024912d-3, &
  62.610505d-3, 62.837393d-3, 63.721774d-3, 63.824265d-3, 64.258455d-3/

data amp_abc/ &
  0.0, 1.2709626d-5, 2.6523662d-5, 3.4000452d-5, 4.1202191d-5, &
  0.0, 2.1414979d-5, 3.0160779d-5, 7.2562722d-5, 11.723375d-5, &
  0.0, 9.0128400d-5, 4.3497037d-5, 84.795348d-5, 170.37206d-5/

a_ht = 2.53d-5
b_ht = 5.49d-3
c_ht = 1.14d-3
a = 0.d0
b = 0.d0
c = 0.d0
doy2rad = 2.d0*3.14159265d0/365.25d0
deg2rad = 3.14159265d0/180.d0
hs_km = height/1000.d0
doy_atm = doy - 28.d0

la = abs(lat)
if (lat .lt. 0.d0) doy_atm = doy_atm + 365.25d0/2.d0

doyr_atm = doy_atm*doy2rad

cost = cos(doyr_atm)

if (la .le. hmf_lats(1)) then
  a = avg_abc(1, 1)
  b = avg_abc(1, 2)
  c = avg_abc(1, 3)
endif


do i = 1, 4
  if (la .gt. hmf_lats(i) .and. la .le. hmf_lats(i + 1)) then
    dl = (la - hmf_lats(i))/(hmf_lats(i + 1) - hmf_lats(i))
    daavg = avg_abc(i + 1, 1) - avg_abc(i, 1)
    daamp = amp_abc(i + 1, 1) - amp_abc(i, 1)
    aavg = avg_abc(i, 1) + dl*daavg
    aamp = amp_abc(i, 1) + dl*daamp
    a = aavg - aamp*cost
    dbavg = avg_abc(i + 1, 2) - avg_abc(i, 2)
    dbamp = amp_abc(i + 1, 2) - amp_abc(i, 2)
    bavg = avg_abc(i, 2) + dl*dbavg
    bamp = amp_abc(i, 2) + dl*dbamp
    b = bavg - bamp*cost
    dcavg = avg_abc(i + 1, 3) - avg_abc(i, 3)
    dcamp = amp_abc(i + 1, 3) - amp_abc(i, 3)
    cavg = avg_abc(i, 3) + dl*dcavg
    camp = amp_abc(i, 3) + dl*dcamp
    c = cavg - camp*cost
  endif
end do


if (la .ge. hmf_lats(5)) then
  a = avg_abc(5, 1) - amp_abc(5, 1)*cost
  b = avg_abc(5, 2) - amp_abc(5, 2)*cost
  c = avg_abc(5, 3) - amp_abc(5, 3)*cost
endif


sine = sin(elev*deg2rad)
cose = cos(elev*deg2rad)
beta = b/(sine + c)
gamma = a/(sine + beta)
contop = (1.d0 + a/(1.d0 + b/(1.d0 + c)))

hmf(1) = contop/(sine + gamma)
hmf(2) = -contop*cose/(sine + gamma)**2*(1.d0 - a/(sine + beta)**2*(1.d0 - b/(sine + c)**2))

beta = b_ht/(sine + c_ht)
gamma = a_ht/(sine + beta)
contop = (1.d0 + a_ht/(1.d0 + b_ht/(1.d0 + c_ht)))
ht_corr_coef = 1/sine - contop/(sine + gamma)
ht_corr = ht_corr_coef*hs_km
hmf(1) = hmf(1) + ht_corr

dhcc_del = -cose/sine**2+ contop*cose/(sine + gamma)**2* &
           (1.d0 - a_ht/(sine + beta)**2*(1.d0 - b_ht/(sine + c_ht)**2))
dht_corr_del = dhcc_del*hs_km
hmf(2) = hmf(2) + dht_corr_del

return
end
