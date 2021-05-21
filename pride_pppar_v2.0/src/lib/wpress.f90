!!   INPUT  : Rhum (Relative humidity)
!!            Temp (Temperature)

real*8 function wpress(Rhum, Temp)
implicit none
!
real*8 Rhum, Temp
wpress = Rhum*6.11d0*10.d0**(7.5d0*Temp/(2.373d2+Temp))
return
end
