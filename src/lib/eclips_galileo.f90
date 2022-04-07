!
!! eclips_galileo.f90
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
!! author: Chen Wang, Qiyuan Zhang, Songfeng Yang
!! 
!!
!!
Subroutine eclips_galileo(TTAG,BLK,SVBCOS,IECLIPS,XSV,SANTXYZ,VSVC,BETA)
implicit none
integer*4 ieclips
integer*4 j,i

real*8 ttag,svbcos
real*8 eclstm,ecletm
real*8 anight,mu
real*8 xsv(*),santxyz(*),vsvc(*)
real*8 beta,yangle,phi
character(len=*) blk

real*8 twohr,halfhr
real*8 dttag,anoon
real*8 murate,det,betadg
real*8 santx,santy,v(3),r(3)
real*8 yawend,sut(3)

logical*4 noon,night
real*8,parameter :: pi=3.14159265358979311599796346854419d0
real*8,parameter :: rad2deg=180.d0/pi
real*8,parameter :: deg2rad=pi/180.d0

ieclips=0
twohr=7200.d0
halfhr=1800.d0


MURATE=DSQRT((VSVC(1)**2+VSVC(2)**2+VSVC(3)**2)/(XSV(1)**2+XSV(2)**2+XSV(3)**2))*RAD2DEG
noon=.false.
night=.false.
betadg=beta*rad2deg
mu=svbcos/dcos(betadg*deg2rad)

select case(trim(blk))
  case('GALILEO-1')
    anoon=15.d0
    if(mu.lt.dcos(195.d0*deg2rad).and.dabs(betadg).le.2.d0)then
      night=.true.
    endif
    if(mu.gt.dcos(15.d0*deg2rad).and.dabs(betadg).le.2.d0)then
      noon=.true.
    endif
  case('GALILEO-2')
    anoon=10.d0
    if(mu.lt.dcos(190.d0*deg2rad).and.dabs(betadg).le.4.1d0)then
      night=.true.
    endif
    if(mu.gt.dcos(10.d0*deg2rad).and.dabs(betadg).le.4.1d0)then
      noon=.true.
    endif
  case default
    goto 1
end select

YANGLE=DACOS((SANTXYZ(1)*VSVC(1)+SANTXYZ(2)*VSVC(2)+SANTXYZ(3)*VSVC(3))/ &
                DSQRT(VSVC(1)**2+VSVC(2)**2+VSVC(3)**2))*RAD2DEG
                
IF (BETADG .GT. 0.d0) YANGLE=-YANGLE

if(night.or.noon)then

  DET=DSQRT((180.d0-DACOS(SVBCOS)*RAD2DEG)**2-BETADG**2)
  PHI=PI/2.d0
  
  if(night)then
    IF (DABS(YANGLE) .LT. 90.d0) DET=-DET
  endif
  
  IF (NOON) THEN
    DET=DSQRT((DACOS(SVBCOS)*RAD2DEG)**2-BETADG**2)
    IF (DABS(YANGLE) .GT. 90.d0) DET=-DET
  END IF
  

  eclstm=ttag+det/murate
  eclstm=eclstm-anoon/murate
  ecletm=eclstm+2.d0*anoon/murate

  if(ttag.ge.eclstm.and.ttag.le.ecletm) then
    do j=1,3
      v(j)=vsvc(j)/dsqrt(vsvc(1)**2+vsvc(2)**2+vsvc(3)**2)
      R(j)=xsv(j)/dsqrt(xsv(1)**2+xsv(2)**2+xsv(3)**2)
    enddo
    
    det=murate*(ecletm-eclstm)/2.d0
    
    if(svbcos.lt.0) then
      if(trim(blk).eq.'GALILEO-1') then
        sut(1)=-dsin(pi-det*deg2rad)*dcos(betadg*deg2rad)
        sut(2)=-dsin(betadg*deg2rad)
        sut(3)=-dcos(pi-det*deg2rad)*dcos(betadg*deg2rad)
        
        yawend=datan2(sut(2)/dsqrt(1.d0-sut(3)**2),sut(1)/dsqrt(1.d0-sut(3)**2))*rad2deg
        
        det=(ttag-eclstm)*murate-anoon+180.d0
        sut(1)=-dsin(det*deg2rad)*dcos(betadg*deg2rad)
        sut(2)=-dsin(betadg*deg2rad)
        sut(3)=-dcos(det*deg2rad)*dcos(betadg*deg2rad)
        
        sut(2)=0.5d0*(dsin(2.d0*deg2rad)*sign(1.d0,yawend)+sut(2))+&
                   0.5d0*(dsin(2.d0*deg2rad)*sign(1.d0,yawend)-sut(2))*dcos(pi*dabs(sut(1))/dsin(15.d0*deg2rad))
        sut(3)=dsqrt(1-sut(1)**2-sut(2)**2)*sign(1.d0,sut(3))
        
        phi=datan2(sut(2)/dsqrt(1.d0-sut(3)**2),sut(1)/dsqrt(1.d0-sut(3)**2))*rad2deg
      endif
      
      if(trim(blk).eq.'GALILEO-2') then
        yawend=datan2(-dtan(betadg*deg2rad),dsin(-det*deg2rad))*rad2deg
        phi=90.d0*sign(1.d0,yawend)+(yawend-90.d0*sign(1.d0,yawend))*dcos(2*pi/5656.d0*(ttag-eclstm))
        yawend=datan2(-dtan(betadg*deg2rad),dsin(det*deg2rad))*rad2deg
        if((yawend/phi).ge.1.d0.or.(phi/yawend).lt.0.d0) then
          phi=yawend
        endif
      endif
      ieclips=1
    else
      if(trim(blk).eq.'GALILEO-1') then
        sut(1)=-dsin(-det*deg2rad)*dcos(betadg*deg2rad)
        sut(2)=-dsin(betadg*deg2rad)
        sut(3)=-dcos(-det*deg2rad)*dcos(betadg*deg2rad)
        
        yawend=datan2(sut(2)/dsqrt(1.d0-sut(3)**2),sut(1)/dsqrt(1.d0-sut(3)**2))*rad2deg
        
        det=(ttag-eclstm)*murate-anoon+180.d0
        sut(1)=-dsin(det*deg2rad)*dcos(betadg*deg2rad)
        sut(2)=-dsin(betadg*deg2rad)
        sut(3)=-dcos(det*deg2rad)*dcos(betadg*deg2rad)
        
        sut(2)=0.5d0*(dsin(2.d0*deg2rad)*sign(1.d0,yawend)+sut(2))+&
                   0.5d0*(dsin(2.d0*deg2rad)*sign(1.d0,yawend)-sut(2))*dcos(pi*dabs(sut(1))/dsin(15.d0*deg2rad))
        sut(3)=dsqrt(1-sut(1)**2-sut(2)**2)*sign(1.d0,sut(3))
        
        phi=datan2(sut(2)/dsqrt(1.d0-sut(3)**2),sut(1)/dsqrt(1.d0-sut(3)**2))*rad2deg
      endif
      
      if(trim(blk).eq.'GALILEO-2') then
        yawend=datan2(-dtan(betadg*deg2rad),dsin(pi-det*deg2rad))*rad2deg
        phi=90.d0*sign(1.d0,yawend)+(yawend-90.d0*sign(1.d0,yawend))*dcos(2*pi/5656.d0*(ttag-eclstm))
        yawend=datan2(-dtan(betadg*deg2rad),dsin(pi+det*deg2rad))*rad2deg
        if((yawend/phi).ge.1.d0.or.(phi/yawend).lt.0.d0) then
          phi=yawend
        endif
      endif
      ieclips=2
    endif
    SANTX=(DCOS((PHI-YANGLE)*DEG2RAD)*(V(2)-V(3)*R(2)/R(3))-DCOS(PHI*DEG2RAD)* &
            (SANTXYZ(2)-SANTXYZ(3)*R(2)/R(3)))/(SANTXYZ(1)*V(2)-SANTXYZ(2)*v(1)+ &
            ((SANTXYZ(2)*V(3)-SANTXYZ(3)*V(2))*R(1)+(SANTXYZ(3)*V(1)- &
            SANTXYZ(1)*V(3))*R(2))/R(3))
    SANTY=(DCOS(PHI*DEG2RAD) - (V(1)-V(3)*R(1)/R(3))*SANTX)/(V(2)-V(3)*R(2)/R(3))
    SANTXYZ(1)=SANTX
    SANTXYZ(2)=SANTY
    SANTXYZ(3)=(-R(1)*SANTX-R(2)*SANTY)/R(3)   
  endif
endif

1 return

end subroutine