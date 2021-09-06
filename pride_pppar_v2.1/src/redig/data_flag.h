!! GOOD   -0- not AMB BAD DEL NODATA.
!! OLDAMB -1- ambiguity flag from LOG file
!! NEWAMB -2- newly found ambiguity  
!! NEWBAD -3- temporerily used in check_jump
!! DELBAD -4- removed because of large residual. 
!! DELSHT -5- removed as short piece
!! DELRMS -6- removed by rmstest 
!! NODATA -9- no data
integer*4, parameter :: GOOD=0, OLDAMB=1, NEWAMB=2
integer*4, parameter :: NEWBAD=3, DELBAD=4
integer*4, parameter :: DELSHT=5, DELRMS=6
integer*4, parameter :: DELORB=7, DELCLK=8
integer*4, parameter :: NODATA=9
