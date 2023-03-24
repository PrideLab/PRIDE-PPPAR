
! Initialize the variables --------------------------------------------------
! init ephemeris ------------------------------------------------------------
subroutine init_eph(eph)
implicit none
include 'file_para.h'
type(eph_t), intent(out) :: eph
eph=eph_t(0,0,0,0,0,0,0,0,&
        gtime_t(0,0.d0),gtime_t(0,0.d0),gtime_t(0,0.d0),&
        0d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
        0d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
        0d0,0.d0,0.d0)
end subroutine

subroutine init_geph(geph)
implicit none
include 'file_para.h'
type(geph_t), intent(out) :: geph
geph=geph_t(0,0,0,0,0,0,&
        gtime_t(0,0.d0),gtime_t(0,0.d0),&
        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
end subroutine

subroutine init_sol(sol)
implicit none
include 'file_para.h'
type(sol_t), intent(out) :: sol
sol=sol_t(gtime_t(0,0.d0),gtime_t(0,0.d0),0.d0,0.d0,0.d0,0,0,0,0.d0)
end subroutine

! initialize station parameter ----------------------------------------------
subroutine init_sta(sta)
implicit none
include 'file_para.h'
type(sta_t), intent(out) :: sta
sta=sta_t('','','','','','','',0,0,0,0.d0,0.d0,0.d0)
end subroutine

subroutine InitGlobal()
implicit none
include 'file_para.h'
type(sta_t) :: sta0
chisqr=(/&  ! chi-sqr(n) (alpha=0.001)
    10.8d0, 13.8d0, 16.3d0, 18.5d0, 20.5d0, 22.5d0, 24.3d0, 26.1d0, 27.9d0, 29.6d0,&
    31.3d0, 32.9d0, 34.5d0, 36.1d0, 37.7d0, 39.3d0, 40.8d0, 42.3d0, 43.8d0, 45.3d0,&
    46.8d0, 48.3d0, 49.7d0, 51.2d0, 52.6d0, 54.1d0, 55.5d0, 56.9d0, 58.3d0, 59.7d0,&
    61.1d0, 62.5d0, 63.9d0, 65.2d0, 66.6d0, 68.0d0, 69.3d0, 70.7d0, 72.1d0, 73.4d0,&
    74.7d0, 76.0d0, 77.3d0, 78.6d0, 80.0d0, 81.3d0, 82.6d0, 84.0d0, 85.4d0, 86.7d0,&
    88.0d0, 89.3d0, 90.6d0, 91.9d0, 93.3d0, 94.7d0, 96.0d0, 97.4d0, 98.7d0, 100.d0,&
    101.d0, 102.d0, 103.d0, 104.d0, 105.d0, 107.d0, 108.d0, 109.d0, 110.d0, 112.d0,&
    113.d0, 114.d0, 115.d0, 116.d0, 118.d0, 119.d0, 120.d0, 122.d0, 123.d0, 125.d0,&
    126.d0, 127.d0, 128.d0, 129.d0, 131.d0, 132.d0, 133.d0, 134.d0, 135.d0, 137.d0,&
    138.d0, 139.d0, 140.d0, 142.d0, 143.d0, 144.d0, 145.d0, 147.d0, 148.d0, 149.d0/)

ura_value=(/&
    2.4d0, 3.4d0, 4.85d0, 6.85d0, 9.65d0, 13.65d0, 24.d0, 48.d0, 96.d0, &
    192.d0, 384.d0, 768.d0, 1536.d0, 3072.d0, 6144.d0/)

timeoffset_=0.d0

! defaults processing options
prcopt_default = prcopt_t(&
    PMODE_SINGLE, 0, 2, or(SYS_GPS,or(SYS_GAL,or(SYS_GLO,SYS_CMP))),&  ! mode, soltype:forward, nf:L1+L2, navsys
    10.d0*D2R, 0, 3, 1, 1,&                                ! elmin, sateph:brdc, ionoopt:dual, tropopt:saas, niter
    30.d0, 30.d0, 0, '', (/0,0,0,0,1,0/), .false.)         ! maxinno, maxgdop, exsats, rnxopt, posopt, lsa

! defaults solution output options
solopt_default = solopt_t(&
    SOLF_XYZ,TIMES_GPST,0,3,&    ! posf,times,timef,timeu
    0,1,1,0,0,0,0,&              ! degf,outhead,outopt,outvel,datum,height,geoid
    0,0,0,0,&                    ! solstatic,sstat,trace,issingle
    0.d0,&                       ! nmeaintv
    " ","",&                     ! separator/program name
    0)                           ! maxsolstd

! global variables in rtklib.f90
solindex_=0
allocate(allsol_(MAXSOLNUM))
allsol_=sol_t(gtime_t(0,0.d0),gtime_t(0,0.d0),0.d0,0.d0,0.d0,0,0,0,0.d0)
headwritten_=0  ! output the header
clkfile_=''

obss=obs_t(0,0,0,null())
navs=nav_t(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
           null(),null(),null(),0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
           0.d0,0.d0,0.d0,0.d0,0,0.d0,0.d0,0.d0,0.d0,0.d0,'')
stas=sta_t('','','','','','','',0,0,0,0.d0,0.d0,0.d0)

nepoch=0; revs=0
iobsu=0; iobsr=0;
end subroutine
