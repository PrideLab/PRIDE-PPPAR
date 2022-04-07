
! Add rinex data ------------------------------------------------------------
! add ephemeris to navigation data ------------------------------------------
subroutine add_eph(nav, eph, stat)
implicit none
include 'file_para.h'
type(nav_t), intent(out) :: nav
type(eph_t), intent(in) :: eph
integer*4, intent(out) :: stat
type(eph_t), pointer :: nav_eph(:)

if(nav%nmax<=nav%n)then
    nav%nmax=nav%nmax+1024
    allocate(nav_eph(nav%nmax))
    nav_eph(1:nav%n)=nav%eph(1:nav%n)
    if(nav%n/=0) deallocate(nav%eph)
    nav%eph=>nav_eph
    nullify(nav_eph)
endif
nav%n=nav%n+1
nav%eph(nav%n)=eph
stat=1
end subroutine

subroutine add_geph(nav, geph, stat)
implicit none
include 'file_para.h'
type(nav_t), intent(out) :: nav
type(geph_t), intent(in) :: geph
integer*4, intent(out) :: stat
type(geph_t), pointer :: nav_geph(:)

if(nav%ngmax<=nav%ng)then
    nav%ngmax=nav%ngmax+1024
    allocate(nav_geph(nav%ngmax))
    nav_geph(1:nav%ng)=nav%geph(1:nav%ng)
    if(nav%ng/=0) deallocate(nav%geph)
    nav%geph=>nav_geph
    nullify(nav_geph)
endif
nav%ng=nav%ng+1
nav%geph(nav%ng)=geph
stat=1
end subroutine

! add obs data --------------------------------------------------------------
subroutine addobsdata(obs, mydata1, stat)
implicit none
include 'file_para.h'
type(obs_t), intent(out) :: obs
type(obsd_t), intent(in) :: mydata1
integer*4, intent(out) :: stat  ! 0-error, 1-normal

type(obsd_t), pointer :: obs_data(:)
if (obs%nmax<=obs%n)then
    if (obs%nmax<=0)then
        obs%nmax=NINCOBS
    else
        !obs%nmax=obs%nmax*2
        obs%nmax=obs%nmax+NINCOBS
    endif
    allocate(obs_data(obs%nmax))
    obs_data(1:obs%n)=obs%mydata(1:obs%n)
    if(obs%n/=0) deallocate(obs%mydata)
    obs%mydata=>obs_data
    nullify(obs_data)
endif

obs%n=obs%n+1
obs%mydata(obs%n)=mydata1
stat=1
end subroutine
