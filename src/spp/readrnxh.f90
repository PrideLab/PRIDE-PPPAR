
! Read rinex header ---------------------------------------------------------
subroutine readrnxh(fp, ver, mytype, sys, tsys, tobs, obs, nav, sta, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp
real*8, intent(out) :: ver
character(*), intent(out) :: mytype
integer*4, intent(out) :: sys, tsys, stat  ! 0-error, 1-normal
character(3), intent(out) :: tobs(NUMSYS,MAXOBSTYPE)  !tobs(:,:)
type(obs_t), intent(out) :: obs
type(nav_t), intent(out) :: nav
type(sta_t), intent(out) :: sta
real*8 bias
character(MAXRNXLEN) buff
character(20) label
integer*4 :: i, myblock, sat, info
i=0; myblock=0; sat=0; info=0
label=buff(61:)
ver=2.10; mytype=''; sys=SYS_GPS

do while(.true.)
    read(fp,"(A)",iostat=info) buff
    label=buff(61:)
    if(info/=0) exit
    if (len_trim(buff)<=60)then
        cycle
    elseif(index(label,'RINEX VERSION / TYPE')/=0)then
        read(buff(1:9),*,iostat=info) ver
        read(buff(21:40),"(A)",iostat=info) mytype
        ! satellite system 
        select case(buff(41:41))
        case(' ','G')
            sys=SYS_GPS; tsys=TSYS_GPS
        case('R')
            sys=SYS_GLO; tsys=TSYS_UTC
        case('E')
            sys=SYS_GAL; tsys=TSYS_GAL  ! v.2.12 
        case('S')
            sys=SYS_SBS; tsys=TSYS_GPS
        case('J')
            sys=SYS_QZS; tsys=TSYS_QZS  ! v.3.02 
        case('C')
            sys=SYS_CMP; tsys=TSYS_CMP  ! v.2.12 
        case('I')
            sys=SYS_IRN; tsys=TSYS_IRN  ! v.3.03 
        case('M')
            sys=SYS_NONE; tsys=TSYS_GPS  ! mixed 
        case default
        end select
        cycle
    elseif(index(label,'PGM / RUN BY / DATE')/=0)then
        cycle
    elseif(index(label,'COMMENT')/=0)then
        cycle
    endif
    ! file type 
    select case(mytype)
    case('O')
        call decode_obsh(fp,buff,ver,tsys,tobs,obs,nav,sta)
    case('N')
        call decode_navh(buff,nav)
    case('G')
        call decode_gnavh(buff,nav)
    end select
    if (index(label,'END OF HEADER')/=0)then
        stat=1; return
    endif
    i=i+1
    if (i>=MAXPOSHEAD .and. mytype==' ') exit  ! no rinex file 
enddo
stat=0
end subroutine
