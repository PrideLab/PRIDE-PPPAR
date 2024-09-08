subroutine get_spp_args(ts,te,ti,prcopt,rnxolist,nrnxo,rnxnlist,nrnxn,outfile)
  implicit none
  include 'file_para.h'
  type(gtime_t),intent(out) :: ts, te
  real*8,intent(out) :: ti
  type(prcopt_t),intent(out) :: prcopt
  character(1024),intent(out) :: rnxolist(FLNUMLIM)
  integer*4,intent(out) :: nrnxo
  character(1024),intent(out) :: rnxnlist(FLNUMLIM)
  integer*4,intent(out) :: nrnxn
  character(1024),intent(out) :: outfile
  
  ! local
  integer*4 :: argc, i, nfile, info, mjd_s, trnxo, trnxn
  character(1024), pointer :: argv(:)
  character(1024) :: buff(6), infile(MAXFILE), flntmp, flntmp2, navdir, obsdir
  real*8 :: es(6), ee(6)

  ! function
  integer*4,external :: ymd2mjd, getrnxtyp, flexist
  logical*1,external :: IsDayRight, IsTimeRight
  type(gtime_t),external :: epoch2time

  ! Initialize variables
  es=(/2000,1,1,0,0,0/); ee=(/2000,12,31,23,59,59/)
  infile=''; flntmp=""; flntmp2=""
  obsdir=""; navdir=""

  ! read cmd arguments
  argc=iargc()
  allocate(argv(argc))
  do i=1,argc
    call getarg(i,argv(i))
  enddo
  if(argc==0)then
    call printhelp(); call exit(0)
  endif
  do i=1,argc
    if(argv(i)=='-?'.or.argv(i)=='-h')then
      call printhelp(); call exit(0)
    endif
  enddo

  ! decode cmd arguments
  nfile=0; i=1
  do while(i<=argc)
    if(argv(i)=='-o'.and.(i+1)<=argc)then
      outfile=argv(i+1); i=i+1
    elseif(argv(i)=='-ts'.and.(i+2)<=argc)then
      call string_split(argv(i+1),'/',buff(1:3),3)
      call string_split(argv(i+2),':',buff(4:6),3)
      read(buff(1:3),*,iostat=info) es(1),es(2),es(3)
      mjd_s=ymd2mjd(int(es(1:3)))
      if(mjd_s<=51666) prcopt%lsa=.true.
      read(buff(4:6),*,iostat=info) es(4),es(5),es(6); i=i+2
      if(.not.IsTimeRight(int(es(4)),int(es(5)),es(6))) call exit(3)
      ts=epoch2time(es)
    elseif(argv(i)=='-te'.and.(i+2)<=argc)then
      call string_split(argv(i+1),'/',buff(1:3),3)
      call string_split(argv(i+2),':',buff(4:6),3)
      read(buff(1:3),*,iostat=info) ee(1),ee(2),ee(3)
      read(buff(4:6),*,iostat=info) ee(4),ee(5),ee(6); i=i+2
      if(.not.IsDayRight (int(ee(1)),int(ee(2)),int(ee(3)))) call exit(3)
      if(.not.IsTimeRight(int(ee(4)),int(ee(5)),ee(6))) call exit(3)
      te=epoch2time(ee)
    elseif(argv(i)=='-ti'.and.(i+1)<=argc)then
      read(argv(i+1),*,iostat=info) ti; i=i+1
    elseif(argv(i)=='-k'.and.(i+1)<=argc)then
      i=i+2; cycle
    elseif(argv(i)=='-trop'.and.(i+1)<=argc)then
      if(argv(i+1)=='non')  prcopt%tropopt=0
      if(argv(i+1)=='saas') prcopt%tropopt=1
      i=i+1
    elseif(argv(i)=='-c'.and.(i+1)<=argc)then
      clkfile_=argv(i+1); i=i+1
    elseif(argv(i)=='-elev'.and.(i+1)<=argc)then
      read(argv(i+1),*,iostat=info) prcopt%elmin
      prcopt%elmin=prcopt%elmin*D2R; i=i+1
    elseif(nfile<MAXFILE)then
      nfile=nfile+1; infile(nfile)=argv(i)
    endif
    i=i+1
  enddo
  deallocate(argv);
  
  ! lose rinex file
  if(nfile/=2)then  ! custom execution format
    if(nfile<=0) write(*,*) 'error : no input file'
    call exit(2)
  endif

  ! decode rinex file name
  do i=1,2
    call getfname(infile(i),flntmp)
    info=getrnxtyp(flntmp)
    if(mod(info,2)==1 .and. info>=1 .and. info<=8)then
      trnxn=info
      call getfdir(infile(i),navdir)
      rnxnlist(1)=flntmp  ! assuming the time is consistent
    endif
    if(mod(info,2)==0 .and. info>=0 .and. info<=8)then
      trnxo=info
      call getfdir(infile(i),obsdir)
      rnxolist(1)=flntmp
    endif
  enddo

  ! get all same format rinex obs files 
  if (rnxolist(1) .ne. "") then
    nrnxo=1
    flntmp=rnxolist(1)
    do i=1,FLNUMLIM-1
      call getrnx_nname(flntmp,flntmp2,i)
      if(flexist(trim(obsdir)//flntmp2)==1)then
        rnxolist(i+1)=flntmp2
        nrnxo=nrnxo+1
      else
        cycle
      endif
    enddo
    if (rnxnlist(1) .eq. "") then
      nrnxn=1
      call getfname(infile(1),flntmp)
      call getfname(infile(2),flntmp2)
      if (rnxolist(1)==flntmp) then
        rnxnlist(1)=flntmp2
        if (navdir .eq. "") call getfdir(infile(2),navdir)
      else if (rnxolist(1)==flntmp2) then
        rnxnlist(1)=flntmp
        if (navdir .eq. "") call getfdir(infile(1),navdir)
      endif
    endif
  endif

  ! get all same format rinex nav files 
  if (rnxnlist(1) .ne. "") then
    nrnxn=1
    flntmp=rnxnlist(1)
    do i=1,FLNUMLIM-1
      call getrnx_nname(flntmp,flntmp2,i)
      if(flexist(trim(navdir)//flntmp2)==1)then
        rnxnlist(i+1)=flntmp2
        nrnxn=nrnxn+1
      else
        cycle
      endif
    enddo
      if (rnxolist(1) .eq. "") then
        nrnxo=1
        call getfname(infile(1),flntmp)
        call getfname(infile(2),flntmp2)
        if (rnxnlist(1)==flntmp) then
          rnxolist(1)=flntmp2
          if (obsdir .eq. "") call getfdir(infile(2),obsdir)
        else if (rnxnlist(1)==flntmp2) then
          rnxolist(1)=flntmp
          if (obsdir .eq. "") call getfdir(infile(1),obsdir)
        endif
      endif
  endif
  if (nrnxo <= 0) nrnxo = 1
  if (nrnxn <= 0) nrnxn = 1
  do i=1,FLNUMLIM
    flntmp=trim(obsdir)//rnxolist(i)
    rnxolist(i)=flntmp
  enddo
  do i=1,FLNUMLIM
    flntmp=trim(navdir)//rnxnlist(i)
    rnxnlist(i)=flntmp
  enddo

end subroutine