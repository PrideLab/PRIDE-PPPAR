
! File relevant function ----------------------------------------------------
! get fname from fpath
subroutine getfname(fpath, fname)
implicit none
include 'file_para.h'
character(*), intent(in) :: fpath
character(*), intent(out) :: fname
integer*4 ind1, ind2
fname=fpath
ind1=index(fpath,'/',BACK = .TRUE.)
ind2=index(fpath,'\\',BACK = .TRUE.)
if(ind1/=0)then
    fname=trim(fpath(ind1+1:)); return
endif
if(ind2/=0)then
    fname=trim(fpath(ind2+2:)); return
endif
end subroutine

! get fdir from fpath
subroutine getfdir(fpath, fdir)
implicit none
include 'file_para.h'
character(*), intent(in) :: fpath
character(*), intent(out) :: fdir
integer*4 ind1, ind2
fdir=""
ind1=index(fpath,'/',BACK = .TRUE.)
ind2=index(fpath,'\\',BACK = .TRUE.)
if(ind1/=0)then
    fdir=trim(fpath(1:ind1)); return
endif
if(ind2/=0)then
    fdir=trim(fpath(1:ind2+1)); return
endif
end subroutine

! confirm file status
! status ==> 0-none, 1-exist
integer*4 function flexist(flname)
implicit none
include 'file_para.h'
character(*), intent(in) :: flname
logical*2 :: alive
inquire(file=flname,exist=alive)
if(alive)then
    flexist=1; return
else
    flexist=0; return
endif
end function

! open output file for append -----------------------------------------------
integer*4 function openfile(outfile)
implicit none
include 'file_para.h'
character(*), intent(in) :: outfile
integer*4 :: fp=FPOUT, info
if(len_trim(outfile)/=0)then
    open(unit=fp,file=outfile,status='old',iostat=info,position='append')
    if(info/=0)then
        openfile=6; return
    else
        openfile=fp; return
    endif
else
    openfile=6; return
endif
end function
