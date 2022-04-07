
! String relevant function --------------------------------------------------
! string to number ----------------------------------------------------------
! convert substring in string to number
! args   : char   *s        I   string ('... nnn.nnn ...')
!          integer*4 i,n    I   substring position and width
! return : converted number (0.d0:error)
!----------------------------------------------------------------------------
real*8 function str2num(s,i,n)
implicit none
include 'file_para.h'
character(*), intent(in) :: s
integer*4, intent(in) :: i,n
integer*4 err
str2num=0.d0
if(i<1 .or. len_trim(s)<i+n-1 .or. n<1) return
read(s(i:i+n-1),*,iostat=err) str2num
if(err/=0) str2num=0.d0
end function

! string splitting ----------------------------------------------------------
subroutine string_split(source, splitter, buff2, size2)  ! string & number (expected to be extracted)
implicit none
include 'file_para.h'
character(*), intent(in) :: source, splitter
integer*4, intent(in) :: size2
character(*), intent(out) :: buff2(size2)
character(1024) str, split1, buff(1024)  ! stores the split string
integer*4 nsize  ! actual number of strings
buff=''; nsize=0

if(splitter=='')then
    str=trim(source)
    call simplify_str(str,'')
    if(len_trim(str)==0) return
    nsize=1
    buff(1)=trim(str)
    do while(index(trim(str),' ')/=0)
        nsize=nsize+1
        buff(nsize-1)=str(1:index(str,' ')-1)
        buff(nsize)  =trim(str(index(str,' ')+1:))
        str=trim(buff(nsize))
    enddo
else
    str=trim(source)
    call simplify_str(str,splitter)
    split1=trim(adjustl(splitter))
    
    if(len_trim(str)<=len_trim(split1))then
        if(str==split1) return
        nsize=1; buff(1)=trim(str); return
    endif
    
    do while(index(str(len_trim(str)-len_trim(split1)+1:),trim(split1))==1)  ! string clipping
        str=str(1:len_trim(str)-len_trim(split1))
    enddo
    do while(index(str,trim(split1))==1)  ! string clipping
        str=str(len_trim(split1)+1:)
        str=adjustl(str)
    enddo
    nsize=1; buff(1)=trim(str)
    do while(index(str,trim(split1))/=0)
        nsize=nsize+1
        buff(nsize-1)=str(1:index(str,trim(split1))-1)
        buff(nsize)  =trim(str(index(str,trim(split1))+len_trim(split1):))
        str=trim(buff(nsize))
    enddo
endif
if(size2<=1024)then
    buff2(1:size2)=buff(1:size2)
endif
end subroutine

!! string splitting ---------------------------------------------------------
!subroutine string_split(source, splitter, nsize, buff)
!implicit none
!include 'file_para.h'
!character(*), intent(in) :: source, splitter
!integer*4, intent(out) :: nsize
!character(*), intent(out) :: buff(:)
!character(1024) str, split1
!buff=''; nsize=0
!
!if(splitter=='')then
!    str=trim(source)
!    call simplify_str(str,'')
!    if(len_trim(str)==0) return
!    nsize=1
!    buff(1)=trim(str)
!    do while(index(trim(str),' ')/=0)
!        nsize=nsize+1
!        buff(nsize-1)=str(1:index(str,' ')-1)
!        buff(nsize)  =trim(str(index(str,' ')+1:))
!        str=trim(buff(nsize))
!    enddo
!else
!    str=trim(source)
!    call simplify_str(str,splitter)
!    split1=trim(adjustl(splitter))
!    
!    if(len_trim(str)<=len_trim(split1))then
!        if(str==split1) return
!        nsize=1; buff(1)=trim(str); return
!    endif
!    
!    do while(index(str(len_trim(str)-len_trim(split1)+1:),trim(split1))==1)  ! string clipping
!        str=str(1:len_trim(str)-len_trim(split1))
!    enddo
!    do while(index(str,trim(split1))==1)  ! string clipping
!        str=str(len_trim(split1)+1:)
!        str=adjustl(str)
!    enddo
!    nsize=1; buff(1)=trim(str)
!    do while(index(str,trim(split1))/=0)
!        nsize=nsize+1
!        buff(nsize-1)=str(1:index(str,trim(split1))-1)
!        buff(nsize)  =trim(str(index(str,trim(split1))+len_trim(split1):))
!        str=trim(buff(nsize))
!    enddo
!endif
!end subroutine

! system options buffer to options ------------------------------------------
subroutine simplify_str(str, splitter)
implicit none
include 'file_para.h'
character(*), intent(inout) :: str
character(*), intent(in) :: splitter
character(1024) split1, split2
integer*4 lensplit
split1=''; split2=''
lensplit=0
if(splitter=='')then
    str=trim(adjustl(str))
    do while(index(trim(str),'  ')/=0)
        str(index(trim(str),'  '):)=str(index(trim(str),'  ')+1:)
    enddo
else
    str=trim(adjustl(str))
    split1=trim(adjustl(splitter))
    split2=trim(split1)//trim(split1)
    lensplit=len_trim(split1)
    do while(index(str,trim(split2))/=0)
        str(index(str,trim(split2)):)=str(index(str,trim(split2))+lensplit:)
    enddo
endif
end subroutine

! set string without tail space ---------------------------------------------
subroutine setstr(dst, src, n)
implicit none
include 'file_para.h'
character(*), intent(out) :: dst
character(*), intent(in) :: src
integer*4, intent(in) :: n
if(n<1.or.n>len_trim(src))then
    dst=''; return
endif
dst=trim(src(1:n))
end subroutine
