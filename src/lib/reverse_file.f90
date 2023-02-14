!purpose: reverse file
!author: Wenyi Li (2023)
!fln: file name
!iunit: number for reading file
subroutine reverse_file(fln,iunit)
    character*(*) fln
    integer*4 iunit
    integer*4 ierr
    integer*4 line_num,j
    character(len=150) line
    character(len=150),dimension(:),allocatable :: lines
    line_num=0
    open(iunit,file=fln,status='old',iostat=ierr,err=100)
    do while(.true.)
       read(iunit,'(a)',end=300)  line
       line_num=line_num+1
    end do
300 rewind(iunit)
    allocate(lines(1:line_num))
    do j=1,line_num
       read(iunit,'(a)',end=200)  line
       lines(j)=trim(line)
    end do
200 close(iunit)
    open(iunit,file=fln,status='replace',iostat=ierr,err=100)
    do j=line_num,1,-1
        write(iunit,'(a)') trim(lines(j))
    end do
    deallocate(lines)
    close(iunit)
    return
100 continue
    write(*,*) '***ERROR(reverse_file):read file ',adjustl(trim(fln)),' failed'
    call exit(1)
end subroutine
