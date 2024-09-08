subroutine otldisp(jd,sod,otlfil,disp)
  implicit none
  integer*4,intent(in) :: jd
  real*8,intent(in) :: sod
  character*(*),intent(in) :: otlfil
  real*8,intent(out) :: disp(3)
  !
  !! function called
  integer*4, external :: get_valid_unit
  real*8, external :: timdif
  !
  !! local
  integer*4 :: i,lfnotl,mjdl,mjdr,stat
  real*8 :: sodl,sodr,enul(3),enur(3) 
  lfnotl = get_valid_unit(10)
  !
  !! real otl file
  open(lfnotl,file=trim(otlfil),access='stream',form='unformatted',status="old")
  
  read(lfnotl,iostat=stat) mjdl,sodl,enul(1:3)
  do
    read(lfnotl,iostat=stat) mjdr,sodr,enur(1:3)
    if (stat .ne. 0) exit
    if (timdif(jd,sod,mjdr,sodr) .ge. 0.d0) then
      mjdl=mjdr; sodl=sodr; enul(1:3)=enur(1:3)
    elseif (timdif(jd,sod,mjdr,sodr) .lt. 0.d0) then
      exit
    endif
  end do
  close(lfnotl)

  if(abs(timdif(jd,sod,mjdl,sodl)) .ge. abs(timdif(jd,sod,mjdr,sodr))) then
    disp(1:3)=enur(1:3)
  else
    disp(1:3)=enul(1:3)
  endif

end subroutine