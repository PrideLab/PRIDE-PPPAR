!
!! added by ranzeng

subroutine ambpenalty(ncad, q22, bias, disall, add)
  implicit none

  integer*4 ncad,m,i,j,k,info
  real*8 q22(1:*), bias(1:*), disall(1:*),add(2)
  real*8, dimension(:),allocatable::a,s
  real*8, dimension(:),allocatable::Q,zn

!
!! local
  real*8 dump

  if (ncad .gt. 1) then
    allocate(a(ncad))
    allocate(Q(ncad*ncad))
  
    m=100
    allocate(s(m))
    allocate(zn(ncad*m))
    
    k=0
    do i=1,ncad
      a(i)=bias(i)
      do j=1,ncad
        if(j.ge.i)then
            k=k+1
            Q(j+(i-1)*ncad)=q22(k)
            Q(i+(j-1)*ncad)=q22(k)
        end if
      enddo
    enddo

    call lambda(ncad,m,a,Q,zn,s, add,info)

    do i=1,ncad
      bias(i)=zn(i)
    end do

    do i=1,m
      disall(i)=s(i)
    end do

    deallocate(a)
    deallocate(Q)
    deallocate(s)
    deallocate(zn)
  else
    dump = bias(1)
    bias(1) = nint(bias(1))*1.d0
    dump = bias(1) - dump
    disall(1) = dump/q22(1)*dump
    dump = 1.d0 - dabs(dump)
    disall(2) = dump/q22(1)*dump

  endif

end
