program reshape_tst
  implicit none
  real,target :: x(2,2)
  real y(4)
  integer s1(1)
  integer s2(2)
  integer i,j
  real,pointer :: xp(:,:),yp(:)
  xp=>x
  write(*,*) shape(x)
  call sub0(x)
  do i=1,2
    write(*,*) (x(i,j),j=1,2)
  end do
  s1(1)=4
  y=reshape(x,s1)
  write(*,*) shape(y)
  call suba(y)
  write(*,*) x(:,:)
  write(*,*) y(:)
  s2(1)=2
  s2(2)=2
  x=reshape(y,s2)
  write(*,*) x(:,:)

contains
  subroutine sub0(x)
    real x(2,2)
    do i=1,2
      do j=1,2
        x(i,j)=j
      end do
    end do
  end subroutine

  subroutine suba(x)
    real x(4)
    write(*,*) x(:)
    do i=1,4
      x(i)=5-i
    end do
  end subroutine
end program
