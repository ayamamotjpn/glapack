! 
!  Permission is hereby granted, free of charge, to any person obtaining
!  a copy of this software and associated documentation files (the "Software"),
!  to deal in the Software without restriction, including without limitation
!  the rights to use, copy, modify, merge, publish, distribute, sublicense,
!  and/or sell copies of the Software, and to permit persons to whom the
!  Software is furnished to do so, subject to the following conditions:
!  
!  The above copyright notice and this permission notice shall be included in
!  all copies or substantial portions of the Software.
!  
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
! 

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
