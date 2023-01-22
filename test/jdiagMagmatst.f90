program jdiagMagmatst
implicit none
real(8),allocatable :: a(:,:),w(:)
integer :: iarg,n,lda,i,j
character(10) :: arg1

iarg=command_argument_count()

if(iarg/=1) then
  call usage()
end if
call getarg(1,arg1)
read(arg1,*) n
lda=n

allocate(a(n,n),w(n))

do i=1,n
  do j=1,i
    a(i,j)=i+1
    a(j,i)=a(i,j)
  end do
end do

write(6,*) 'a'
do i=1,n
  write(6,'(20f10.4)') a(i,:)
end do

call magmaf_init() 
call jdiagMAGMA(lda,n,A,W)

write(6,*) 'w'
write(6,'(20f10.4)') w(:)

stop

contains

subroutine usage()
  write(6,'(a)') 'Usage : jdiagMagmatst n'
  write(6,'(a)') ' n    : dimension of a symmetric matrix for diagonalization'
  stop
end subroutine

subroutine jdiagMAGMA(lda,n,A,W) 
! A: input matrix, and eigenvectors will be stored here
! W: vector for eigenvalues
! n: matrix size   (4*NCC, 4*az atomok szama)
! lda: leading dimension
  use magma

  implicit none 

  integer :: info, lda, n 
  real(8) :: A(lda,n), W(n) 
  real(8), allocatable :: WORK(:) 
  integer, allocatable :: IWORK(:) 
  integer(4) :: lwork = -1  
  integer :: ngpu = 1 
  integer :: liwork = -1 
  integer ::status 

allocate (work(10000), IWORK(10000) )    !, stat=status) 
call magmaf_dsyevd_m (ngpu, 'V', 'L', n, A, lda, W, WORK, lwork, IWORK, liwork, info) 

lwork = int(work(1)) 
liwork = int(iwork(1)) 
deallocate (work, iwork)     

! allocate with the new, proper size:
allocate (work(lwork), IWORK(liwork))        !, stat=status) 
call magmaf_dsyevd_m (ngpu, 'V', 'L', n, A, lda, W, WORK, lwork, IWORK, liwork, info) 
call magmaf_finalize() 
return 
end

end program