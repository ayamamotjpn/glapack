program jdiagMagmatst
implicit none
real(8),allocatable :: ad(:,:),wd(:)
real,allocatable :: as(:,:),ws(:)
complex(8),allocatable :: az(:,:),wz(:)
complex,allocatable :: ac(:,:),wc(:)

integer :: iarg,n,ngpu,lda
integer :: i,j
character(10) :: arg1,arg2

iarg=command_argument_count()

if(iarg/=2) then
  call usage()
end if
call getarg(1,arg1)
call getarg(2,arg2)
read(arg1,*) n
read(arg2,*) ngpu
lda=n

allocate(as(n,n),ws(n))
allocate(ad(n,n),wd(n))
allocate(ac(n,n),wc(n))
allocate(az(n,n),wz(n))

do i=1,n
  do j=1,i
    as(i,j)=i+1
    as(j,i)=as(i,j)
  end do
end do

ad(:,:)=as(:,:)
ac(:,:)=as(:,:)
az(:,:)=ac(:,:)

write(6,*) 'as'
do i=1,n
  write(6,'(20f10.4)') as(i,:)
end do

call magmaf_init() 

call jdiagMAGMAd(ngpu,lda,n,ad,wd)

write(6,*) 'wd'
write(6,'(20f10.4)') wd(:)

stop

contains

subroutine usage()
  write(6,'(a)') 'Usage : jdiagMagmatst n'
  write(6,'(a)') ' n    : dimension of a symmetric matrix for diagonalization'
  write(6,'(a)') ' ngpu : number of GPUs used'
  stop
end subroutine

subroutine jdiagMAGMAd(ngpu,lda,n,A,W) 
! A: input matrix, and eigenvectors will be stored here
! W: vector for eigenvalues
! n: matrix size   (4*NCC, 4*az atomok szama)
! lda: leading dimension
! ngpu : number of GPUs used
  use magma

  implicit none 
  integer :: ngpu
  integer :: info, lda, n
  real(8) :: A(lda,n), W(n) 
  real(8), allocatable :: WORK(:) 
  integer, allocatable :: IWORK(:) 
  integer(4) :: lwork = -1  
  integer :: liwork = -1 
  integer ::status 

allocate (work(1), IWORK(1) )    !, stat=status)
! for estimating necessary work arrays
if(ngpu==1) then
  call magmaf_dsyevd ('V', 'L', n, A, lda, W, WORK, lwork, IWORK, liwork, info) 
else
  call magmaf_dsyevd_m (ngpu, 'V', 'L', n, A, lda, W, WORK, lwork, IWORK, liwork, info) 
end if

lwork = int(work(1)) 
liwork = int(iwork(1)) 
deallocate (work, iwork)     

! allocate with the new, proper size:
allocate (work(lwork), IWORK(liwork))        !, stat=status) 
if(ngpu==1) then
  call magmaf_dsyevd ( 'V', 'L', n, A, lda, W, WORK, lwork, IWORK, liwork, info)
else
  call magmaf_dsyevd_m (ngpu, 'V', 'L', n, A, lda, W, WORK, lwork, IWORK, liwork, info) 
end if

call magmaf_finalize() 
return 
end

end program
