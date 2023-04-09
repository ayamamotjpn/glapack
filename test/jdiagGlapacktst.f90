program jdiagGlapacktst
  use glapack
#ifdef MAGMA
  use magma
#endif
  implicit none
  real(8),allocatable :: az(:,:),wz(:)
  real,allocatable :: ac(:,:),wc(:)
  real(8),allocatable :: ad(:,:),wd(:)
  real,allocatable :: as(:,:),ws(:)

  integer :: iarg,ngpu,n,lda,i,j
  character(10) :: arg1,arg2

  iarg=command_argument_count()

  if(iarg/=2) then
    call usage()
  end if
  call getarg(1,arg1)
  call getarg(1,arg2)
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

  if(ngpu==1) then
    call glapack_init()
  else
    call glapack_init(ngpu)
  end if
  call jdiagGlapackd(lda,n,ad,wd)

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

  subroutine jdiagGlapackd(lda,n,A,W)
    use glapack
    ! A: input matrix, and eigenvectors will be stored here
    ! W: vector for eigenvalues
    ! n: matrix size   (4*NCC, 4*az atomok szama)
    ! lda: leading dimension
#ifdef MAGMA
    use magma
#endif
    implicit none

    integer :: info, lda, n
    real(8) :: A(lda,n), W(n)
    real(8), allocatable :: WORK(:)
    integer, allocatable :: IWORK(:)
    integer(4) :: lwork = -1
    integer :: ngpu = 1
    integer :: liwork = -1
    integer :: lrwork
    integer :: status

    call glapack_ev_prm('V', 'L', n, A, lda, lwork, lrwork, liwork) ! error occures here
!#ifdef MAGMA
!    !call magmaf_dsyevd_m (ngpu, 'V', 'L', n, A, lda, W, WORK, lwork, IWORK, liwork, info)
!#else
!    !call dsyev('V', 'L', n, A, lda, lwork, lrwork, liwork)
!#endif
    write(6,'(a,3i10)') 'lwork,lrwork,liwork=',lwork,lrwork,liwork ! for test
    ! allocate with the new, proper size:
    allocate (work(lwork), IWORK(liwork))        !, stat=status)
    call glapack_ev('V', 'L', n, A, lda, W, WORK, lwork, IWORK, liwork, info)
!#ifdef MAGMA
!    !call magmaf_dsyevd_m (ngpu, 'V', 'L', n, A, lda, W, WORK, lwork, IWORK, liwork, info)
!#else
!    !call dsyev('V', 'L', n, A, lda, lwork, lrwork, liwork)
!#endif

!#ifdef MAGMA 
!    call magmaf_finalize()  !only for gpu version
!call magmaf_init() #endif 
    return
  end

end program
