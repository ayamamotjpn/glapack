program jdiagGlapacktst
  use glapack
#ifdef MAGMA
  use magma
#endif
  implicit none
  real(8),allocatable :: a(:,:),w(:)
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

  if(ngpu>1) then
    call glapack_init()
  else
    call glapack_init(ngpu)
  end if
  call jdiagGlapack(lda,n,A,W)

  write(6,*) 'w'
  write(6,'(20f10.4)') w(:)

  stop

contains

  subroutine usage()
    write(6,'(a)') 'Usage : jdiagMagmatst n'
    write(6,'(a)') ' n    : dimension of a symmetric matrix for diagonalization'
    write(6,'(a)') ' ngpu : number of GPUs used'
    stop
  end subroutine

  subroutine jdiagGlapack(lda,n,A,W)
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
#ifdef MAGMA
    !call magmaf_dsyevd_m (ngpu, 'V', 'L', n, A, lda, W, WORK, lwork, IWORK, liwork, info)
#else
    !call dsyev('V', 'L', n, A, lda, lwork, lrwork, liwork)
#endif
    write(6,'(a,3i10)') 'lwork,lrwork,liwork=',lwork,lrwork,liwork ! for test
    ! allocate with the new, proper size:
    allocate (work(lwork), IWORK(liwork))        !, stat=status)
    call glapack_ev('V', 'L', n, A, lda, W, WORK, lwork, IWORK, liwork, info)
#ifdef MAGMA
    !call magmaf_dsyevd_m (ngpu, 'V', 'L', n, A, lda, W, WORK, lwork, IWORK, liwork, info)
#else
    !call dsyev('V', 'L', n, A, lda, lwork, lrwork, liwork)
#endif

#ifdef MAGMA 
    call magmaf_finalize()  !only for gpu version
#endif 
    return
  end

end program
