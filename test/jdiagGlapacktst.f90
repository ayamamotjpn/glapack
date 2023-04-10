program jdiagGlapacktst
  use glapack
  use kinds

#ifdef MAGMA
  use magma
#endif
  implicit none
  real(8),allocatable :: az(:,:),wz(:)
  real,allocatable :: ac(:,:),wc(:)
  real(8),allocatable :: ad(:,:),wd(:)
  real,allocatable :: as(:,:),ws(:)
  real(DP),allocatable :: a(:,:),w(:)

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
  allocate(a(n,n),w(n))   ! a(DP)=ad

  do i=1,n
    do j=1,i
      as(i,j)=i+1
      as(j,i)=as(i,j)
    end do
  end do

  write(6,*) 'as'
  do i=1,n
    write(6,'(20f10.4)') as(i,:)
  end do

  ad(:,:)=as(:,:)
  ac(:,:)=as(:,:)
  az(:,:)=ac(:,:)
  a(:,:)=as(:,:)

  if(ngpu==1) then
    call glapack_init()
  else
    call glapack_init(ngpu)
  end if

  print *,'before jdiagGlapackd'
  call jdiagGlapackd(lda,n,ad,wd)
  write(6,*) 'wd'
  write(6,'(a,20f10.4)') 'jdiagGlapackd', wd(:)
  ad(:,:)=as(:,:)

  print *,'before jdiagGlapack(w)'
  ! w should be equal to wd for DP=8
  call jdiagGlapack(lda,n,a,w)
  write(6,*) 'w'
  write(6,'(a,20f10.4)') 'jdiagGlapack(w)', w(:)

  print *,'before jdiagGlapack(ws)'
  call jdiagGlapack(lda,n,as,ws)
  write(6,*) 'ws'
  write(6,'(a,20f10.4)') 'jdiagGlapack(ws)', ws(:)

  print *,'before jdiagGlapack(wd)'
  call jdiagGlapack(lda,n,ad,wd)
  write(6,*) 'wd'
  write(6,'(a,20f10.4)') 'jdiagGlapack(ws)', wd(:)

  print *,'before jdiagGlapack(wc)'
  call jdiagGlapack(lda,n,ac,wc)
  write(6,*) 'wc'
  write(6,'(a,20f10.4)') 'jdiagGlapack(wc)', wc(:)

  print *,'before jdiagGlapack(wz)'
  call jdiagGlapack(lda,n,az,wz)
  write(6,*) 'wz'
  write(6,'(a,20f10.4)') 'jdiagGlapack(wz)', wz(:)

  stop

contains

  subroutine usage()
    write(6,'(a)') 'Usage : jdiagMagmatst n'
    write(6,'(a)') ' n    : dimension of a symmetric matrix for diagonalization'
    write(6,'(a)') ' ngpu : number of GPUs used'
    stop
  end subroutine

  ! genelic function
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
    class(*) :: A(lda,n),W(n)

    class(*), allocatable :: work(:)
    class(*), allocatable :: rwork(:)
    integer, allocatable :: iwork(:)
    integer(4) :: lwork = -1
    integer :: ngpu = 1
    integer :: liwork = -1
    integer :: lrwork
    integer :: status
    write(6,'(a,2i5)') 'n,lda=',n,lda  ! for test

    select type(A)
      type is(real)
      write(6,'(a)') 'A is real'
      call glapack_ev_prm('V', 'L', n, A, lda, lwork, lrwork, liwork) ! error occures here
      allocate(real::work(lwork))
      type is(real(8))
      write(6,'(a)') 'A is real(8)'
      call glapack_ev_prm('V', 'L', n, A, lda, lwork, lrwork, liwork) ! error occures here
      allocate(real(8)::work(lwork))
      type is(complex)
       write(6,'(a)') 'A is complex'
      call glapack_ev_prm('V', 'L', n, A, lda, lwork, lrwork, liwork) ! error occures here
      allocate(complex::work(lwork))
      allocate(real::rwork(lrwork))
      type is(complex(8))
      write(6,'(a)') 'A is complex(8)'
      call glapack_ev_prm('V', 'L', n, A, lda, lwork, lrwork, liwork) ! error occures here
      allocate(complex(8)::work(lwork))
      allocate(real(8)::rwork(lrwork))
    end select

    write(6,'(a,3i10)') 'lwork,lrwork,liwork=',lwork,lrwork,liwork ! for test
    allocate(iwork(liwork))        !, stat=status)

    select type(A)
      type is(real)
      select type(work)
        type is(real)
        select type(w)
          type is(real)
          call glapack_ev('V', 'L', n, A, lda, W, work, lwork, iwork, liwork, info)
        end select
      end select
      type is(real(8))
      select type(work)
        type is(real(8))
        select type(w)
          type is(real(8))
          call glapack_ev('V', 'L', n, A, lda, W, work, lwork, iwork, liwork, info)
        end select
      end select
      type is(complex)
      select type(work)
        type is(complex)
        select type(w)
          type is(real)
          select type(rwork)
            type is(real)
            call glapack_ev('V', 'L', n, A, lda, W, work, lwork, rwork, lrwork, iwork, liwork, info)
          end select
        end select
      end select
      type is(complex(8))
      select type(work)
        type is(complex(8))
        select type(w)
          type is(real(8))
          select type(rwork)
            type is(real(8))
            call glapack_ev('V', 'L', n, A, lda, W, work, lwork, rwork, lrwork, iwork, liwork, info)
          end select
        end select
      end select
    end select
    return
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
    real(8), allocatable :: work(:)
    integer, allocatable :: iwork(:)
    integer(4) :: lwork = -1
    integer :: ngpu = 1
    integer :: liwork != -1
    integer :: lrwork
    integer :: status

    call glapack_ev_prm('V', 'L', n, A, lda, lwork, lrwork, liwork) ! error occures here?
    write(6,'(a,3i10)') 'lwork,lrwork,liwork=',lwork,lrwork,liwork ! for test
    allocate (work(lwork), iwork(liwork))        !, stat=status)
    call glapack_ev('V', 'L', n, A, lda, W, work, lwork, iwork, liwork, info)

    return
  end subroutine

end program
