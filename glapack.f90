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

module glapack
#ifdef MAGMA
  use magma
  use cublasXt
#endif
  implicit none

  integer, private :: ngpu

  interface glapack_init
    module procedure glapack_init0,glapack_init1
  end interface

  interface glapack_ev_prm
    !module procedure glapack_gev_prm
    module procedure glapack_ssyev_prm,glapack_dsyev_prm,glapack_cheev_prm,glapack_zheev_prm
  end interface

  interface glapack_gemm
    module procedure glapack_sgemm,glapack_dgemm,glapack_cgemm,glapack_zgemm
    !module procedure glapack_ggemm
  end interface

  interface glapack_potri
    module procedure glapack_spotri,glapack_dpotri,glapack_cpotri,glapack_zpotri
    !module procedure glapack_gpotri
  end interface

  interface glapack_tri2x
    module procedure glapack_ssytri2x,glapack_dsytri2x,glapack_chetri2x,glapack_zhetri2x
    !module procedure glapack_gpotri
  end interface

  interface glapack_potrf
    module procedure glapack_spotrf,glapack_dpotrf,glapack_cpotrf,glapack_zpotrf
    !module procedure glapack_gpotrf
  end interface

  interface glapack_trf
    module procedure glapack_ssytrf,glapack_dsytrf,glapack_chetrf,glapack_zhetrf
    !module procedure glapack_gpotrf
  end interface

  interface glapack_ev
    module procedure glapack_ssyev,glapack_dsyev,glapack_cheev,glapack_zheev,glapack_ssyev0
    !module procedure  glapack_gev
  end interface

  interface glapack_gst
    module procedure glapack_ssygst,glapack_dsygst,glapack_chegst,glapack_zhegst
    !module procedure glapack_ggst
  end interface

  interface glapack_qp3
    module procedure glapack_sgeqp3,glapack_dgeqp3,glapack_cgeqp3,glapack_zgeqp3
    !module procedure glapack_gqp3
  end interface

  interface glapack_gqr
    module procedure glapack_sorgqr,glapack_dorgqr,glapack_cungqr,glapack_zungqr
    !module procedure glapack_ggqr
  end interface

  interface glapack_trsm
    module procedure glapack_strsm,glapack_dtrsm,glapack_ctrsm,glapack_ztrsm
    !module procedure glapack_gtrsm
  end interface

  private wt_info
  private glapack_getTrans,glapack_getSide, glapack_getUplo, glapack_getDiag

contains

  subroutine glapack_init0()
    ngpu=1
#ifdef MAGMA
    call magmaf_init()
#endif
  end subroutine

  subroutine glapack_init1(ngpu_)
    integer ngpu_
    ngpu=ngpu_
#ifdef MAGMA
    call magmaf_init()
#endif
  end subroutine

  integer function glapack_getTrans(trans)
    character trans   ! 'N', 'T' or 'C'
    if(trans=='N'.or.trans=='n') then
      glapack_getTrans=0
    else if(trans=='T'.or.trans=='t') then
      glapack_getTrans=1
    else if(trans=='C'.or.trans=='c') then
      glapack_getTrans=2
    end if
  end function

  integer function glapack_getUplo(uplo)
    character uplo   ! 'L' or 'U'
    if(uplo=='L'.or.uplo=='l') then
      glapack_getUplo=0
    else if(uplo=='U'.or.uplo=='u') then
      glapack_getUplo=1
    end if
  end function

  integer function glapack_getSide(side)
    character side  ! 'L' or 'R'
    if(side=='L'.or.side=='l') then
      glapack_getSide=0
    else if(side=='R'.or.side=='r') then
      glapack_getSide=1
    end if
  end function

  integer function glapack_getDiag(diag)
    character diag   ! 'N' or 'U'
    if(diag=='N'.or.diag=='n') then
      glapack_getDiag=0
    else if(diag=='U'.or.diag=='u') then
      glapack_getDiag=1
    end if
  end function

  !  subroutine glapack_set_classp(A,sap,dap,cap,zap)
  !    class(*),target:: A
  !    real,pointer:: sap
  !    real(8),pointer:: dap
  !    complex,pointer :: cap
  !    complex(8),pointer :: zap
  !    select type (A)
  !      type is (real)
  !      sap=>A
  !      type is (real(8))
  !      dap=>A
  !      type is (complex)
  !      cap=>A
  !      type is (complex(8))
  !      zap=>A
  !    end select
  !  end subroutine
  !
  !  subroutine glapack_set_class1dp(A,lda,sap,dap,cap,zap)
  !    class(*),target:: A(lda)
  !    integer:: lda
  !    real,pointer:: sap(:)
  !    real(8),pointer:: dap(:)
  !    complex,pointer :: cap(:)
  !    complex(8),pointer :: zap(:)
  !    select type (A)
  !      type is (real)
  !      sap=>A
  !      type is (real(8))
  !      dap=>A
  !      type is (complex)
  !      cap=>A
  !      type is (complex(8))
  !      zap=>A
  !    end select
  !  end subroutine
  !
  !  subroutine glapack_set_class1drp(A,lda,sap,dap)
  !    class(*),target:: A(lda)
  !    integer:: lda
  !    real,pointer:: sap(:)
  !    real(8),pointer:: dap(:)
  !    select type (A)
  !      type is (real)
  !      sap=>A
  !      type is (real(8))
  !      dap=>A
  !    end select
  !  end subroutine
  !
  !  subroutine glapack_set_class2dp(A,lda,k,sap,dap,cap,zap)
  !    class(*),target:: A(lda,k)
  !    integer:: lda,k
  !    real,pointer:: sap(:,:)
  !    real(8),pointer:: dap(:,:)
  !    complex,pointer :: cap(:,:)
  !    complex(8),pointer :: zap(:,:)
  !    select type (A)
  !      type is (real)
  !      sap=>A
  !      type is (real(8))
  !      dap=>A
  !      type is (complex)
  !      cap=>A
  !      type is (complex(8))
  !      zap=>A
  !    end select
  !  end subroutine
  !
  !  subroutine glapack_ggemm( transA, transB, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
  !    dC, lddc )
  !    character           :: transA
  !    character           :: transB
  !    integer               :: m
  !    integer               :: n
  !    integer               :: k
  !    class(*)   :: alpha
  !    class(*)   :: dA(ldda,n)   ! magma_devptr_t  :: dA
  !    integer                :: ldda
  !    class(*)   :: dB(lddb,n)   !magma_devptr_t   :: dB
  !    integer                :: lddb
  !    class(*)   :: beta
  !    class(*)   :: dC(lddc,n)   !magma_devptr_t   :: dC
  !    integer                :: lddc
  !    real,pointer :: sdap(:,:),sdbp(:,:),sdcp(:,:),salphap,sbetap
  !    real(8),pointer :: ddap(:,:),ddbp(:,:),ddcp(:,:),dalphap,dbetap
  !    complex,pointer ::cdap(:,:),cdbp(:,:),cdcp(:,:),calphap,cbetap
  !    complex(8),pointer :: zdap(:,:),zdbp(:,:),zdcp(:,:),zalphap,zbetap
  !    call glapack_set_class2dp(dA,ldda,k,sdap,ddap,cdap,zdap)
  !    call glapack_set_class2dp(dB,lddb,n,sdbp,ddbp,cdbp,zdbp)
  !    call glapack_set_class2dp(dC,lddc,n,sdcp,ddcp,cdcp,zdcp)
  !    call glapack_set_classp(alpha,salphap,dalphap,calphap,zalphap)
  !    call glapack_set_classp(beta,sbetap,dbetap,cbetap,zbetap)
  !
  !    select type (dA)
  !      type is (real)
  !      call glapack_sgemm( transA, transB, m, n, k, salphap, sdAp, ldda, sdBp, lddb, sbetap,  &
  !        sdCp, lddc )
  !      type is (real(8))
  !      call glapack_dgemm( transA, transB, m, n, k, dalphap, ddap, ldda, ddbp, lddb, dbetap,  &
  !        ddcp, lddc )
  !      type is (complex)
  !      call glapack_cgemm( transA, transB, m, n, k, calphap, cdap, ldda, cdbp, lddb, cbetap,  &
  !        cdcp, lddc )
  !      type is (complex(8))
  !      call glapack_zgemm( transA, transB, m, n, k, zalphap, zdap, ldda, zdbp, lddb, zbetap,  &
  !        zdcp, lddc )
  !    end select
  !
  !  end subroutine

  subroutine glapack_sgemm( transA, transB, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
    dC, lddc )
    character           :: transA
    character           :: transB
    integer             :: m
    integer             :: n
    integer             :: k
    real                :: alpha
    real :: dA(ldda,*)   ! magma_devptr_t  :: dA
    integer             :: ldda
    real  :: dB(lddb,*)   !magma_devptr_t   :: dB
    integer             :: lddb
    real                :: beta
    real  :: dC(lddc,*)   !magma_devptr_t   :: dC
    integer             :: lddc
    integer:: tra,trb,istat
#ifdef MAGMA
    !external magmablas_sgemm
    !#    ifdef CUDAFOR
    integer :: h1
    call cublasXt_Create(h1)
    !#    endif
#endif


#ifdef MAGMA
    !#   ifdef CUDAFOR
    tra=glapack_getTrans(transA); trb=glapack_getTrans(transB)
    istat=cublasXt_Sgemm(h1, tra, trb, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
      dC, lddc )
    !#    else
    !istat = cublasXt_sgemm(h1, transA, transB, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
    !  dC, lddc )
    !#    endif
#else
    call sgemm( transA, transB, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
      dC, lddc )
#endif
  end subroutine

  subroutine glapack_cgemm( transA, transB, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
    dC, lddc )
    character        :: transA
    character        :: transB
    integer            :: m
    integer            :: n
    integer            :: k
    complex             :: alpha
    complex :: dA(ldda,*)    ! magma_devptr_t  :: dA
    integer          :: ldda
    complex :: dB(lddb,*)    ! magma_devptr_t   :: dB
    integer             :: lddb
    complex             :: beta
    complex   :: dC(lddc,*)   !magma_devptr_t   :: dC
    integer             :: lddc
    integer:: tra,trb,istat
#ifdef MAGMA
    !external magmablas_cgemm
    !#ifdef CUDAFOR
    integer :: h1
    call cublasXt_Create(h1)
    !#endif
#endif
#ifdef MAGMA
    !#   ifdef CUDAFOR
    tra=glapack_getTrans(transA); trb=glapack_getTrans(transB)
    istat=cublasXt_Cgemm(h1, tra, trb, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
      dC, lddc )
    !#   else
    !istat = cublasXt_cgemm(h1, transA, transB, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
    !  dC, lddc )
    !#   endif
#else
    call cgemm( transA, transB, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
      dC, lddc )
#endif
  end subroutine

  subroutine glapack_dgemm( transA, transB, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
    dC, lddc )
    character        :: transA
    character        :: transB
    integer          :: m
    integer          :: n
    integer          :: k
    real(8)             :: alpha
    real(8) :: dA(ldda,*)   !magma_devptr_t  :: dA
    integer          :: ldda
    real(8) :: dB(lddb,*)   !magma_devptr_t   :: dB
    integer          :: lddb
    real(8)             :: beta
    real(8)  :: dC(lddc,*)   !magma_devptr_t   :: dC
    integer          :: lddc
    integer:: tra,trb,istat
#ifdef MAGMA
    !external magmablas_dgemm
    !#    ifdef CUDAFOR
    integer :: h1
    call cublasXt_Create(h1)
    !#    endif
#endif

#ifdef MAGMA
    !#   ifdef CUDAFOR
    tra=glapack_getTrans(transA); trb=glapack_getTrans(transB)
    istat=cublasXt_Dgemm(h1, tra, trb, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
      dC, lddc )
    !#    else
    !istat = cublasXt_dgemm(h1, transA, transB, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
    !  dC, lddc )
    !#    endif
#else
    call dgemm( transA, transB, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
      dC, lddc )

#endif
  end subroutine

  subroutine glapack_zgemm( transA, transB, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
    dC, lddc )
    character        :: transA
    character        :: transB
    integer          :: m
    integer          :: n
    integer          :: k
    complex(8)  :: alpha
    complex(8)   :: dA(ldda,*)   !magma_devptr_t  :: dA
    integer          :: ldda
    complex(8)  :: dB(lddb,*)   !magma_devptr_t   :: dB
    integer          :: lddb
    complex(8)   :: beta
    complex(8)   :: dC(lddc,*)   !magma_devptr_t  :: dC
    integer          :: lddc
    integer:: tra,trb,istat
#ifdef MAGMA
    !external magmablas_zgemm
    !#    ifdef CUDAFOR
    integer :: h1
    call cublasXt_Create(h1)
    !#    endif
#endif

#ifdef MAGMA
    !#   ifdef CUDAFOR
    tra=glapack_getTrans(transA); trb=glapack_getTrans(transB)
    istat=cublasXt_Zgemm(h1, tra, trb, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
      dC, lddc )
    !#    else
    !istat = cublasXt_zgemm(h1, transA, transB, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
    !  dC, lddc )
    !#    endif
#else
    call zgemm( transA, transB, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
      dC, lddc )
#endif
  end subroutine

  !  subroutine glapack_gpotri( uplo, n, A, lda, info )
  !    character        :: uplo
  !    integer          :: n
  !    class(*),target :: A(lda,n)
  !    integer          :: lda
  !    integer          :: info
  !
  !#ifdef MAGMA
  !    !call magmaf_spotri( uplo, n, reshape(A,asize), lda, info )
  !    !call magmaf_spotri( uplo, n, A, lda, info )
  !    select type (A)
  !      type is (real)
  !      call magmaf_spotri( uplo, n, A, lda, info )
  !      type is (real(8))
  !      call magmaf_dpotri( uplo, n, A, lda, info )
  !      type is (complex)
  !      call magmaf_cpotri( uplo, n, A, lda, info )
  !      type is (complex(8))
  !      call magmaf_zpotri( uplo, n, A, lda, info )
  !    end select
  !#else
  !    select type (A)
  !      type is (real)
  !      call spotri( uplo, n, A, lda, info )
  !      type is (real(8))
  !      call dpotri( uplo, n, A, lda, info )
  !      type is (complex)
  !      call cpotri( uplo, n, A, lda, info )
  !      type is (complex(8))
  !      call zpotri( uplo, n, A, lda, info )
  !    end select
  !#endif
  !  end subroutine

  subroutine glapack_spotri( uplo, n, A, lda, info )
    character        :: uplo
    integer          :: n
    real             :: A(lda,n)
    integer          :: lda
    integer          :: info
    !integer :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_spotri( uplo, n, reshape(A,asize), lda, info )
    call magmaf_spotri( uplo, n, A, lda, info )
#else
    call spotri( uplo, n, A, lda, info )
#endif
  end subroutine

  subroutine glapack_cpotri( uplo, n, A, lda, info )
    character        :: uplo
    integer          :: n
    complex          :: A(lda,n)
    integer          :: lda
    integer          :: info
    !integer :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_cpotri( uplo, n, reshape(A,asize), lda, info )
    call magmaf_cpotri( uplo, n, A, lda, info )
#else
    call cpotri( uplo, n, A, lda, info )
#endif
  end subroutine

  subroutine glapack_dpotri( uplo, n, A, lda, info )
    character        :: uplo
    integer          :: n
    real(8)           :: A(lda,n)
    integer          :: lda
    integer          :: info
    !integer :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_dpotri( uplo, n, reshape(A,asize), lda, info )
    call magmaf_dpotri( uplo, n, A, lda, info )
#else
    call dpotri( uplo, n, A, lda, info )
#endif
  end subroutine

  subroutine glapack_zpotri( uplo, n, A, lda, info )
    character        :: uplo
    integer          :: n
    complex(8)    :: A(lda,n)
    integer          :: lda
    integer          :: info
    !integer :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_zpotri( uplo, n, reshape(A,asize), lda, info )
    call magmaf_zpotri( uplo, n, A, lda, info )
#else
    call zpotri( uplo, n,A, lda, info )
#endif
  end subroutine

  subroutine glapack_ssytrf( uplo, n, A, lda, ipiv, work, lwork, info )
    character        :: uplo
    integer          :: n
    real(4)          :: A(lda,n)
    integer          :: lda
    integer          :: ipiv(n)
    real(4)          :: work(lwork)
    integer          :: lwork
    integer          :: info
#ifdef MAGMA
    call magmaf_ssytrf( uplo, n, A, lda, ipiv, info )
#else
    call ssytrf( uplo, n, A, lda, ipiv, work, lwork, info )
#endif
  end subroutine

  subroutine glapack_dsytrf( uplo, n, A, lda, ipiv, work, lwork, info )
    character        :: uplo
    integer          :: n
    real(8)          :: A(lda,n)
    integer          :: lda
    integer          :: ipiv(n)
    real(8)          :: work(lwork)
    integer          :: lwork
    integer          :: info
#ifdef MAGMA
    call magmaf_dsytrf( uplo, n, A, lda, ipiv, info )
#else
    call dsytrf( uplo, n, A, lda, ipiv, work, lwork, info )
#endif
  end subroutine


  subroutine glapack_chetrf( uplo, n, A, lda, ipiv, work, lwork, info )
    character        :: uplo
    integer          :: n
    complex(4)       :: A(lda,n)
    integer          :: lda
    integer          :: ipiv(n)
    complex(4)       :: work(lwork)
    integer          :: lwork
    integer          :: info
#ifdef MAGMA
    call magmaf_chetrf( uplo, n, A, lda, ipiv, info )
#else
    call chetrf( uplo, n, A, lda, ipiv, work, lwork, info )
#endif
  end subroutine

  subroutine glapack_zhetrf( uplo, n, A, lda, ipiv, work,lwork,info )
    character        :: uplo
    integer          :: n
    complex(8)       :: A(lda,n)
    integer          :: lda
    integer          :: ipiv(n)
    complex(8)       :: work(lwork)
    integer          :: lwork
    integer          :: info
#ifdef MAGMA
    call magmaf_zhetrf( uplo, n, A, lda, ipiv, info )
#else
    call zhetrf( uplo, n, A, lda, ipiv, work, lwork, info )
#endif
  end subroutine


  !  subroutine glapack_gpotrf( uplo, n, A, lda, info )
  !    character        :: uplo
  !    integer          :: n
  !    class(*):: A(lda,n)  !class(*),target :: A(lda,n)
  !    integer          :: lda
  !    integer          :: info
  !
  !#ifdef MAGMA
  !    !call magmaf_spotrf( uplo, n, reshape(A,asize), lda, info )
  !    select type (A)
  !      type is (real)
  !      call magmaf_spotrf( uplo, n, A, lda, info )
  !      type is (real(8))
  !      call magmaf_dpotrf( uplo, n, A, lda, info )
  !      type is (complex)
  !      call magmaf_cpotrf( uplo, n, A, lda, info )
  !      type is (complex(8))
  !      call magmaf_zpotrf( uplo, n, A, lda, info )
  !    end select
  !#else
  !    !call spotrf( uplo, n, A, lda, info )
  !    select type (A)
  !      type is (real)
  !      call spotrf( uplo, n, A, lda, info )
  !      type is (real(8))
  !      call dpotrf( uplo, n, A, lda, info )
  !      type is (complex)
  !      call cpotrf( uplo, n, A, lda, info )
  !      type is (complex(8))
  !      call zpotrf( uplo, n, A, lda, info )
  !    end select
  !#endif
  !  end subroutine

  subroutine glapack_spotrf( uplo, n, A, lda, info )
    character        :: uplo
    integer          :: n
    real :: A(lda,n)
    integer          :: lda
    integer          :: info
    !integer :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_spotrf( uplo, n, reshape(A,asize), lda, info )
    call magmaf_spotrf( uplo, n, A, lda, info )
#else
    call spotrf( uplo, n, A, lda, info )
#endif
  end subroutine

  subroutine glapack_cpotrf( uplo, n, A, lda, info )
    character        :: uplo
    integer          :: n
    complex        :: A(lda,n)
    integer          :: lda
    integer          :: info
    !integer :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_cpotrf( uplo, n, reshape(A,asize), lda, info )
    call magmaf_cpotrf( uplo, n, A, lda, info )
#else
    call cpotrf( uplo, n, A, lda, info )
#endif
  end subroutine

  subroutine glapack_dpotrf( uplo, n, A, lda, info )
    character        :: uplo
    integer          :: n
    real(8)            :: A(lda,n)
    integer          :: lda
    integer          :: info
    !integer :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_dpotrf( uplo, n, reshape(A,asize), lda, info )
    call magmaf_dpotrf( uplo, n, A, lda, info )
#else
    call dpotrf( uplo, n, A, lda, info )
#endif
  end subroutine

  subroutine glapack_zpotrf( uplo, n, A, lda, info )
    character      :: uplo
    integer          :: n
    complex(8)    :: A(lda,n)
    integer          :: lda
    integer          :: info
    !integer :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_zpotrf( uplo, n, reshape(A,asize), lda, info )
    call magmaf_zpotrf( uplo, n, A, lda, info )
#else
    call zpotrf( uplo, n, A, lda, info )
#endif
  end subroutine

  subroutine glapack_ssytri2x(uplo,n,A,lda,ipiv,work,lwork,info)
    character        :: uplo
    integer          :: n
    real(4)          :: A(lda,n)
    integer          :: lda
    integer          :: ipiv(n)
    real(4)          :: work(lwork)
    integer          :: lwork
    integer          :: info
#ifdef MAGMA
    call magmaf_ssytri2x( uplo, n, A, lda, ipiv, work, lwork, info )
#else
    call ssytri2x( uplo, n, A, lda, ipiv, work, lwork, info )
#endif
  end subroutine

  subroutine glapack_dsytri2x(uplo,n,A,lda,ipiv,work,lwork,info)
    character        :: uplo
    integer          :: n
    real(8)          :: A(lda,n)
    integer          :: lda
    integer          :: ipiv(n)
    real(8)          :: work(lwork)
    integer          :: lwork
    integer          :: info
#ifdef MAGMA
    call magmaf_dsytri2x( uplo, n, A, lda, ipiv, work, lwork, info )
#else
    call dsytri2x( uplo, n, A, lda, ipiv, work, lwork, info )
#endif
  end subroutine

  subroutine glapack_chetri2x(uplo,n,A,lda,ipiv,work,lwork,info)
    character        :: uplo
    integer          :: n
    complex(4)       :: A(lda,n)
    integer          :: lda
    integer          :: ipiv(n)
    complex(4)       :: work(lwork)
    integer          :: lwork
    integer          :: info
#ifdef MAGMA
    call magmaf_chetri2x( uplo, n, A, lda, ipiv, work, lwork, info )
#else
    call chetri2x( uplo, n, A, lda, ipiv, work, lwork, info )
#endif
  end subroutine

  subroutine glapack_zhetri2x(uplo,n,A,lda,ipiv,work,lwork,info)
    character        :: uplo
    integer          :: n
    complex(8)       :: A(lda,n)
    integer          :: lda
    integer          :: ipiv(n)
    complex(8)       :: work(lwork)
    integer          :: lwork
    integer          :: info
#ifdef MAGMA
    call magmaf_zhetri2x( uplo, n, A, lda, ipiv, work, lwork, info )
#else
    call zhetri2x( uplo, n, A, lda, ipiv, work, lwork, info )
#endif
  end subroutine


  subroutine glapack_ssyev0( jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info )
    character        :: jobz
    character        :: uplo
    integer          :: n
    real             :: A(lda,n)
    integer          :: lda
    real             :: w(*)
    real             :: work(*)
    integer          :: lwork
    real             :: rwork(*)    ! dummy
    integer          :: lrwork      ! dummy
    integer          :: iwork(*)
    integer          :: liwork
    integer          :: info
    !integer :: asize(1)
    !asize(1)=size(A)
    call glapack_ssyev( jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info )
  end subroutine

  subroutine glapack_ssyev_prm(jobz,uplo,n,A,lda,lwork,lrwork,liwork)
    character       :: jobz
    character       :: uplo
    integer :: n
    real,target:: A(lda,n)
    integer :: lda
    integer :: lwork
    integer :: lrwork
    integer :: liwork
    real :: work(1)   !real,allocatable:: work(:)
    real :: rwork(1)  !real,allocatable:: rwork(:)
    integer :: iwork(1)  !integer,allocatable :: iwork(:)
    real,allocatable :: w(:)
    integer :: info
    allocate(w(n))
    lwork  = -1; liwork=-1; lrwork=-1  ! calculate work array size
    ! calculate work array size used in ev
    call glapack_ssyev( jobz, uplo, n, A, lda, w, work, lwork,iwork,liwork,info )
    if(info/=0) then
      write(6,*) 'info in glapack_ssyev ',info; stop
    end if 
    lwork = int(work(1)); liwork=iwork(1)  ! added
    lwork=max(1,lwork)
    lrwork=1
    liwork=max(1,liwork)
    !print *,'lwork=',lwork,' liwork=',liwork; stop  ! for test
  end subroutine

  subroutine glapack_dsyev_prm(jobz,uplo,n,A,lda,lwork,lrwork,liwork)
    character       :: jobz
    character       :: uplo
    integer :: n
    real(8),target:: A(lda,n)
    integer :: lda
    integer :: lwork
    integer :: lrwork
    integer :: liwork
    real(8):: work(1)   !    real(8),allocatable:: work(1)
    real(8):: rwork(1)  !    real(8),allocatable:: rwork(1)
    real(8),allocatable :: w(:)      !    real(8),allocatable:: w(:)
    integer ::iwork(1)  !    integer,allocatable :: iwork(:)
    integer :: info
    allocate(w(n))
    lwork  = -1; liwork=-1; lrwork=-1  ! calculate work array size
    ! calculate work array size used in ev
    call glapack_dsyev( jobz, uplo, n, A, lda, w, work, lwork,iwork,liwork,info )
    if(info/=0) then
      write(6,*) 'info in glapack dsyev ',info; stop
    end if
    lwork = int(work(1)); liwork=iwork(1)  ! added
    lwork=max(1,lwork)
    lrwork=1
    liwork=max(1,liwork)
    !print *,'lwork=',lwork,' liwork=',liwork; stop  ! for test
  end subroutine

  subroutine glapack_cheev_prm(jobz,uplo,n,A,lda,lwork,lrwork,liwork)
    character       :: jobz
    character       :: uplo
    integer :: n
    complex,target:: A(lda,n)
    integer :: lda
    integer :: lwork
    integer :: lrwork
    integer :: liwork
    complex:: work(1)  !complex,allocatable:: work(:)
    real:: rwork(1)    !real,allocatable:: rwork(:)
    real,allocatable :: w(:) !real,allocatable:: w(:)
    integer:: iwork(1) !integer,allocatable :: iwork(:)
    integer :: info
    allocate(w(n))
    lwork  = -1; liwork=-1; lrwork=-1  ! calculate work array size
    ! calculate work array size used in ev
    call glapack_cheev( jobz, uplo, n, A, lda, w, work, lwork,rwork,lrwork,iwork,liwork,info )
    if(info/=0) then
      write(6,*) 'info in glapack cheev ',info; stop
    end if
    lwork = int(work(1)); lrwork=int(rwork(1)); liwork=iwork(1)  ! added
    lwork=max(1,lwork)
    lrwork=max(1,lrwork)
    liwork=max(1,liwork)
    !print *,'lwork=',lwork,' liwork=',liwork; stop  ! for test
  end subroutine

  subroutine glapack_zheev_prm(jobz,uplo,n,A,lda,lwork,lrwork,liwork)
    character       :: jobz
    character       :: uplo
    integer :: n
    complex(8),target:: A(lda,n)
    integer :: lda
    integer :: lwork
    integer :: lrwork
    integer :: liwork
    complex(8) :: work(1)  !complex(8),allocatable:: work(:)
    real(8):: rwork(1)    !real(8),allocatable:: rwork(:)
    real(8),allocatable :: w(:)  !real(8),allocatable:: w(:)
    integer:: iwork(1)    !integer,allocatable :: iwork(:)
    integer :: info
    allocate(w(n))
    lwork  = -1; liwork=-1; lrwork=-1  ! calculate work array size
    ! calculate work array size used in ev
    call glapack_zheev( jobz, uplo, n, A, lda, w, work, lwork,rwork,lrwork,iwork,liwork,info )
    if(info/=0) then
      write(6,*) 'info in glapack zheev ',info; stop
    end if
    lwork = int(work(1)); lrwork=int(rwork(1)); liwork=iwork(1)  ! added
    lwork=max(1,lwork)
    lrwork=max(1,lrwork)
    liwork=max(1,liwork)
    !print *,'lwork=',lwork,' liwork=',liwork; stop  ! for test
  end subroutine

  !  ! this estimate waork array size lwork, liwork
  !  subroutine glapack_gev_prm(jobz,uplo,n,A,lda,lwork,lrwork,liwork)
  !    character       :: jobz
  !    character       :: uplo
  !    integer :: n
  !    class(*),target:: A(lda,n)
  !    integer :: lda
  !    integer :: lwork
  !    integer :: lrwork
  !    integer :: liwork
  !    !integer :: lwmax
  !    integer:: iwork(1)
  !    !real:: rtmp(1)
  !    class(*),allocatable:: work(:)
  !    class(*),allocatable:: rwork(:)
  !    class(*),allocatable:: w(:)
  !    integer :: info
  !            !allocate(work(1),rwork(1),w(1))
  !    select type (A)
  !      type is (real)
  !      allocate(real::work(1),rwork(1),w(1))
  !      type is (real(8))
  !      allocate(real(8)::work(1),rwork(1),w(1))
  !      type is (complex)
  !      allocate(complex::work(1))
  !      allocate(real::rwork(1),w(1))
  !      type is (complex(8))
  !      allocate(complex(8)::work(1))
  !      allocate(real(8)::rwork(1),w(1))
  !    end select
  !    lwork  = -1
  !    lrwork = -1
  !    liwork = -1
  !    !lrwork=1
  !    ! calculate work array size used in ev
  !    call glapack_gev( jobz, uplo, n, A, lda, w, work, lwork,rwork,lrwork,iwork, liwork, info )
  !    if(info/=0) call wt_info('magmaf_ssyevd',info)
  !    select type (work)
  !      type is (real)
  !      lwork = work(1)
  !      type is (real(8))
  !      lwork = work(1)
  !      type is (complex)
  !      lwork = work(1)
  !      type is (complex(8))
  !      lwork = work(1)
  !    end select
  !    select type (rwork)
  !      type is (real)
  !      lrwork=rwork(1)
  !      type is (real(8))
  !      lrwork=rwork(1)
  !    end select
  !    liwork = iwork(1)
  !    lwork=max(1,lwork)
  !    lrwork=max(1,lrwork)
  !    liwork=max(1,liwork)
  !    !lworkt =1 + 6*nb + 2*nb**2
  !    !lwork = max(lwork,lworkt)    ! for strange lwork
  !    !write(6,'(a,4i15)') 'n,lwork,lrwork,liwork in glapack=',n,lwork,lrwork,liwork  ! for test
  !    if(info/=0) then
  !      write(6,'(a,3i10)') 'cannt estimate lapack_param, lwork,lrwork,liwork=',lwork,lrwork,liwork
  !      stop
  !    end if
  !  end subroutine

  !  ! calculate work array size
  !  subroutine glapack_gev( jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info )
  !    character        :: jobz
  !    character        :: uplo
  !    integer          :: n
  !    class(*),target:: A(lda,n)
  !    integer          :: lda
  !    ! use work array in this module
  !    class(*),target:: w(n)  !w(:)          ! real or real(8)
  !    class(*),target:: work(lwork)   !work(*)
  !    integer          :: lwork
  !    class(*),target:: rwork(lrwork)   ! only for complex
  !    integer          ::  lrwork       ! only for complex
  !    integer          :: iwork(liwork)
  !    integer          :: liwork
  !    integer          :: info
  !
  !    real,pointer:: sap(:,:),sworkp(:)
  !    real(8),pointer:: dap(:,:),dworkp(:)
  !    complex,pointer:: cap(:,:),cworkp(:)
  !    complex(8),pointer:: zap(:,:),zworkp(:)
  !    real,pointer:: swp(:)
  !    real(8),pointer:: dwp(:)
  !    real,pointer::srworkp(:)
  !    real(8),pointer::drworkp(:)
  !
  !    call glapack_set_class2dp(A,lda,n,sap,dap,cap,zap)
  !    call glapack_set_class1dp(work,1,sworkp,dworkp,cworkp,zworkp)
  !    call glapack_set_class1drp(w,1,swp,dwp)
  !    call glapack_set_class1drp(rwork,1,srworkp,drworkp)
  !
  !    select type (A)
  !      type is (real)
  !      !if(same_type_as(A,W).and.same_type_as(A,work)) then
  !      !    call glapack_zungqr( m, n, k, A, lda, tau, dt, nb, info )
  !      !end if
  !      call glapack_ssyev( jobz, uplo, n, sAp, lda, sWp, sworkp, lwork, iwork, liwork, info )
  !      type is (real(8))
  !      !if(same_type_as(A,W).and.same_type_as(A,work)) then
  !      !    call glapack_dsyev( jobz, uplo, n, A, lda, W, work, lwork, iwork, liwork, info )
  !      !end if
  !      call glapack_dsyev( jobz, uplo, n, dAp, lda, dWp, dworkp, lwork, iwork, liwork, info )
  !      type is (complex)
  !      !if(same_type_as(A,W).and.same_type_as(A,work)) then
  !      !    select type (rwork)
  !      !        type is (real)
  !      !        call glapack_cheev( jobz, uplo, n, A, lda, W, work, lwork, rwork, lrwork, iwork, liwork, info )
  !      !    end select
  !      call glapack_cheev( jobz, uplo, n, cAp, lda, sWp, cworkp, lwork, srworkp, lrwork, iwork, liwork, info )
  !      !end if
  !      type is (complex(8))
  !      !if(same_type_as(A,W).and.same_type_as(A,work)) then
  !      !    select type (rwork)
  !      !        type is (real(8))
  !      !        call glapack_zheev( jobz, uplo, n, A, lda, W, work, lwork, rwork, lrwork, iwork, liwork, info )
  !      !    end select
  !      call glapack_zheev( jobz, uplo, n, zAp, lda, dWp, zworkp, lwork, drworkp, lrwork, iwork, liwork, info )
  !        !end if
  !    end select
  !  end subroutine

  subroutine glapack_ssyev( jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info )
    character        :: jobz
    character        :: uplo
    integer          :: n
    real             :: A(lda,n)  ! real symmetric matrix -> eigenvector
    integer          :: lda
    real             :: w(*)      ! eigenvalue
    real             :: work(*)
    integer          :: lwork
    integer          :: iwork(*)
    integer          :: liwork
    integer          :: info
      !integer          :: ngpu=1
      !integer :: asize(1)
      !asize(1)=size(A)
#ifdef MAGMA
    if(ngpu==1) then
      call magmaf_ssyevd( jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info )
    else
      call magmaf_ssyevd_m(ngpu, jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info )
    end if
#else
    call ssyevd( jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info )
#endif
  end subroutine

  subroutine glapack_cheev( jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info )
    character        :: jobz
    character        :: uplo
    integer          :: n
    complex          :: A(lda,n)  ! Hermite matrix -> eigenvector
    integer          :: lda
    real             :: w(*)      ! eigenvalue
    complex          :: work(*)
    integer          :: lwork
    real             :: rwork(*)
    integer          :: lrwork
    integer          :: iwork(*)
    integer          :: liwork
    integer          :: info
    integer          :: ngpu = 1
      !integer          :: asize(1)
      !asize(1)=size(A)
#ifdef MAGMA
    if(ngpu==1) then
      call magmaf_cheevd( jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info )
    else
      call magmaf_cheevd_m(ngpu, jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info )
    end if
#else
    call cheevd( jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info )
#endif
  end subroutine

  subroutine glapack_dsyev( jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info )
    character        :: jobz  ! 'V' or 'N'
    character        :: uplo  ! 'L' or 'U'
    integer          :: n
    real(8)          :: A(lda,n)  ! real symmetric matrix -> eigenvector
    integer          :: lda
    real(8)          :: w(*)      ! eigenvalue
    real (8)         :: work(*)
    integer          :: lwork
    integer          :: iwork(*)
    integer          :: liwork
    integer          :: info
    integer          :: ngpu = 1
    !integer          :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    if(ngpu==1) then
      call magmaf_dsyevd( jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info )
    else
      call magmaf_dsyevd_m(ngpu, jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info )
    end if
#else
    call dsyevd( jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info )
    liwork=1
#endif
  end subroutine

  subroutine glapack_zheev( jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info )
    character        :: jobz  ! 'V' or 'N'
    character        :: uplo  ! 'L' or 'U'
    integer          :: n
    complex(8)       :: A(lda,n)  ! Hermite matrix -> eigenvector
    integer          :: lda
    real(8)          :: w(*)      ! eigenvalue
    complex(8)       :: work(*)
    integer          :: lwork
    real(8)          :: rwork(*)
    integer          :: lrwork
    integer          :: iwork(*)
    integer          :: liwork
    integer          :: info
    integer          :: ngpu = 1
     !integer          :: asize(1)
     !asize(1)=size(A)
#ifdef MAGMA
    if(ngpu==1) then
      call magmaf_zheevd( jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info )
    else
      call magmaf_zheevd_m(ngpu, jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info )
    end if
#else
    call zheevd( jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info )
#endif
  end subroutine

  !  subroutine glapack_ggqr( m, n, k, A, lda, tau, dT, nb, info )
  !    integer          :: m
  !    integer          :: n
  !    integer          :: k
  !    class(*),target:: A(lda,n)
  !    integer          :: lda
  !    class(*),target :: tau(:)
  !    class(*),target :: dT(:)   !magma_devptr_t   :: dT
  !    integer          :: nb
  !    integer          :: info
  !    real,pointer:: sap(:,:),staup(:),sdtp(:)
  !    real(8),pointer::dap(:,:),dtaup(:),ddtp(:)
  !    complex,pointer::cap(:,:),ctaup(:),cdtp(:)
  !    complex(8),pointer::zap(:,:),ztaup(:),zdtp(:)
  !    call glapack_set_class2dp(A,lda,n,sap,dap,cap,zap)
  !    call glapack_set_class1dp(tau,1,staup,dtaup,ctaup,ztaup)
  !    call glapack_set_class1dp(dt,1,sdtp,ddtp,cdtp,zdtp)
  !    select type (A)
  !      type is (real)
  !      !if(same_type_as(A,tau).and.same_type_as(A,dt)) then
  !      !    call glapack_sorgqr( m, n, k, A, lda, tau, dt, nb, info )
  !      !end if
  !      call glapack_sorgqr( m, n, k, sAp, lda, staup, sdtp, nb, info )
  !      type is (real(8))
  !      !if(same_type_as(A,tau).and.same_type_as(A,dt)) then
  !      !    call glapack_dorgqr( m, n, k, A, lda, tau, dt, nb, info )
  !      !end if
  !      call glapack_dorgqr( m, n, k, dAp, lda, dtaup, ddtp, nb, info )
  !      type is (complex)
  !      !if(same_type_as(A,tau).and.same_type_as(A,dt)) then
  !       !   call glapack_cungqr( m, n, k, A, lda, tau, dt, nb, info )
  !      !end if
  !      call glapack_cungqr( m, n, k, cAp, lda, ctaup, cdtp, nb, info )
  !      type is (complex(8))
  !      !if(same_type_as(A,tau).and.same_type_as(A,dt)) then
  !      !    call glapack_zungqr( m, n, k, A, lda, tau, dt, nb, info )
  !      !end if
  !      call glapack_zungqr( m, n, k, zAp, lda, ztaup, zdtp, nb, info )
  !    end select
  !
  !  end subroutine

  ! get Q matrix in QR decomposition
  subroutine glapack_sorgqr( m, n, k, A, lda, tau, dT, nb, info )
    integer          :: m
    integer          :: n
    integer          :: k
    real             :: A(lda,n)
    integer          :: lda
    real               :: tau(*)
    real    :: dT(*)   !magma_devptr_t   :: dT
    integer          :: nb
    integer          :: info
    !integer      :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_sorgqr_m( m, n, k, reshape(A,asize), lda, tau, dT, nb, info )
    call magmaf_sorgqr_m( m, n, k, A, lda, tau, dT, nb, info )
#else
    call sorgqr( m, n, k, A, lda, tau, dT, nb, info )
#endif
  end subroutine

  subroutine glapack_dorgqr( m, n, k, A, lda, tau, dT, nb, info )
    integer          :: m
    integer          :: n
    integer          :: k
    real(8)             :: A(lda,n)
    integer          :: lda
    real(8)            :: tau(*)
    real(8)   :: dT(*)   !  magma_devptr_t   :: dT
    integer          :: nb
    integer          :: info
    !integer      ::asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_dorgqr_m( m, n, k, reshape(A,asize), lda, tau, dT, nb, info )
    call magmaf_dorgqr_m( m, n, k, A, lda, tau, dT, nb, info )
#else
    call dorgqr( m, n, k, A, lda, tau, dT, nb, info )
#endif
  end subroutine

  subroutine glapack_cungqr( m, n, k, A, lda, tau, dT, nb, info )
    integer          :: m
    integer          :: n
    integer          :: k
    complex        :: A(lda,n)
    integer          :: lda
    complex        :: tau(*)
    complex  :: dT(*)   !magma_devptr_t  :: dT
    integer          :: nb
    integer          :: info
    !integer     :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_cungqr_m( m, n, k, reshape(A,asize), lda, tau, dT, nb, info )
    call magmaf_cungqr_m( m, n, k, A, lda, tau, dT, nb, info )
#else
    call cungqr( m, n, k, A, lda, tau, dT, nb, info )
#endif
  end subroutine

  subroutine glapack_zungqr( m, n, k, A, lda, tau, dT, nb, info )
    integer          :: m
    integer          :: n
    integer          :: k
    complex(8)    :: A(lda,n)
    integer          :: lda
    complex(8)    :: tau(*)
    complex(8)  :: dT(*)   ! magma_devptr_t   :: dT
    integer          :: nb
    integer          :: info
    !integer     :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_zungqr_m( m, n, k, reshape(A,asize), lda, tau, dT, nb, info )
    call magmaf_zungqr_m( m, n, k, A, lda, tau, dT, nb, info )
#else
    call zungqr( m, n, k, A, lda, tau, dT, nb, info )
#endif
  end subroutine

  !  subroutine glapack_ggst( itype, uplo, n, A, lda, B, ldb, info )
  !    integer          :: itype
  !    character      :: uplo
  !    integer          :: n
  !    class(*),target:: A(lda,n)
  !    integer          :: lda
  !    class(*),target :: B(ldb,n)
  !    integer          :: ldb
  !    integer          :: info
  !    real,pointer:: sap(:,:),sbp(:,:)
  !    real(8),pointer::dap(:,:),dbp(:,:)
  !    complex,pointer::cap(:,:),cbp(:,:)
  !    complex(8),pointer::zap(:,:),zbp(:,:)
  !    call glapack_set_class2dp(A,lda,n,sap,dap,cap,zap)
  !    call glapack_set_class2dp(B,ldb,n,sbp,dbp,cbp,zbp)
  !    select type (A)
  !      type is (real)
  !      !if(same_type_as(A,B)) then
  !      !    call glapack_ssygst( itype, uplo, n,A, lda, B, ldb, info )
  !      !end if
  !      call glapack_ssygst( itype, uplo, n,sAp, lda, sBp, ldb, info )
  !      type is(real(8))
  !      !if(same_type_as(A,B)) then
  !      !    call glapack_dsygst( itype, uplo, n, A, lda, B, ldb, info )
  !      !end if
  !      call glapack_dsygst( itype, uplo, n, dAp, lda, dBp, ldb, info )
  !      type is (complex)
  !      !if(same_type_as(A,B)) then
  !      !    call glapack_chegst( itype, uplo, n, A, lda, B, ldb, info )
  !      !end if
  !      call glapack_chegst( itype, uplo, n, cAp, lda, cBp, ldb, info )
  !      type is (complex(8))
  !      !if(same_type_as(A,B)) then
  !      !    call  glapack_zhegst( itype, uplo, n, A, lda, B, ldb, info )
  !      !end if
  !      call  glapack_zhegst( itype, uplo, n, zAp, lda, zBp, ldb, info )
  !    end select
  !  end subroutine

  subroutine glapack_ssygst( itype, uplo, n, A, lda, B, ldb, info )
    integer          :: itype
    character      :: uplo
    integer          :: n
    real                :: A(lda,n)
    integer          :: lda
    real                :: B(ldb,n)
    integer          :: ldb
    integer          :: info
    !integer         ::asize(1)
    !integer         ::bsize(1)
    !asize(1)=size(A)
    !bsize(1)=size(B)
#ifdef MAGMA
    !call magmaf_ssygst( itype, uplo, n, reshape(A,asize), lda,reshape(B,bsize), ldb, info )
    call magmaf_ssygst( itype, uplo, n, A, lda, B, ldb, info )
#else
    call ssygst( itype, uplo, n, A, lda, B, ldb, info )
#endif
  end subroutine

  subroutine glapack_chegst( itype, uplo, n, A, lda, B, ldb, info )
    integer          :: itype
    character        :: uplo
    integer           :: n
    complex         :: A(lda,n)
    integer          :: lda
    complex        :: B(ldb,n)
    integer          :: ldb
    integer          :: info
    !integer         ::asize(1)
    !integer         ::bsize(1)
    !asize(1)=size(A)
    !bsize(1)=size(B)
#ifdef MAGMA
    !call magmaf_chegst( itype, uplo, n, reshape(A,asize), lda, reshape(B,bsize), ldb, info )
    call magmaf_chegst( itype, uplo, n, A, lda, B, ldb, info )
#else
    call chegst( itype, uplo, n, A, lda, B, ldb, info )
#endif
  end subroutine

  subroutine glapack_dsygst( itype, uplo, n, A, lda, B, ldb, info )
    integer          :: itype
    character        :: uplo
    integer          :: n
    real(8)            :: A(lda,n)
    integer          :: lda
    real(8)            :: B(ldb,n)
    integer          :: ldb
    integer          :: info
    !integer         ::asize(1)
    !integer         ::bsize(1)
    !asize(1)=size(A)
    !bsize(1)=size(B)
#ifdef MAGMA
    !call magmaf_dsygst( itype, uplo, n, reshape(A,asize), lda, reshape(B,bsize), ldb, info )
    call magmaf_dsygst( itype, uplo, n, A, lda, B, ldb, info )
#else
    call dsygst( itype, uplo, n, A, lda, B, ldb, info )
#endif
  end subroutine

  subroutine glapack_zhegst( itype, uplo, n, A, lda, B, ldb, info )
    integer          :: itype
    character       :: uplo
    integer          :: n
    complex(8)    :: A(lda,n)
    integer          :: lda
    complex(8)    :: B(ldb,n)
    integer          :: ldb
    integer          :: info
    !integer         ::asize(1)
    !integer         ::bsize(1)
    !asize(1)=size(A)
    !bsize(1)=size(B)
#ifdef MAGMA
    !call magmaf_zhegst( itype, uplo, n, reshape(A,asize), lda, reshape(B,bsize), ldb, info )
    call magmaf_zhegst( itype, uplo, n, A, lda, B, ldb, info )
#else
    call zhegst( itype, uplo, n, A, lda, B, ldb, info )
#endif
  end subroutine

  !  ! for column pivoting QR factorization
  !  subroutine glapack_gqp3( m, n, A, lda, jpvt, tau, work, lwork,rwork,info )
  !    integer          :: m
  !    integer          :: n
  !    class(*),target:: A(lda,n)
  !    integer          :: lda
  !    integer          :: jpvt(*)
  !    class(*),target:: tau(:)
  !    class(*) :: work(:)   ! use a type of work for swork etc
  !    class(*) :: rwork(:)  ! use a type of rwork for srwork etc
  !    integer          :: lwork
  !    integer          :: info
  !    real,pointer:: sap(:,:),staup(:),sworkp(:),srworkp(:)
  !    real(8),pointer::dap(:,:),dtaup(:),dworkp(:),drworkp(:)
  !    complex,pointer::cap(:,:),ctaup(:),cworkp(:)
  !    complex(8),pointer::zap(:,:),ztaup(:),zworkp(:)
  !    call glapack_set_class2dp(A,lda,n,sap,dap,cap,zap)
  !    call glapack_set_class1dp(tau,1,staup,dtaup,ctaup,ztaup)
  !    call glapack_set_class1dp(work,1,sworkp,dworkp,cworkp,zworkp)
  !    call glapack_set_class1drp(rwork,1,srworkp,drworkp)
  !    select type (A)
  !      type is (real)
  !      !if(same_type_as(A,tau).and.same_type_as(A,work)) then
  !      !    call glapack_sgeqp3( m, n, A, lda, jpvt, tau, work, lwork, info )
  !      !end if
  !      call glapack_sgeqp3( m, n, sAp, lda, jpvt, staup, sworkp, lwork, info )
  !      type is (real(8))
  !      !if(same_type_as(A,tau).and.same_type_as(A,work)) then
  !      !    call glapack_dgeqp3( m, n, A, lda, jpvt, tau, work, lwork, info )
  !      !end if
  !      call glapack_dgeqp3( m, n, dAp, lda, jpvt, dtaup, dworkp, lwork, info )
  !      type is (complex)
  !      !if(same_type_as(A,tau).and.same_type_as(A,work)) then
  !          !select type (rwork)
  !          !    type is (real)
  !          !    call glapack_cgeqp3( m, n, A, lda, jpvt, tau, work, lwork, rwork, info )
  !          !end select
  !      !end if
  !      call glapack_cgeqp3( m, n, cAp, lda, jpvt, ctaup, cworkp, lwork, srworkp, info )
  !      type is (complex(8))
  !      !if(same_type_as(A,tau).and.same_type_as(A,work)) then
  !      !    select type (rwork)
  !      !        type is (real(8))
  !       !       call glapack_zgeqp3( m, n, A, lda, jpvt, tau, work, lwork, rwork, info )
  !       !   end select
  !      !end if
  !      call glapack_zgeqp3( m, n, zAp, lda, jpvt, ztaup, zworkp, lwork, drworkp, info )
  !    end select
  !  end subroutine

  ! for column pivoting QR factorization
  subroutine glapack_sgeqp3( m, n, A, lda, jpvt, tau, work, lwork, info )
    integer          :: m
    integer          :: n
    real               :: A(lda,n)
    integer          :: lda
    integer          :: jpvt(*)
    real               :: tau(*)
    real               :: work(*)
    integer          :: lwork
    integer          :: info
    !integer     :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_sgeqp3( m, n, reshape(A,asize), lda, jpvt, tau, work, lwork, info )
    call magmaf_sgeqp3( m, n, A, lda, jpvt, tau, work, lwork, info )
#else
    call sgeqp3( m, n, A, lda, jpvt, tau, work, lwork, info )
#endif
  end subroutine

  ! for column pivoting QR factorization
  subroutine glapack_cgeqp3( m, n, A, lda, jpvt, tau, work, lwork, rwork, info )
    integer          :: m
    integer          :: n
    complex        :: A(lda,n)
    integer          :: lda
    integer          :: jpvt(*)
    complex        :: tau(*)
    complex        :: work(*)
    integer          :: lwork
    real        :: rwork(*)
    integer          :: info
    !integer     :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_cgeqp3( m, n, reshape(A,asize), lda, jpvt, tau, work, lwork,  rwork, info )
    call magmaf_cgeqp3( m, n, A, lda, jpvt, tau, work, lwork,  rwork, info )
#else
    call cgeqp3( m, n, A, lda, jpvt, tau, work, lwork, rwork, info )
#endif
  end subroutine

  ! for column pivoting QR factorization
  subroutine glapack_dgeqp3( m, n, A, lda, jpvt, tau, work, lwork, info )
    integer          :: m
    integer          :: n
    real(8)             :: A(lda,n)
    integer          :: lda
    integer          :: jpvt(*)
    real (8)            :: tau(*)
    real(8)             :: work(*)
    integer          :: lwork
    integer          :: info
    !integer     :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_dgeqp3( m, n, reshape(A,asize), lda, jpvt, tau, work, lwork, info )
    call magmaf_dgeqp3( m, n, A, lda, jpvt, tau, work, lwork, info )
#else
    call dgeqp3( m, n, A, lda, jpvt, tau, work, lwork, info )
#endif
  end subroutine

  ! for column pivoting QR factorization
  subroutine glapack_zgeqp3( m, n, A, lda, jpvt, tau, work, lwork, rwork, info )
    integer          :: m
    integer          :: n
    complex(8)  :: A(lda,n)
    integer          :: lda
    integer          :: jpvt(*)
    complex(8)  :: tau(*)
    complex(8)  :: work(*)
    integer          :: lwork
    real(8)            :: rwork(*)
    integer          :: info
    !integer     :: asize(1)
    !asize(1)=size(A)
#ifdef MAGMA
    !call magmaf_zgeqp3( m, n, reshape(A,asize), lda, jpvt, tau, work, lwork, rwork, info )
    call magmaf_zgeqp3( m, n, A, lda, jpvt, tau, work, lwork, rwork, info )
#else
    call zgeqp3( m, n, A, lda, jpvt, tau, work, lwork, rwork, info )
#endif
  end subroutine

  !    subroutine glapack_set_class2dp1(A,lda,k,sap,dap,cap,zap)
  !      class(*),target:: A(lda,k)
  !      integer:: lda,k
  !      real,pointer:: sap(:,:)
  !      real(8),pointer:: dap(:,:)
  !      complex,pointer:: cap(:,:)
  !      complex(8),pointer:: zap(:,:)
  !      call glapack_set_class2dp(A,lda,k,sap,dap,cap,zap)
  !    end subroutine

  !  ! this part does not work in gfortran any version including ver.7
  !  subroutine glapack_gtrsm (SIDE, UPLO,TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
  !    character   SIDE
  !    character   UPLO
  !    character   TRANSA
  !    character   DIAG
  !    integer     M
  !    integer     N
  !    class(*),target::ALPHA
  !    class(*),target :: A(LDA,M)   ! A(lda,k) k is m or n for side='L' or 'R'
  !    integer :: LDA
  !    class(*),target :: B(LDB,N)
  !    integer     LDB
  !    integer :: iside,iuplo,tra,idiag
  !    real,pointer:: sap(:,:),sbp(:,:),salphap
  !    real(8),pointer:: dap(:,:),dbp(:,:),dalphap
  !    complex,pointer:: cap(:,:),cbp(:,:),calphap
  !    complex(8),pointer:: zap(:,:),zbp(:,:),zalphap
  !    integer k
  !    if(SIDE=='L'.or.SIDE=='l') then
  !      k=M
  !    else if(SIDE=='R'.or.SIDE=='r') then
  !      k=N
  !        !write(6,'(a)') 'side=r is not supported in glapack'; stop
  !    end if
  !    call glapack_set_class2dp1(A,lda,k,sap,dap,cap,zap)
  !    call glapack_set_class2dp(B,ldb,n,sbp,dbp,cbp,zbp)
  !    call glapack_set_classp(alpha,salphap,dalphap,calphap,zalphap)
  !    select type (A)
  !      type is (real)
  !      !if(same_type_as(A,B).and.same_type_as(A,alpha)) then
  !      !    call glapack_strsm (SIDE, UPLO,TRANSA, DIAG, M, N, alpha, A, LDA, B, LDB)
  !      !end if
  !      call glapack_strsm (SIDE, UPLO,TRANSA, DIAG, M, N, salphap, sAp, LDA, sBp, LDB)
  !      type is (real(8))
  !      !if(same_type_as(A,B).and.same_type_as(A,alpha)) then
  !      !    call glapack_dtrsm (SIDE, UPLO,TRANSA, DIAG, M, N, alpha, A, LDA, B, LDB)
  !      !end if
  !      call glapack_dtrsm (SIDE, UPLO,TRANSA, DIAG, M, N, dalphap, dAp, LDA, dBp, LDB)
  !      type is (complex(4))
  !      !if(same_type_as(A,B).and.same_type_as(A,alpha)) then
  !      !    call glapack_ctrsm (SIDE, UPLO,TRANSA, DIAG, M, N, alpha, A, LDA, B, LDB)
  !      !end if
  !      call glapack_ctrsm (SIDE, UPLO,TRANSA, DIAG, M, N, calphap, cAp, LDA, cBp, LDB)
  !      type is (complex(8))
  !      !if(same_type_as(A,B).and.same_type_as(A,alpha)) then
  !      !    call glapack_ztrsm (SIDE, UPLO,TRANSA, DIAG, M, N, alpha, A, LDA, B, LDB)
  !      !end if
  !      call glapack_ztrsm (SIDE, UPLO,TRANSA, DIAG, M, N, zalphap, zAp, LDA, zBp, LDB)
  !    end select
  !
  !  end subroutine

  subroutine glapack_strsm (SIDE, UPLO,TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
    character   SIDE
    character   UPLO
    character   TRANSA
    character   DIAG
    integer     M
    integer     N
    real      ALPHA
    real, dimension(lda,*) ::   A
    integer     LDA
    real, dimension(ldb,*) ::   B
    integer     LDB
    integer :: iside,iuplo,tra,idiag,istat
#ifdef MAGMA
    !external magmablas_strsm
    !#   ifdef CUDAFOR
    integer :: h1
    call cublasXt_Create(h1)
    !#    endif
#endif

#ifdef MAGMA
    !#    ifdef CUDAFOR
    iside=glapack_getSide(SIDE); iuplo=glapack_getUplo(UPLO); tra=glapack_getTrans(TRANSA); idiag=glapack_getDiag(DIAG)
    istat=cublasXt_Strsm (h1, iside, iuplo,tra, idiag, M, N, ALPHA, A, LDA, B, LDB)
    !#    else
    !istat = cublasXt_strsm (h1, SIDE, UPLO,TRANSA, idiag, M, N, ALPHA, A, LDA, B, LDB)
    !#    endif
#else
    call strsm (SIDE, UPLO,TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
#endif
  end subroutine

  subroutine glapack_dtrsm (SIDE, UPLO,TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
    character   SIDE
    character   UPLO
    character   TRANSA
    character   DIAG
    integer     M
    integer     N
    real(8)      ALPHA
    real(8), dimension(lda,*) ::   A
    integer     LDA
    real(8), dimension(ldb,*)  ::  B
    integer     LDB
    integer :: iside,iuplo,tra,idiag,istat
#ifdef MAGMA
    !external cublasXt_dtrsm
    !#    ifdef CUDAFOR
    integer :: h1
    call cublasXt_Create(h1)
    !#    endif
#endif

#ifdef MAGMA
    !#    ifdef CUDAFOR
    iside=glapack_getSide(SIDE); iuplo=glapack_getUplo(UPLO); tra=glapack_getTrans(TRANSA); idiag=glapack_getDiag(DIAG)
    istat=cublasXt_Dtrsm (h1, iside, iuplo, tra, idiag, M, N, ALPHA, A, LDA, B, LDB)
    !#    else
    !istat = cublasXt_dtrsm (h1, SIDE, UPLO,TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
    !#    endif
#else
    call dtrsm (SIDE, UPLO,TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
#endif
  end subroutine

  subroutine glapack_ctrsm (SIDE, UPLO,TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
    character   SIDE
    character   UPLO
    character   TRANSA
    character   DIAG
    integer     M
    integer     N
    complex      ALPHA
    complex, dimension(lda,*) ::   A
    integer     LDA
    complex, dimension(ldb,*) ::   B
    integer     LDB
    integer :: iside,iuplo,tra,idiag,istat
#ifdef MAGMA
    !external cublasXt_ctrsm
    !#    ifdef CUDAFOR
    integer :: h1
    call cublasXt_Create(h1)
    !#    endif
#endif


#ifdef MAGMA
    !#    ifdef CUDAFOR
    iside=glapack_getSide(SIDE); iuplo=glapack_getUplo(UPLO); tra=glapack_getTrans(TRANSA); idiag=glapack_getDiag(DIAG)
    istat=cublasXt_Ctrsm (h1, iside, iuplo, tra, idiag, M, N, ALPHA, A, LDA, B, LDB)
    !#    else
    !istat = cublasXt_ctrsm (h1, SIDE, UPLO,TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
    !#    endif
#else
    call ctrsm (SIDE, UPLO,TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
#endif
  end subroutine

  subroutine glapack_ztrsm (SIDE, UPLO,TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
    character   SIDE
    character   UPLO
    character   TRANSA
    character   DIAG
    integer     M
    integer     N
    complex(8)      ALPHA
    complex(8), dimension(lda,*) ::   A
    integer     LDA
    complex(8), dimension(ldb,*)  ::  B
    integer     LDB
    integer :: iside,iuplo,tra,idiag,istat
#ifdef MAGMA
    !external cublasXt_Ztrsm
    !#ifdef CUDAFOR
    integer :: h1
    call cublasXt_Create(h1)
    !#endif
#endif

#ifdef MAGMA
    !#    ifdef CUDAFOR
    iside=glapack_getSide(SIDE); iuplo=glapack_getUplo(UPLO); tra=glapack_getTrans(TRANSA); idiag=glapack_getDiag(DIAG)
    istat=cublasXt_Ztrsm(h1, iside, iuplo, tra, idiag, M, N, ALPHA, A, LDA, B, LDB)
    !#    else
    !istat = cublasXt_ztrsm (h1, SIDE, UPLO,TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
    !#    endif
#else
    call ztrsm (SIDE, UPLO,TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
#endif
  end subroutine

  subroutine glapack_slarnv(idist,iseed,n,x)
    integer idist,iseed,n
    real x(n)
    !#ifdef MAGMA
    !    call magmaf_slarnv(idist,iseed,n,x)
    !#else
    call slarnv(idist,iseed,n,x)
  !#endif
  end subroutine

  subroutine glapack_dlarnv(idist,iseed,n,x)
    integer idist,iseed,n
    real(8) x(n)
    !#ifdef MAGMA
    !    call magmaf_dlarnv(idist,iseed,n,x)
    !#else
    call dlarnv(idist,iseed,n,x)
  !#endif
  end subroutine

  subroutine glapack_clarnv(idist,iseed,n,x)
    integer idist,iseed,n
    complex x(n)
    !#ifdef MAGMA
    !    call magmaf_clarnv(idist,iseed,n,x)
    !#else
    call clarnv(idist,iseed,n,x)
  !#endif
  end subroutine


  subroutine glapack_zlarnv(idist,iseed,n,x)
    integer idist,iseed,n
    complex(8) x(n)
    !#ifdef MAGMA
    !    call magmaf_zlarnv(idist,iseed,n,x)
    !#else
    call zlarnv(idist,iseed,n,x)
  !#endif
  end subroutine

  subroutine glapack_clacgv(n,x,incx)
    integer:: n
    complex:: x(:)
    integer incx
    !#ifdef MAGMA
    !    call magmaf_clacgv(n,x,incx)
    !    !call clacgv(n,x,incx)
    !#else
    call clacgv(n,x,incx)
  !#endif
  end subroutine

  subroutine glapack_zlacgv(n,x,incx)
    integer:: n
    complex(8):: x(:)
    integer incx
    !#ifdef MAGMA
    !    call magmaf_zlacgv(n,x,incx)
    !    !call zlacgv(n,x,incx)
    !#else
    call zlacgv(n,x,incx)
  !#endif
  end subroutine


  subroutine wt_info(str,info)
    character(*):: str
    integer:: info
    if(info/=0) then
      write(6,'(2a,a,i7)') 'non-zero info is obtained in ',str,' info=',info
      stop
    end if
  end subroutine
end module

           ! if(same_type_as(dA,dB).and.same_type_as(dA,dC).and. &
           !     & same_type_as(dA,alpha).and.same_type_as(dA,beta)) then
           !     call glapack_zgemm( transA, transB, m, n, k, alpha, dA, ldda, dB, lddb, beta,  &
           !         dC, lddc )
           ! end if


    !    subroutine glapack_evprm(h,ng,nb,lwork,lrwork,liwork)
    !        class(*) :: h(ng,nb)
    !        integer :: ng
    !        integer :: nb
    !        integer :: lwork
    !        integer :: lrwork
    !        integer :: liwork
    !        !integer :: lwmax
    !        if(rorc=='r') then
    !            call glapack_evprm_r(h,ng,lwork,lrwork,liwork)
    !        else if(rorc=='c') then
    !            call glapack_evprm_c(h,ng,lwork,lrwork,liwork)
    !        end if
    !    end subroutine
    !
    !    ! this estimate waork array size lwork, liwork
    !    subroutine glapack_evprm_s(h,ng,nb,lwork,lrwork,liwork)
    !        real:: h(ng,ng)
    !        integer :: ng
    !        integer :: nb
    !        integer :: lwork
    !        integer :: lrwork
    !        integer :: liwork
    !        !integer :: lwmax
    !        integer:: iwork(1)
    !        real:: rtmp(1)
    !        real:: workr(1)
    !        real:: rwork(1)
    !        real:: w(1)
    !        real :: h(1,1)
    !        integer :: info
    !        !integer :: lworkt
    !        lwork  = -1
    !        liwork = 0
    !        call glapack_gev( 'V', 'L', nb, h, ng, w, workr, lwork,rwork,lrwork,iwork, liwork, info )
    !        if(info/=0) call wt_info('magmaf_ssyevd',info)
    !        lwork = workr(1)
    !        lrwork=1
    !        liwork = iwork(1)
    !        !lworkt =1 + 6*nb + 2*nb**2
    !        !lwork = max(lwork,lworkt)    ! for strange lwork
    !        write(6,'(a,3i10)') 'lwork,lrwork,liwork=',lwork,lrwork,liwork
    !        if(info/=0) then
    !            write(6,'(a,3i10)') 'cannt estimate lapack_param, lwork,lrwork,liwork=',lwork,lrwork,liwork
    !            stop
    !        end if
    !    end subroutine
    !
    !    subroutine glapack_evprm_c(h,ng,nb,lwork,lrwork,liwork)
    !        complex:: h(ng,ng)
    !        integer :: ng
    !        integer :: nb
    !        integer :: lwork
    !        integer :: lrwork
    !        integer :: liwork
    !        !integer :: lwmax
    !        integer:: iwork(1)
    !        real:: rtmp(1)
    !        real:: rwork(1)
    !        complex:: workc(1)
    !        real:: w(1)
    !        complex :: h(1,1)
    !        integer :: info
    !        !integer lworkt
    !        lwork  = -1
    !        liwork = 0
    !        call glapack_ev( 'V', 'L', nb, h, ng, w, workc, lwork, rwork,lrwork,iwork, liwork, info )
    !        if(info/=0) call wt_info('magmaf_cheev',info)
    !        lwork = workc(1)
    !        lrwork = rwork(1)
    !        liwork = iwork(1)
    !        !lworkt =1 + 6*nb + 2*nb**2
    !        !lwork = max(lwork,lworkt)    ! for strange lwork
    !        write(6,'(a,3i10)') 'lwork,lrwork,liwork=',lwork,lrwork,liwork
    !        if(info/=0) then
    !            write(6,'(a,3i10)') 'cannt estimate lapack_param, lwork,lrwork,liwork=',lwork,lrwork,liwork
    !            stop
    !        end if
    !    end subroutine
