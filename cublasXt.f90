module cublasXt
  use iso_c_binding

  interface cublasXt_Create
    subroutine cublasXt_Create(handle) bind(C,name='cublasXtCreate')
      integer handle
    end subroutine
  end interface

  interface cublasXt_sgemm
    integer function cublasXt_sgemm(handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) bind(C,name='cublasXtSgemm')
      integer handle
      integer transa, transb
      integer m, n, k
      real    alpha
      real    A(lda,*)
      integer lda
      real    B(ldb,*)
      integer ldb
      real    beta
      real    C(ldc,*)
      integer ldc
    end function
  end interface

  interface cublasXt_dgemm
    integer function cublasXt_dgemm(handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) bind(C,name='cublasXtDgemm')
      integer handle
      integer transa, transb
      integer m, n, k
      real(8) alpha
      real(8) A(lda,*)
      integer lda
      real(8) B(ldb,*)
      integer ldb
      real(8) beta
      real(8) C(ldc,*)
      integer ldc
    end function
  end interface

  interface cublasXt_cgemm
    integer function cublasXt_cgemm(handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) bind(C,name='cublasXtCgemm')
      integer handle
      integer transa, transb
      integer m, n, k
      complex alpha
      complex A(lda,*)
      integer lda
      complex B(ldb,*)
      integer ldb
      complex beta
      complex C(ldc,*)
      integer ldc
    end function
  end interface

  interface cublasXt_zgemm
    integer function cublasXt_zgemm(handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) bind(C,name='cublasXtZgemm')
      integer handle
      integer transa, transb
      integer m, n, k
      complex(8) alpha
      complex(8) A(lda,*)
      integer lda
      complex(8) B(ldb,*)
      integer ldb
      complex(8) beta
      complex(8) C(ldc,*)
      integer ldc
    end function
  end interface


  interface cublasXt_strsm
    integer function cublasXt_strsm(handle,side,uplo,trans,diag,m,n,alpha,A,lda,B,ldb) bind(C,name='cublasXtStrsm')
      integer handle
      integer side
      integer uplo
      integer trans
      integer diag
      integer m, n
      real    alpha
      real    A(lda,*)
      integer lda
      real    B(ldb,*)
      integer ldb
    end function
  end interface

  interface cublasXt_dtrsm
    integer function cublasXt_dtrsm(handle,side,uplo,trans,diag,m,n,alpha,A,lda,B,ldb) bind(C,name='cublasXtDtrsm')
      integer handle
      integer side
      integer uplo
      integer trans
      integer diag
      integer m, n
      real(8) alpha
      real(8) A(lda,*)
      integer lda
      real(8) B(ldb,*)
      integer ldb
    end function
  end interface

  interface cublasXt_ctrsm
    integer function cublasXt_ctrsm(handle,side,uplo,trans,diag,m,n,alpha,A,lda,B,ldb) bind(C,name='cublasXtCtrsm')
      integer handle
      integer side
      integer uplo
      integer trans
      integer diag
      integer m, n
      complex alpha
      complex A(lda,*)
      integer lda
      complex B(ldb,*)
      integer ldb
    end function
  end interface

  interface cublasXt_ztrsm
    integer function cublasXt_ztrsm(handle,side,uplo,trans,diag,m,n,alpha,A,lda,B,ldb) bind(C,name='cublasXtZtrsm')
      integer handle
      integer side
      integer uplo
      integer trans
      integer diag
      integer m, n
      complex(8) alpha
      complex(8) A(lda,*)
      integer lda
      complex(8) B(ldb,*)
      integer ldb
    end function
  end interface


end module
