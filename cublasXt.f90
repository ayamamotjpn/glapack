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

module cublasXt
  use iso_c_binding

  interface cublasXt_Create
    subroutine cublasXt_Create(handle) bind(C,name='cublasXtCreate')
      use iso_c_binding
      integer handle
    end subroutine
  end interface

  interface cublasXt_sgemm
    integer function cublasXt_sgemm(handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) bind(C,name='cublasXtSgemm')
      use iso_c_binding
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
      use iso_c_binding
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
      use iso_c_binding
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
      use iso_c_binding
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
      use iso_c_binding
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
      use iso_c_binding
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
      use iso_c_binding
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
      use iso_c_binding
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
