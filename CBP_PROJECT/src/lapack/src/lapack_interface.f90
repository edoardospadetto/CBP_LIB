

#ifdef CBPCMPLX
#define CTYPE complex
#define CTYPE0 cmplx(0.0,0.0)
#define CTYPE1 cmplx(1.0,0.0)
#define MPI_CTYPE MPI_COMPLEX

#define pxelset pcelset
#define pxelget pcelget

#define xgesvd cgesvd
#define xgeqrf cgeqrf
#define xxxmqr cunmqr

#endif

#ifdef CBPDOUBLE
#define CTYPE double precision
#define CTYPE0 0.0
#define CTYPE1 1.0
#define MPI_CTYPE MPI_DOUBLE_PRECISION

#define pxelset pdelset
#define pxelget pdelget

#define xgeqrf dgeqrf
#define xxxmqr dormqr
#define xgesvd dgesvd
#endif

#ifdef CBPZMPLX
#define CTYPE0 dcmplx(0.0,0.0)
#define CTYPE1 dcmplx(1.0,0.0)
#define CTYPE double complex
#define MPI_CTYPE MPI_DOUBLE_COMPLEX

#define pxelset pzelset
#define pxelget pzelget
#define xgesvd zgesvd

#define xgeqrf zgeqrf
#define xxxmqr zunmqr

#endif

#ifdef CBPREAL
#define CTYPE0 0.0
#define CTYPE1 1.0
#define CTYPE real
#define MPI_CTYPE MPI_REAL

#define pxelset pselset
#define pxelget pselget
#define xgesvd sgesvd
#define xgeqrf sgeqrf
#define xxxmqr sormqr

#endif

!#include "misc_mod.f90"

 !call xgeqrf(dimA(1), dimA(2), Am , dimA(1), tau, query, lwork, info)
  !call xxxmqr('L', 'N',dimA(1), minval(dimA), minval(dimA), Am , dimA(1), tau, TensQ%elem,dimA(1), query, ii, info)


module lapack_interface
implicit none
contains
  subroutine QR_decompose_lapack(A,Q,R,info)

    integer :: lwork , info, row, col, ii
    CTYPE , dimension(:) , allocatable :: work
    CTYPE , dimension(:) , allocatable :: tau
    CTYPE , dimension(:,:) , allocatable :: A, A0, Q , R


    row = size(A,dim=1)
    col = size(A,dim=2)
    allocate (A0(row , col))
    allocate (tau(min(row , col)))
    allocate (work(1))
    A0 = A


    lwork = -1

    call xgeqrf(row , col, A , row, tau, work, lwork, info)
    if(info .ne.0) return
    lwork = ceiling(abs(work(1)))
    deallocate (work)
    allocate (work(lwork))
    call xgeqrf(row , col, A , row, tau, work, lwork, info)
    if(info .ne.0) return

    Q = CTYPE0
    do ii = 1 , min(row , col)
      Q(ii,ii) = CTYPE1
    end do

    lwork = -1
    !maybe the third parameter should me min(row,col)-1 but not sure.. 
    call xxxmqr('L', 'N', row , min(row , col), min(row , col), A , row , tau, Q , row , work , lwork , info)
    if(info .ne.0) return
    lwork = ceiling(abs(work(1)))
    deallocate (work)
    allocate (work(lwork))
    call xxxmqr('L', 'N', row , min(row , col), min(row , col), A , row , tau, Q , row , work , lwork , info)
    if(info .ne.0) return

    R=CTYPE0
    do ii = 1, col
      R(:min(ii,min(row,col)),ii) = A(:min(row,col),ii)
    end do
    A = A0
    print*, "check llapack qr" , abs(sum(matmul(conjg(transpose(Q)),Q))), sum(abs(A - matmul(Q,R)))/(row*col), row, col

    deallocate(A0)
    deallocate(tau)
    deallocate(work)

  end subroutine

  subroutine SVD_decompose_lapack(A,U,D,VT)
    CTYPE , dimension(:,:) , allocatable :: A, U, VT, A0
    double precision, dimension(:) :: D
    integer :: row, col,minrc, lwork, info
    double precision , dimension(:), allocatable :: rwork
    CTYPE , dimension(:), allocatable :: work


    row = size(A, dim=1)
    col = size(A, dim=2)
    minrc  = min(row,col)
    lwork = -1
    allocate(A0(row,col))
    A0 = A
    allocate(work(1))
    allocate(rwork(5*minrc))


    call zgesvd 	( 'S','S',row,col,A0,row,D,U,row,VT, minrc, work, lwork, rwork, info	)

    lwork = ceiling(abs(work(1)))
    deallocate(work)
    allocate(work(lwork))

    call zgesvd 	( 'S','S',row,col,A0,row,D,U,row,VT, minrc, work, lwork, rwork, info	)


    deallocate(rwork)
    deallocate(work)
    deallocate(A0)
  end subroutine
end module
