#ifdef CBPREAL
#define CTYPE   REAL
#define CTYPE1 1.0
#define CTYPE0 0.0
#endif

#ifdef CBPCMPLX
#define CTYPE complex
#define CTYPE1 cmplx(1.0,0.0)
#define CTYPE0 cmplx(0.0,0.0)
#endif

#ifdef CBPZMPLX
#define CTYPE double complex
#define CTYPE1 dcmplx(1.0,0.0)
#define CTYPE0 dcmplx(0.0,0.0)
#endif

#ifdef CBPDOUBLE
#define CTYPE double precision
#define CTYPE1 1.d0
#define CTYPE0 0.d0
#endif



module cbpfor_parallel
        use iso_c_binding
        implicit none

        contains

        subroutine destroy_xt_context(h)
                type(c_ptr) :: h

                call destroy_xt_context_wrapped(h)

        end subroutine

        subroutine init_xt_context(ngpu, h)
               integer:: ngpu
               type(c_ptr) :: h
               call init_xt_context_wrapped(ngpu, h)
        end subroutine

        !ngpu is useless since h (cublasxt context contains the info)
        subroutine cbpfor_auto_gemm(A,B,C,add,ngpu,cpuratio,h, dimAi, dimBi, dimCi)
            CTYPE, dimension(:,:) :: A,B,C
            integer :: dimA(2), dimB(2), dimC(2)
            integer :: ngpu
            integer :: add
            integer, dimension(2) , optional :: dimAi, dimBi, dimCi
            real:: cpuratio
            type(c_ptr):: h

            if (.not. present(dimAi) ) then
            dimA = shape(A)
            else
            dimA = dimAi
            end if

            if (.not. present(dimBi) ) then
            dimB = shape(B)
            else
            dimB = dimBi
            end if

            if (.not. present(dimCi) ) then
            dimC = shape(C)
            else
            dimC = dimCi
            end if


            if (dimA(2) .ne. dimB(1)) then
                print*, "STOP, error size ", __FILE__, __LINE__
                STOP
            else if (dimA(1) .ne. dimC(1)) then
            	print*, "STOP, error size ", __FILE__, __LINE__
                STOP
            else if (dimB(2) .ne. dimC(2)) then
            	print*, "STOP, error size ", __FILE__, __LINE__
                STOP
            end if

            call cbpfor_auto_gemm_wrapped(A,B,C,dimA,dimB,dimC,add,ngpu, cpuratio, h)

        end subroutine

  subroutine cbpfor_tensorout(v,C,add,ngpu,cpuratio,h, dimvi, dimCi)
	CTYPE, dimension(:) :: v
	CTYPE, dimension(:,:) :: C

	integer :: ngpu, ii
	integer :: add, dimv
	integer, dimension(2) :: dimC
	integer, dimension(2) , optional :: dimCi
	integer, optional :: dimvi
	real:: cpuratio
	type(c_ptr):: h

	if (.not. present(dimvi) ) then
	dimv = size(v)
	else
	dimv = dimvi
	end if


	if (.not. present(dimCi) ) then
	dimC = shape(C)
	else
	dimC = dimCi
	end if


	if (dimv.ne. dimC(1)) then
	print*, "STOP, error size", __FILE__, __LINE__
	STOP
	else if (dimv .ne. dimC(2)) then
	print*, "STOP, error size", __FILE__, __LINE__
	STOP
	endif

#if defined(CBPCMPLX) | defined(CBPZMPLX)
	    !do ii = 1, dimv
	    !	v(ii) = cmplx(ii,1.0)
	    !end do
	    !C = CTYPE0

            call cbpfor_auto_gemm_wrapped(v,conjg(v),C(:dimC(1),:dimC(2)),[dimv, 1],[1,dimv],dimC,add,ngpu, cpuratio, h)

#endif
#if defined(CBPREAL) | defined(CBPDOUBLE)
            call cbpfor_auto_gemm_wrapped(v,v,C,[dimv, 1],[1,dimv],dimC,add,ngpu, cpuratio, h)
#endif
        end subroutine

        subroutine cbpfor_outer(u,v,C,add,ngpu,cpuratio,h, dimui, dimvi, dimCi)
        	CTYPE, dimension(:) :: v, u
		CTYPE, dimension(:,:) :: C

		integer :: ngpu, ii
		integer :: add, dimv , dimu
		integer, dimension(2) :: dimC
		integer, dimension(2) , optional :: dimCi
		integer, optional :: dimvi, dimui
		real:: cpuratio
		type(c_ptr):: h


		if (.not. present(dimvi) ) then
		dimv = size(v)
		else
		dimv = dimvi
		end if
		if (.not. present(dimui) ) then
		dimu = size(u)
		else
		dimu = dimui
		end if


		if (.not. present(dimCi) ) then
		dimC = shape(C)
		else
		dimC = dimCi
		end if


		if (dimu.ne. dimC(1)) then
		print*, "STOP, error size u ne C"
		STOP
		else if (dimv .ne. dimC(2)) then
		print*, "STOP, error size v ne C"
		STOP
		endif

        	call cbpfor_auto_gemm_wrapped(u,v,C(:dimC(1),:dimC(2)),[dimu, 1],[1,dimv],dimC,add,ngpu,cpuratio, h)

        end subroutine
        !here ngpu is usefull
        !subroutine cbpfor_auto_svd(A,S,V,D, ngpu)

        !	CTYPE , dimension(:,:) :: A, S, V
        !	CTYPE, dimension(:) :: D
        !	integer :: ngpu

        	!split A in blocks

        !	call cbpfor_auto_svd_wrapped (A, S, V, D, size(A,dim=1) , size(A, dim=2))

        	!each block to different processes

        	!recover full svd.







        !end subroutine

end module
