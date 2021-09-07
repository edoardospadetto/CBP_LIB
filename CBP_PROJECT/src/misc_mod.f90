#ifdef CBPZMPLX
#define CTYPE0 dcmplx(0.0,0.0)
#define CTYPE1 dcmplx(1.0,0.0)
#define CTYPE double complex
#define MPI_CTYPE MPI_DOUBLE_COMPLEX

#define p_gemm pzgemm
#define p_geqrf pzgeqrf
#define p_elset pzelset
#define p_elget pzelget
#endif

#ifdef CBPCMPLX
#define CTYPE complex
#define CTYPE0 cmplx(0.0,0.0)
#define CTYPE1 cmplx(1.0,0.0)
#define MPI_CTYPE MPI_COMPLEX

#define p_gemm pcgemm
#define p_geqrf pcgeqrf
#define p_elset pcelset
#define p_elget pcelget
#endif

#ifdef CBPDOUBLE
#define CTYPE double precision
#define CTYPE0 0.0
#define CTYPE1 1.0
#define MPI_CTYPE MPI_DOUBLE_PRECISION

#define p_gemm pdgemm
#define p_geqrf pdgeqrf
#define p_elset pdelset
#define p_elget pdelget
#endif

#ifdef CBPREAL
#define CTYPE0 0.0
#define CTYPE1 1.0
#define CTYPE real
#define MPI_CTYPE MPI_REAL

#define p_gemm psgemm
#define p_geqrf psgeqrf
#define p_elset pselset
#define p_elget pselget
#endif


module misc_mod
	implicit none


	interface c_random
			module procedure cm_random
			module procedure ca_random
	end interface

	contains

	function lcg(s)
        integer :: lcg
        integer*8, intent(in out) :: s
        if (s == 0) then
            s = 104729
        else
            s = mod(s, 4294967 )
        end if
        s = mod(s * 279470273, 4294291)
        lcg = int(mod(s, int(huge(0), 8)), kind(0))
    end function lcg


	subroutine init_random_seed()
        implicit none
        integer, allocatable :: seed(:)
        integer :: i, n, istat, dt(8), pid
        integer*8 :: t

        integer, parameter :: un=703

        call random_seed(size = n)
        allocate(seed(n))
        ! First try if the OS provides a random number generator
        open(unit=un, file="/dev/urandom", access="stream", &
            form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
            read(un) seed
            close(un)
        else
            ! The PID is
            ! useful in case one launches multiple instances of the same
            ! program in parallel.
            call system_clock(t)
            if (t == 0) then
                call date_and_time(values=dt)
                t = (dt(1) - 1970) * 365 * 24 * 60 * 60 * 1000 &
                + dt(2) * 31 * 24 * 60 * 60 * 1000 &
                + dt(3) * 24 * 60 * 60 * 1000 &
                + dt(5) * 60 * 60 * 1000 &
                + dt(6) * 60 * 1000 + dt(7) * 1000 &
                + dt(8)
            end if
            pid = 53
            t = ieor( t, int(pid, kind(t)) )
            do i = 1, n
                seed(i) = lcg(t)
            end do
        end if
        call random_seed(put=seed)
        !print*, "optimal seed = ", seed
    end subroutine init_random_seed


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine ca_random(A)
	implicit none
	CTYPE, dimension(:) :: A

#ifdef CBPCMPLX
        real , dimension (:), allocatable  ::  rai, rar
#endif
#ifdef CBPZMPLX
        double precision , dimension(:) , allocatable :: rar, rai
#endif
	integer, allocatable :: seed(:)
	integer :: seed_size
    	 call random_seed(size=seed_size)
	 allocate(seed(seed_size))

	 call init_random_seed()
#if defined(CBPCMPLX) |defined (CBPZMPLX)
	    allocate(rar(size(A)))
	    allocate(rai(size(A)))

	    call random_number(rar)
	    call random_number(rai)
#endif
#ifdef CBPCMPLX
	  A = cmplx(rar, rai)
#endif

#ifdef CBPZMPLX
	 A = dcmplx(rar, rai)
#endif

#if defined(CBPCMPLX) |defined (CBPZMPLX)
	   deallocate(rar)
	   deallocate(rai)
#endif


#if defined(CBPREAL) | defined(CBPDOUBLE)
    !Fill A, B with random numbers

    call random_number(A)

#endif

end subroutine
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine cm_random(A)
implicit none
CTYPE, dimension(:,:) :: A

#ifdef CBPCMPLX
			real , dimension (:,:), allocatable  ::  rai, rar
#endif
#ifdef CBPZMPLX
			double precision , dimension(:,:) , allocatable :: rar, rai
#endif
integer, allocatable :: seed(:)
integer :: seed_size
		 call random_seed(size=seed_size)
 allocate(seed(seed_size))

 call init_random_seed()
#if defined(CBPCMPLX) |defined (CBPZMPLX)
		allocate(rar(size(A, dim =1) , size(A, dim=2) ))
		allocate(rai(size(A, dim =1) , size(A, dim=2) ))

		call random_number(rar)
		call random_number(rai)
#endif
#ifdef CBPCMPLX
	A = cmplx(rar, rai)
#endif

#ifdef CBPZMPLX
 A = dcmplx(rar, rai)
#endif

#if defined(CBPCMPLX) |defined (CBPZMPLX)
	 deallocate(rar)
	 deallocate(rai)
#endif


#if defined(CBPREAL) | defined(CBPDOUBLE)
	!Fill A, B with random numbers

	call random_number(A)

#endif

end subroutine
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!Find first value equal to the requested in an array, if not found resturns -1
function findval(array, val) result(ii)
    integer , dimension(:) :: array
    integer :: val,ii
    logical :: good

    good = .false.
    do ii = 1, size(array)
        if(array(ii) .eq. val) then
        good = .true.
        return
        exit
        end if
    end do

    ii = -1

    return


    end function

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




end module misc_mod
