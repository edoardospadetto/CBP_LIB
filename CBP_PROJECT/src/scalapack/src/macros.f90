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


