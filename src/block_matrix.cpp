/***************************************************************************
*     CLASS db_matrix                                                      *
*                                                                          *
*     MINATEC-INPG, 3, parvis Luis Neel                                    *
*     38016, Grenoble CEDEX, France                                        *
*                                                                          *
*     bucci@minatec.grenoble-inp.fr                                        *
*                                                                          *
****************************************************************************
*     Provides a class of methods for easily manipulating double complex   *
*     matrices easily used with Fortran methods.                           *
*     This class makes use of LAPACK. This library should be installed     *
*     on the system, along with a BLAS library. I strongly recommend using *
*     a version of BLAS optimized for the target system. In most cases,    *
*     ATLAS should achieve good results in comparison with the reference   *
*     BLAS included in the standard LAPACK distribution.                   *
****************************************************************************/

/* This file is part of cRCWA.

    cRCWA is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.
    
    cRCWA is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
    FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
    details.
    
    You should have received a copy of the GNU General Public License along
    with cRCWA. If not, see <https://www.gnu.org/licenses/>. 

    Davide Bucci, 2008-2025
    Jérôme Michallon, 2012-2014
*/


/* Version history

    1.0    Dec, 13 2008  First version with the correct Toeplitz calculation
    1.0.1  Jan,    2009  A few memory leaks corrected
    1.0.2  Feb, 17 2009  Eigenvectors calculation in eig function
    1.0.3  Feb, 26 2009  corrected a bug with eig(..., NULL)
                         eig() does not allocate memory for eigenvectors if they
                         are not calculated
    1.0.4  Apr, 1  2009  operator * between matrices does not return a const
    1.0.5  Apr, 27 2009  throw an exception when unable to allocate memory
                         corrected a bug in matrix sizes when calculating fft
    1.0.6  May, 27 2009  corrected a few compilation issues on newest gcc's
    1.0.7  Aug, 25 2009  corrected a possible issue with kill.
    1.0.8  Nov, 15 2009  make use of ZGEEVX instead of ZGEEV.
    1.0.9  Apr, 12 2010  corrected a few bugs.
    1.0.10 Mar, 17 2011  explored a few possibilities about circulant matrices
    1.0.11 Mar, 28 2011  a little more tolerant while copying empty matrices
    1.0.12 Nov, 10 2011  bug corrected in the zero padding routine.

    ... then SVN has been adopted!

*/


// Comment to avoid using FFTW. The fft2 and ifft2 routines will not work
// in that case
#define USE_FFTW
//#define USE_KISSFFT

// Define this as true to use the expert version of the LAPACK eigenvalue
// driver routine for eigenvalue calculation
#define ZGEEVX true

// Define this to use BLAS for basic operations such as copying matrices or
// performing element-wise sums and differences.
#define BLASBASIC true

// You may have a description of each memory allocation and deallocation
// done by this class by uncommenting the following line. This tends to become
// extremely verbose, but it can be useful for finding memory leaks
// #define MEMORY_TRACK

#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include <fcntl.h>
#include <semaphore.h>
#include <errno.h>

#include "block_matrix.h"
#include "fortran_types.h"

// Unnamed semaphores are the simplest solutions on Linux, but they are not
// available on MacOSX, so we must provide both named and unnamed semaphores.
// TODO: check if this works as well in Linux
// #define USE_UNNAMED

// During the configuration operation, the configure script will try to
// determine if the BLAS and LAPACK libraries have been compiled with a
// Fortran compiler which adds an underscore before of after the name of
// the function to be called.

#if defined (BLAS_ADD_FINAL_UNDERSCORE)
    #define zgetri zgetri_
    #define zcopy zcopy_
    #define zaxpy zaxpy_
    #define zscal zscal_
#elif defined (BLAS_ADD_LEADING_UNDERSCORE)
    #define zgetri _zgetri
    #define zcopy _zcopy
    #define zaxpy _zaxpy
    #define zscal _zscal
    #warning BLAS_ADD_LEADING_UNDERSCORE
#elif defined (BLAS_ADD_BOTH_UNDERSCORES)
    #define zgetri _zgetri_
    #define zcopy _zcopy_
    #define zaxpy _zaxpy_
    #define zscal _zscal_
    #warning BLAS_ADD_BOTH_UNDERSCORES
#endif


#if defined (LAPACK_ADD_FINAL_UNDERSCORE)
    #define zgeev zgeev_
    #define zgeevx zgeevx_
    #define zggevx zggevx_
    #define zgemm zgemm_
    #define zgetrf zgetrf_
#elif defined (LAPACK_ADD_LEADING_UNDERSCORE)
    #define zgeev _zgeev
    #define zgeevx _zgeevx
    #define zggevx _zggevx
    #define zgemm _zgemm
    #define zgetrf _zgetrf
#elif defined (LAPACK_ADD_BOTH_UNDERSCORES)
    #define zgeev _zgeev_
    #define zgeevx _zgeevx_
    #define zggevx _zggevx_
    #define zgemm _zgemm_
    #define zgetrf _zgetrf_
#endif

int total_elements;
int max_elements;

extern "C" {
    // LAPACK and BLAS functions used in this class

    // LAPACK functions explicitely called
    int zgeev(char *jobvl, const char *jobvr, integer *n,
        doublecomplex *a, integer *lda, doublecomplex *w, doublecomplex *vl,
        integer *ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work,
        integer *lwork, doublereal *rwork, integer *info);

    void zgeevx(const char *balanc, const char *jobvl, const char *jobvr,
        const char *sense,
        integer *n,
        doublecomplex *a, integer *lda, doublecomplex *w, doublecomplex *vl,
        integer *ldvl, doublecomplex *vr, integer *ldvr, integer *ilo,
        integer *ihi,
        double *scale,  double *abnrm, double *rconde, double *rcondv,
        doublecomplex *work,
        integer *lwork, doublereal *rwork, integer *info);

    void zggevx(const char *balanc, const char *jobvl, const char *jobvr,
        const char *sense, integer *n,
        doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
        doublecomplex *alpha, doublecomplex *beta,
        doublecomplex *vl,
        integer *ldvl, doublecomplex *vr, integer *ldvr, integer *ilo,
        integer *ihi,
        double *lscale, double *rscale, double *abnrm, double *bbnrm,
        double *rconde, double *rcondv,
        doublecomplex *work,
        integer *lwork, doublereal *rwork,
        integer *iwork, bool *bwork,
        integer *info);

    int zgemm(char *TRANSA, char *TRANSB, integer *M, integer *N, integer *K,
        doublecomplex *ALPHA, doublecomplex *A, integer *LDA, doublecomplex *B,
        integer *LDB, doublecomplex *BETA, doublecomplex *C,integer *LDC);

    int zgetrf(integer *m, integer *n, doublecomplex *a,
        integer *lda, integer *ipiv, integer *info);

    // BLAS functions explicitely called
    int zgetri(integer *n, doublecomplex *a, integer *lda,
        integer *ipiv, doublecomplex *work, integer *lwork, integer *info);

    int zcopy(integer *n, doublecomplex *zx, integer *incx,
        doublecomplex *zy, integer *incy);

    int zaxpy(integer *n, doublecomplex *za, doublecomplex *zx, integer *incx,
        doublecomplex *zy, integer *incy);

    int zscal(integer *n, doublecomplex *za, doublecomplex *zx, integer *incx);
}

// For using with FFTW ver. 3
#ifdef USE_FFTW
    #include <fftw3.h>
#endif

// For using with KISSFFT
#ifdef USE_KISSFFT
    #include <kiss_fftnd.h>
#endif

using namespace std;

#include <semaphore.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "compat.h"

sem_t *mutex_fftw; // Semaphore for fftw
sem_t unnamed_mutex_fftw; // Unnamed semaphore (if applicable).
bool mutex_fftw_created=false;
#define SEM_NAME_SIZE 255
char semaphore_name[SEM_NAME_SIZE+1];


void init_semaphore_FFTW(void)
{
    #ifdef USE_UNNAMED
        // Code for the unnamed semaphore: the best solution for Linux, but
        // which does not work on MacOSX :-(
        sem_init(&unnamed_mutex_fftw, 0, 1);
        mutex_fftw=&unnamed_mutex_fftw;
        mutex_fftw_created=true;
    #else
        // Code for a named semaphore, which seems unnecessary here and
        // error-prone
        int p=getpid();
        snprintf(semaphore_name,SEM_NAME_SIZE,"/semaphore_fftw_%d", p);
        mutex_fftw=sem_open(semaphore_name, O_CREAT, 0644, 1);
        if(mutex_fftw==SEM_FAILED) {
            cerr<<"Could not set up a semaphore (FFTW3). "
                "DO NOT USE MULTIPLE THREADS!";
            cerr<<endl<<"errno="<<errno<<endl;
            cout<<"FFTW3 semaphore not set."<<endl;
            exit(1);
            mutex_fftw_created=false;
        } else {
            cout<<"Init semaphore OK."<<endl;
            mutex_fftw_created=true;
        }
    #endif
}

void delete_semaphore_FFTW(void)
{
    #ifdef USE_UNNAMED
        sem_destroy(mutex_fftw);
    #else
        if(mutex_fftw_created) {
            sem_post(mutex_fftw);
            // close semaphore
            sem_close(mutex_fftw);
            sem_unlink(semaphore_name);
            mutex_fftw_created=false;
            // This means that a call to sem_post or sem_wait would show that 
            // the semaphore has been closed.
            mutex_fftw=SEM_FAILED;
        }
    #endif
}
// Standard constructor: fabricate a null matrix in which the A pointer is not
// defined.
db_matrix::db_matrix(void)
{
    ncol = 0;
    nrow = 0;
    A=NULL;
}

// Standard constructor: fabricate a matrix of the given size
db_matrix::db_matrix(int row, int col)
{
    nrow=row;
    ncol=col;

    A =  new complex<double>[nrow*ncol];
    if (A==NULL)
        throw db_matrix_error("Error: could not allocate memory!");

    #ifdef MEMORY_TRACK
        //cout<<"File "<<__FILE__<<" at line "<<__LINE__<<"\n";
        //cout<<"Matrix allocated of size "<<nrow<<"x"<<ncol<<"\n";
        total_elements+=nrow*ncol;
        if(total_elements>max_elements) max_elements=total_elements;
        //cout << "Total number of elements in memory: "<<total_elements<<"\n";
    #endif

    int nelem=nrow*ncol;
    for(int i=0; i<nelem; ++i)
        A[i]=complex<double>(0,0);
}

// Copy constructor: element by element copy of the given matrix
db_matrix::db_matrix(const db_matrix& d)
{
    nrow=d.nrow;
    ncol=d.ncol;

    if(d.A == NULL) {
        A=NULL;
        return;
    }

    A =  new complex<double>[nrow*ncol];

    if (A==NULL)
        throw db_matrix_error("Error: could not allocate memory!");

    #ifdef MEMORY_TRACK
        //cout<<"File "<<__FILE__<<" at line "<<__LINE__<<"\n";
        //cout<<"Matrix allocated of size "<<nrow<<"x"<<ncol<<"\n";
        total_elements+=nrow*ncol;
        if(total_elements>max_elements) max_elements=total_elements;
        //cout << "Total number of elements in memory: "<<total_elements<<"\n";
    #endif
    #ifdef BLASBASIC
        integer nelem=nrow*ncol;
        integer incx = 1;
        integer incy = 1;
        // A memory hack based to the fact that we actually store
        // information in a Fortran compatible way
        doublecomplex *X = (doublecomplex *) d.A;
        doublecomplex *Z = (doublecomplex *) A;

        zcopy(&nelem, X, &incx, Z, &incy);
    #else
        for(int i=0;i<nrow;++i)
            for(int j=0;j<ncol;++j)
                A[i+j*nrow]=d.A[i+j*nrow];
    #endif
}

// Destructor
db_matrix::~db_matrix()
{
    if (A!=NULL) {
        delete[] A;
        A=NULL;
        #ifdef MEMORY_TRACK
          total_elements-=nrow*ncol;
        #endif
    }
    #ifdef MEMORY_TRACK
      else {
      }
    #endif
}

void db_matrix::kill(void)
{
    if(A==NULL) {
        nrow = 0;
        ncol = 0;
        return;
    }
    delete[] A;
    A=NULL;

    #ifdef MEMORY_TRACK
        total_elements-=nrow*ncol;
    #endif
    // Prevents using the matrix when it has been deallocated (1.0.7)
    nrow = 0;
    ncol = 0;
}

// Element by element copy of a given matrix
db_matrix& db_matrix::operator=(const db_matrix &d)
{
    // Same object?
    if (this == &d)
      return *this;
    #ifdef MEMORY_TRACK
        //cout<<"File "<<__FILE__<<" at line "<<__LINE__<<"\n";
        //cout << "Copying matrix\n";
    #endif

    if(A!=NULL) {
        #ifdef MEMORY_TRACK
            //cout<<"File "<<__FILE__<<" at line "<<__LINE__<<"\n";
            //cout<<"Matrix deallocated (copy over) of size "<<nrow<<"x"<<
             //   ncol<<"\n";
            total_elements-=nrow*ncol;
            //cout << "Total number of elements in memory: "<<total_elements
            //    <<"\n";
        #endif
        delete[] A;
        A=NULL;
    }
    nrow=d.nrow;
    ncol=d.ncol;

    //assert(d.A != NULL);

    if (d.A != NULL) {
        A =  new complex<double>[nrow*ncol];
        if (A==NULL)
            throw db_matrix_error("Error: could not allocate memory!");

        #ifdef MEMORY_TRACK
            //cout<<"File "<<__FILE__<<" at line "<<__LINE__<<"\n";
            //cout<<"Matrix allocated of size "<<nrow<<"x"<<ncol<<"\n";
            total_elements+=nrow*ncol;
            if(total_elements>max_elements) max_elements=total_elements;
           //cout << "Total number of elements in memory: "<<
            //   total_elements<<"\n";
        #endif

        #ifdef BLASBASIC
            integer nelem=nrow*ncol;
            integer incx = 1;
            integer incy = 1;
            // A memory hack based to the fact that we actually store
            // information in a Fortran compatible way
            doublecomplex *X = (doublecomplex *) d.A;
            doublecomplex *Z = (doublecomplex *) A;

            zcopy(&nelem, X, &incx, Z, &incy);
        #else
            for(int i=0;i<nrow;++i)
                for(int j=0;j<ncol;++j)
                    A[i+j*nrow]=d.A[i+j*nrow];
        #endif
    } else
        A=NULL;

    return *this;
}

// Multiplication by a complex coefficient
db_matrix& db_matrix::operator*=(const complex<double> &c)
{
    assert(A!=NULL);

    #ifdef BLASBASIC
        integer nelem=nrow*ncol;
        integer incx = 1;
        // A memory hack based to the fact that we actually store
        // information in a Fortran compatible way
        doublecomplex *X = (doublecomplex *) A;
        doublecomplex *za = (doublecomplex *) &c;

        zscal(&nelem, za, X, &incx);
    #else
        int nelem=nrow*ncol;
        for(int i=0; i<nelem; ++i)
            A[i]*=c;
    #endif

    return *this;
}

// Unary minus
db_matrix db_matrix::operator-()
{
    assert(A!=NULL);

    #ifdef BLASBASIC
        db_matrix R= *this;
        integer nelem=nrow*ncol;
        integer incx = 1;
        // A memory hack based to the fact that we actually store
        // information in a Fortran compatible way
        complex<double> sc = complex<double>(-1.0,0);

        doublecomplex *X = (doublecomplex *) R.A;
        doublecomplex *za = (doublecomplex *) &sc;

        zscal(&nelem, za, X, &incx);
    #else
        db_matrix R(nrow, ncol);
        int nelem = nrow*ncol;
        for(int i=0; i<nelem; ++i)
            R.A[i]=-A[i];
    #endif

    return R;
}

// Divide by a complex coefficient
db_matrix& db_matrix::operator/=(const complex<double> &c)
{
    assert(A!=NULL);

    return operator*=(1.0/c);
}

// Sum of two  matrices
db_matrix &db_matrix::operator+=(const db_matrix &d)
{
    assert(d.ncol == ncol);
    assert(d.nrow == nrow);
    assert(A!=NULL);
    assert(d.A!=NULL);

    #ifdef BLASBASIC
        integer nelem=nrow*ncol;
        integer incx = 1;
        integer incy = 1;
        // A memory hack based to the fact that we actually store
        // information in a Fortran compatible way
        doublecomplex *X = (doublecomplex *) d.A;
        doublecomplex *Y = (doublecomplex *) A;
        doublecomplex c = {1.0, 0.0};

        // SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
        // ZAXPY constant times a vector plus a vector
        // y = a*x + y

        zaxpy(&nelem, &c, X, &incx, Y, &incy);
    #else
        int nelem=nrow*ncol;
        for(int i=0; i<nelem; ++i)
            A[i]+=d.A[i];
    #endif


    return *this;
}

// Substraction of two matrices
db_matrix &db_matrix::operator-=(const db_matrix &d)
{
    assert(d.ncol == ncol);
    assert(d.nrow == nrow);
    assert(A!=NULL);
    assert(d.A!=NULL);

    #ifdef BLASBASIC
        integer nelem=nrow*ncol;
        integer incx = 1;
        integer incy = 1;
        // A memory hack based to the fact that we actually store
        // information in a Fortran compatible way
        doublecomplex *X = (doublecomplex *) d.A;
        doublecomplex *Y = (doublecomplex *) A;
        doublecomplex c = {-1.0, 0.0};

        // SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
        // ZAXPY constant times a vector plus a vector
        // y = a*x + y

        zaxpy(&nelem, &c, X, &incx, Y, &incy);
    #else
        int nelem=nrow*ncol;
        for(int i=0; i<nelem; ++i)
            A[i]-=d.A[i];
    #endif

    return *this;
}

// Row by column multiplication of two matrices
db_matrix &db_matrix::operator*=(const db_matrix &d)
{
    assert(d.nrow == ncol);
    assert(A!=NULL);
    assert(d.A!=NULL);

    // Use the level 3 BLAS ZGEMM routine to multiply matrices
    char TransA = 'N';
    char TransB = 'N';
    doublecomplex alpha = {1,0};
    doublecomplex beta = {0,0};

    doublecomplex *Ar;
    doublecomplex *Br;
    doublecomplex *Cr;
    integer M = nrow;
    integer N = d.ncol;
    integer K = ncol;
    integer LDA = M;
    integer LDB = K;
    integer LDC = M;

    // A memory hack based to the fact that we actually store information
    // in a Fortran compatible way
    Ar = (doublecomplex*) A;
    Br = (doublecomplex*) d.A;
    Cr = new doublecomplex[nrow*d.ncol];

    // Probably this is not needed. However, valgrind complains that there
    // is some "Conditional jump or move depends on uninitialised value(s)"
    // if this code is not executed. I trust it much more than my naive
    // intuition!
    int nelem=nrow*d.ncol;
    for(int i=0; i< nelem; ++i)
        Cr[i]=beta; 

    #ifdef MEMORY_TRACK
            total_elements+=nrow*d.ncol;
            if(total_elements>max_elements) max_elements=total_elements;
    #endif

    zgemm(&TransA, &TransB, &M, &N, &K, &alpha, Ar, &LDA,
        Br, &LDB, &beta,  Cr, &LDC);
    delete[] A;
    A=NULL;
    #ifdef MEMORY_TRACK
        total_elements-=nrow*ncol;
    #endif
    A=(complex<double>*)Cr;
    ncol = d.ncol;

    return *this;
}

// Row by column multiplication of two matrices. Scale the first matrix:
// X *= alpha*d
// where d is a second matrix and alpha is a constant
db_matrix &db_matrix::scalemult(const complex<double> alpha, const db_matrix &d)
{
    assert(d.nrow == ncol);
    assert(A!=NULL);
    assert(d.A!=NULL);

    // Use the level 3 BLAS ZGEMM routine to multiply matrices
    char TransA = 'N';
    char TransB = 'N';
    doublecomplex alpha_f = {alpha.real(),alpha.imag()};
    doublecomplex beta = {0,0};

    doublecomplex *Ar;
    doublecomplex *Br;
    doublecomplex *Cr;
    integer M = nrow;
    integer N = d.ncol;
    integer K = ncol;
    integer LDA = M;
    integer LDB = K;
    integer LDC = M;

    // A memory hack based to the fact that we actually store information
    // in a Fortran compatible way
    Ar = (doublecomplex*) A;
    Br = (doublecomplex*) d.A;
    Cr = new doublecomplex[nrow*d.ncol];

    // Probably this is not needed. However, valgrind complains that there
    // is some "Conditional jump or move depends on uninitialised value(s)"
    // if this code is not executed. I trust it much more than my naive
    // intuition!
    for(int i=0; i< nrow*d.ncol; ++i)
        Cr[i]=beta;

    #ifdef MEMORY_TRACK
            total_elements+=nrow*d.ncol;
            if(total_elements>max_elements) max_elements=total_elements;
    #endif

    zgemm(&TransA, &TransB, &M, &N, &K, &alpha_f, Ar, &LDA,
        Br, &LDB, &beta,  Cr, &LDC);
    delete[] A;
    A=NULL;
    #ifdef MEMORY_TRACK
        total_elements-=nrow*ncol;
    #endif
    A=(complex<double>*)Cr;
    ncol = d.ncol;

    return *this;
}


/** Row by column multiplication of two matrices and sum:
  C = A*B+C
  C matrix is loaded with the contents of the operation
*/
db_matrix &multsum(const db_matrix &A, const db_matrix &B, db_matrix &C)
{
    return multsumscale(complex<double>(1.0,0), A, B,
        complex<double>(1.0,0), C);
}

/** Row by column multiplication of two matrices and sum:
  C = aleph*AA*B+beth*C
  C matrix is loaded with the contents of the operation
*/
db_matrix &multsumscale(const complex<double> aleph, const db_matrix &AA,
    const db_matrix &B, const complex<double> beth, db_matrix &C)
{

    assert(B.nrow == AA.ncol);
    assert(AA.A!=NULL);
    assert(B.A!=NULL);
    assert(C.A!=NULL);

    // Use the level 3 BLAS ZGEMM routine to multiply matrices
    char TransA = 'N';
    char TransB = 'N';
    doublecomplex alpha = {aleph.real(), aleph.imag()};
    doublecomplex beta = {beth.real(), beth.imag()};

    doublecomplex *Ar;
    doublecomplex *Br;
    doublecomplex *Cr;
    integer M = AA.nrow;
    integer N = B.ncol;
    integer K = AA.ncol;
    integer LDA = M;
    integer LDB = K;
    integer LDC = M;

    // A memory hack based to the fact that we actually store information
    // in a Fortran compatible way
    Ar = (doublecomplex*) AA.A;
    Br = (doublecomplex*) B.A;
    Cr = (doublecomplex*) C.A;

    zgemm(&TransA, &TransB, &M, &N, &K, &alpha, Ar, &LDA,
        Br, &LDB, &beta,  Cr, &LDC);

    return C;
}

// Create a test matrix of the given size
db_matrix db_matrix::createTestMatrix(const int nrow, const int ncol)
{
    int i, j;
    int a, b;

    db_matrix d(nrow, ncol);


    for(i=0; i<nrow; ++i) {
        for(j=0; j<ncol; ++j) {

            if(j>(ncol/2))
                a=j-ncol;
            else
                a=j;

            if(i>(nrow/2))
                b=i-nrow;
            else
                b=i;

            d(i,j) = complex<double>(double(a), double(b));
        }
    }

    return d;
}

// Create a matrix full of ones, of the given size
db_matrix db_matrix::createOnesMatrix(const int nrow, const int ncol)
{
    int i, j;
    int a, b;

    db_matrix d(nrow, ncol);
    int nelem=nrow*ncol;

    for(i=0; i<nrow; ++i) {
        for(j=0; j<ncol; ++j) {
            d(i,j) = complex<double>(1.0, 0.0);
        }
    }

    return d;
}

// Create an unit matrix
db_matrix db_matrix::createUnitMatrix(const int nrow, const int ncol)
{
    db_matrix d(nrow, ncol);

    for(int i=0; i<nrow && i<ncol; ++i) {
        d(i,i) = complex<double>(1.0, 0);
    }
    return d;
}

// Calculate the Fourier transform of a rectangle and store it in the usual
// FFT way.
db_matrix db_matrix::rectangleFT(double sx, // Total size in the x direction
    double sy,              // Total size in the y direction
    double periodx,         // x-period (size) of the window
    double periody,         // y-period (size) of the window
    int Nx,                 // number of x Fourier coefficients
    int Ny,                 // number of y Fourier coefficients
    double centerx,         // center of the rectangle in the window: x
    double centery)         // center of the rectangle in the window: y
{
    double nux=2.0*M_PI/periodx;
    double nuy=2.0*M_PI/periody;

    int ii, jj;
    int i, j;

    db_matrix d(Ny, Nx);

    for(ii=0;ii<Ny;++ii) {
        for(jj=0;jj<Nx;++jj) {

            if ((ii+1)<Ny/2.0+1)
               i=ii;
            else
               i=-(Ny-ii);


            if ((jj+1)<Nx/2.0+1)
               j=jj;
            else
               j=-(Nx-jj);

            d(ii,jj)=(sx)/periodx*(sy)/periody*sinc(sx*double(j)*nux/2.0)*
                sinc(sy*double(i)*nuy/2.0);
            d(ii,jj)=d(ii,jj)*exp(complex<double>(0,1)*(nux*double(j)*
                centerx+nuy*double(i)*centery));
            // For test purpouses only
            //d(ii,jj)=complex<double>(i,j);
       }
    }
    return d;
}

/** Truncation of a FFT matrix
*/
db_matrix db_matrix::windowing(void)
{
    int ii, jj;
    int i, j;

    int Nx = ncol;
    int Ny = nrow;

    for(ii=0;ii<Ny;++ii) {
        for(jj=0;jj<Nx;++jj) {

            if ((ii+1)<Ny/2.0+1)
               i=ii;
            else
               i=-(Ny-ii);


            if ((jj+1)<Nx/2.0+1)
               j=jj;
            else
               j=-(Nx-jj);

            if(j>Nx/4 || j<-Nx/4 || i>Ny/4 || i<-Ny/4)
                this->operator()(ii,jj)=0;
       }
    }
    return *this;
}

// Bloc-merge two matrices on the same row. They must have the same number of
// rows
db_matrix db_matrix::mergeMatrixOnRow(const db_matrix &M1,
    const db_matrix &M2)
{
    int i,j;
    int ncol=M1.ncol + M2.ncol;

    assert(M1.A!=NULL);
    assert(M2.A!=NULL);
    assert(M1.nrow==M2.nrow);

    db_matrix d(M1.nrow, ncol);

    // Loop on all points

    for(i=0; i<M1.nrow; ++i) {
        for(j=0; j<M1.ncol; ++j) {
            d(i,j) = M1(i,j);
        }
        for(j=0; j<M2.ncol; ++j) {
            d(i,j+M1.ncol) = M2(i,j);
        }
    }

    return d;
}

// Bloc-merge two matrices on the same column. They must have the same number
// of columns
db_matrix db_matrix::mergeMatrixOnColumn(const db_matrix &M1,
    const db_matrix &M2)
{
    int i,j;
    int nrow=M1.nrow + M2.nrow;

    assert(M1.A!=NULL);
    assert(M2.A!=NULL);

    assert(M1.ncol==M2.ncol);

    db_matrix d(nrow, M1.ncol);

    for(i=0; i<M1.nrow; ++i) {
        for(j=0; j<M1.ncol; ++j) {
            d(i,j) = M1(i,j);
        }
    }
    for(i=0; i<M2.nrow; ++i) {
        for(j=0; j<M2.ncol; ++j) {
            d(i+M1.nrow,j) = M2(i,j);
        }
    }
    return d;
}

// Block-merge four matrices:
//    M1  M2
//    M3  M4

db_matrix db_matrix::mergeMatrixQuad(const db_matrix &M1,
    const db_matrix &M2,const db_matrix &M3,
    const db_matrix &M4)
{
    int i,j;
    int ncol=M1.ncol + M2.ncol;

    assert(M1.A!=NULL);
    assert(M2.A!=NULL);
    assert(M3.A!=NULL);
    assert(M4.A!=NULL);
    assert(M1.nrow==M2.nrow);
    assert(M3.nrow==M4.nrow);
    assert(M1.ncol+M2.ncol==M3.ncol+M4.ncol);

    db_matrix d(M1.nrow+M3.nrow, ncol);

    // Loop on all points

    for(i=0; i<M1.nrow; ++i) {
        for(j=0; j<M1.ncol; ++j) {
            d(i,j) = M1(i,j);
        }
        for(j=0; j<M2.ncol; ++j) {
            d(i,j+M1.ncol) = M2(i,j);
        }
    }
    for(i=0; i<M3.nrow; ++i) {
        for(j=0; j<M3.ncol; ++j) {
            d(i+M1.nrow,j) = M3(i,j);
        }
    }
    for(i=0; i<M4.nrow; ++i) {
        for(j=0; j<M4.ncol; ++j) {
            d(i+M1.nrow,j+M3.ncol) = M4(i,j);
        }
    }
    return d;
}

// Create the Hadamard product between two matrices
db_matrix db_matrix::hadamard(const db_matrix &M1, const db_matrix &M2)
{
    int i,j; 
    int ncol=M1.ncol;

    assert(M1.nrow==M2.nrow);

    db_matrix d(M1.nrow, ncol);

    // Loop on all points
    for(i=0; i<M1.nrow; ++i) {
        for(j=0; j<M1.ncol; ++j) {
            d(i,j) = M1(i,j)*M2(i,j);
        }
    }
    return d;
}

// Create the Hadamard product between two matrices
// Overwrite the actual matrix with the result
db_matrix &db_matrix::hadamard(const db_matrix &M2)
{
    int i,j;

    assert(nrow==M2.nrow);
    assert(ncol==M2.ncol);

    // Loop on all points
    for(i=0; i<nrow; ++i) {
        for(j=0; j<ncol; ++j) {
            operator()(i,j) *= M2(i,j);
        }
    }
    return *this;
}

// Print the given matrix
void db_matrix::printMatrix(void)
{
    int i, j;

    for(i=0; i<nrow; ++i) {
        for(j=0; j<ncol; ++j) {
            cout <<" "<< setw(14)<<operator()(i,j)<<"  ";
        }
        cout <<"\n";
        cout.flush();
    }
}

/** Print the given matrix on a file whose name is indicated.
*/
void db_matrix::printMatrixOnFile(string filename)
{
    int i, j;
    ofstream fout;
    fout.open(filename.c_str());

    for(i=0; i<nrow; ++i) {
        for(j=0; j<ncol; ++j) {
            fout <<" "<< setw(14)<<operator()(i,j)<<"  ";
        }
        fout <<"\n";
    }
    fout.close();
}


// The same as Matlab!
// DEVELOPMENT NOT COMPLETE!!!
/*db_matrix db_matrix::fftshift()
{

    // Matrix size calculation.
    // The variables dimx and dimy are the number of harmonics taken in x
    // and in y, by taking into account only zero and positive frequencies.

    int size_x=ncol;
    int size_y=nrow;
    int dimx = int(floor(double(size_x)/2.0));
    int dimy = int(floor(double(size_y)/2.0));
    int shiftx = dimx;
    int shifty = dimy;

    db_matrix R(nrow, ncol);

    int i,j;

    ////
    for(i=0; i<dimy; ++i) {
        for(j=0; j<dimx; ++j) {
            // Swap the first and the third quadrant
            R(i,j)=operator()(i+shifty,j+shiftx);
            R(i+shifty,j+shiftx)=operator()(i,j);
            // Swap the second and the fourth quadrant
            R(i+shifty,j)=operator()(i,j+shiftx);
            R(i,j+shiftx)=operator()(i+shifty,j);

        }
    }
    ///
    if(size_x & 1)  {
        for(i=0; i<dimy; ++i)  {
            R(i+shifty, shiftx)= operator()(i,0);
            R(i, shiftx)= operator()(i+shifty+1,0);
        }
        if (size_y & 1) {
            R(dimy+shifty, shiftx)= operator()(dimy,0);
        }
    }
    if(size_y & 1) {
        for(j=0; j<dimx; ++j) {
            R(shifty,j+shiftx)= operator()(0,j);
            R(shifty,j)= operator()(0,j+shiftx+1);
        }
        if (size_x & 1) {
            R(shifty,dimx+shiftx)= operator()(0,dimx);
        }
    }

    return R;
}///
// The same as Matlab!
// DEVELOPMENT NOT COMPLETE!!!
////db_matrix db_matrix::ifftshift()
{

    // Matrix size calculation.
    // The variables dimx and dimy are the number of harmonics taken in x
    // and in y, by taking into account only zero and positive frequencies.

    int size_x=ncol;
    int size_y=nrow;
    int dimx = int(floor(double(size_x)/2.0));
    int dimy = int(floor(double(size_y)/2.0));
    int shiftx = dimx;
    int shifty = dimy;

    bool xodd=size_x & 1;
    bool yodd=size_y & 1;

    db_matrix R(nrow, ncol);
    int i,j;

    if(yodd) {
        for(j=0; j<dimx; ++j) {
            R(shifty,j+shiftx+xodd)= operator()(shifty,j);
            R(shifty,j+xodd)= operator()(shifty,j+shiftx+1);
        }
        if (xodd) {
            R(shifty,0)= operator()(shifty,shiftx);
        }
    }

    if(xodd) {
        for(i=0; i<dimy; ++i) {
            R(i+shifty+yodd, shiftx)= operator()(i, shiftx);
            R(i+yodd, shiftx)= operator()(i+shifty+1, shiftx);
        }
        if (yodd) {
            R(0, shiftx)= operator()(shifty,shiftx);
        }
    }

    return R;
}*/

long db_matrix::getSize()
{
    if(A==NULL) {
        return sizeof(db_matrix);
    } else {
        return sizeof(db_matrix)+ncol*nrow*sizeof(complex<double>);
    }
}

// Transpose the matrix
db_matrix db_matrix::transpose()
{
    int i,j;

    db_matrix d(ncol, nrow);

    // Loop on all points

    for(i=0; i<nrow; ++i) {
        for(j=0; j<ncol; ++j) {
            d(j,i) = operator()(i,j);
        }
    }

    return d;
}

// Invert the matrix. Call LAPACK routines ZGETRF and ZGETRI
// The original matrix is overwritten with the new results.
// It returns the inverted matrix.
db_matrix &db_matrix::invert()
{
    integer info;
    integer M=nrow;
    integer N=ncol;
    integer *pivot = new integer[(nrow<ncol)?nrow:ncol];

    /* WARNING! This code makes use of violent casts between the C++
       complex<double> template and the double complex type of Fortran.
       Even if this can be dangerous, it is normally relatively safe to admit
       that the STL template is coded in memory as the real and imaginary part
       in double.
       This works well on x86 machines, gcc and g77 compilers.
    */
    doublecomplex *B;
    doublecomplex *w = new doublecomplex[nrow];
    B= (doublecomplex*) A;
    assert(ncol==nrow);

    // At first, perform a LU decomposition
    zgetrf(&M, &N, B, &M, pivot, &info);

    // Then, obtain the inverse matrix
    zgetri(&N, B, &N, pivot, w, &N, &info);
    delete[] pivot;
    delete[] w;

    return *this;
}

// Calculate eigenvalues and eigenvectors. Call the LAPACK ZGEEV routine
// If copy is true, the original matrix is unchanged and the eigenvalues only
// are given.
// The original matrix is overwritten with the eigenvectors.
// It returns the eigenvalues vectors.
// If eigvect!=NULL, calculates also the eigenvectors and store them (in the
// same order as eigenvalues) in the eigvect matrix, which will be initizialized
db_matrix &db_matrix::eig(db_matrix *eigvect)
{
    integer info;
    integer N=ncol;

    /* WARNING! This code makes use of violent casts between the C++
       complex<double> template and the double complex type of Fortran.
       Even if this can be dangerous, it is normally relatively safe to admit
       that the STL template is coded in memory as the real and imaginary part
       in double.
       This works well on x86 machines, gcc and g77, gfortran compilers.
    */
    doublecomplex *B;
    doublecomplex *V = new doublecomplex[nrow];

    if (V==NULL) {
        throw db_matrix_error("Error: could not allocate memory!");
    }
    doublecomplex zero={0.0,0.0};

    // Probably this is not needed. However, valgrind complains that there
    // is some "Conditional jump or move depends on uninitialised value(s)"
    // if this code is not executed. I trust it much more than my naive
    // intuition!
    for(int i=0; i< nrow; ++i)
        V[i]=zero;

    char jobvl = 'N';
    char jobvr = 'V';

    assert (ncol==nrow);

    B= (doublecomplex*) A;
    assert(ncol==nrow);

    integer LDVR = nrow;
    integer LDVL = 1;
    integer NB = 64;
    integer LWORK = (NB+1)*nrow;
    doublecomplex *WORK = new doublecomplex[(NB+1)*nrow];

    if (WORK==NULL) {
        // V has been allocated.
        delete[] V;
        throw db_matrix_error("Error: could not allocate memory!");
    }

    for(int i=0; i< (NB+1)*nrow; ++i)
        WORK[i]=zero;

    doublereal *RWORK = new doublereal[3*nrow];
    if (RWORK==NULL) {
        // V and WORK have been allocated.
        delete[] V;
        delete[] WORK;
        throw db_matrix_error("Error: could not allocate memory!");
    }

    for(int i=0; i< 3*nrow; ++i)
        RWORK[i]=0.0;

    doublecomplex *VL=NULL;
    doublecomplex *VR=NULL;

    if(eigvect==NULL)
        jobvr ='N';
    else {
        VR= new doublecomplex[nrow*ncol];

        #ifdef MEMORY_TRACK
            //cout<<"File "<<__FILE__<<" at line "<<__LINE__<<"\n";
            //cout<<"Matrix allocated of size "<<nrow<<"x"<<ncol<<"\n";
            total_elements+=nrow*ncol;
            if(total_elements>max_elements) max_elements=total_elements;
            //cout << "Total number of elements in memory: "
            //  <<total_elements<<"\n";
        #endif
        if (VR==NULL) {
            // V and WORK have been allocated.
            delete[] RWORK;
            delete[] V;
            delete[] WORK;
            throw db_matrix_error("Error: could not allocate memory!");
        }
        for(int i=0; i< nrow*ncol; ++i)
            VR[i]=zero;
    }
    if(!ZGEEVX) {
        zgeev(&jobvl, &jobvr, &N, B, &N, V, VL, &LDVL, VR, &LDVR, WORK, &LWORK,
            RWORK, &info);
    } else {
        integer ilo, ihi;
        double *scale = new double[N];
        double *rconde = new double[N];
        double *rconv = new double[N];

        for(int i=0; i< N; ++i) {
            scale[i]=0.0;
            rconde[i]=0.0;
            rconv[i]=0.0;
        }

        double abnrm;

        char cn = 'N';

        zgeevx(&cn,&jobvl, &jobvr, &cn, &N,
            B, &N, V, VL, &LDVL, VR, &LDVR, &ilo, &ihi,
            scale,  &abnrm, rconde, rconv,
            WORK, &LWORK,   RWORK, &info);
        delete[] scale;
        delete[] rconde;
        delete[] rconv;
    }

    if (eigvect!=NULL) {
        if(eigvect->A!=NULL) {
            delete[] eigvect->A;
            eigvect->A=NULL;

            #ifdef MEMORY_TRACK
                total_elements-=eigvect->nrow*eigvect->ncol;
            #endif
        }
        eigvect->A=(complex<double>*) VR;
        eigvect->ncol = ncol;
        eigvect->nrow = nrow;

    } else {
        delete [] VR;
    }

    delete[] V;
    delete[] RWORK;
    delete[] WORK;

    // In Fortran, ints are coded on 32 bits. If this code is compiled in
    // a 64 bits system, ints will be 64 bits, thus displaying rubbish.
    info=info & 0xFFFFFFFFl;

    cout << "LAPACK ZGEEVX executed, info = " << info << "\n";

    return *this;
}

// Calculate generalized eigenvalues and eigenvectors, as defined by the
// equation A*u=lambda*B*u, where A and B are two square matrix, lambda is the
// associated eigenvalue, u is a column vector representing the eigenvector.
// Calls the LAPACK routine ZGGEVX.
// eigval will contain the generalized eigenvalues, whereas a matrix
// containing eivenvectors as column vectors stored in the same order as
// eigenvalues is returned.

db_matrix &db_matrix::eigGen(db_matrix matrixB, db_matrix *eigval)
{
    integer info;
    integer N=ncol;

    /* WARNING! This code makes use of violent casts between the C++
       complex<double> template and the double complex type of Fortran.
       Even if this can be dangerous, it is normally relatively safe to admit
       that the STL template is coded in memory as the real and imaginary part
       in double.
       This works well on x86 machines, gcc and g77 compilers.
    */

    doublecomplex *AA;
    doublecomplex *V = new doublecomplex[nrow];

    if (V==NULL) {
        throw db_matrix_error("Error: could not allocate memory!");
    }

    doublecomplex zero={0.0,0.0};

    doublecomplex *alpha = new doublecomplex[nrow];

    if (alpha==NULL) {
        delete[] V;
        throw db_matrix_error("Error: could not allocate memory!");
    }
    doublecomplex *beta = new doublecomplex[nrow];
    if (beta==NULL) {
        delete[] V;
        delete[] alpha;
        throw db_matrix_error("Error: could not allocate memory!");
    }

    // Probably this is not needed. However, valgrind complains that there
    // is some "Conditional jump or move depends on uninitialised value(s)"
    // if this code is not executed. I trust it much more than my naive
    // intuition!
    for(int i=0; i< nrow; ++i) {
        V[i]=zero;
        alpha[i]=beta[i]=zero;
    }

    assert (ncol==nrow);

    AA = (doublecomplex*) A;

    doublecomplex *BB =  (doublecomplex*) matrixB.A;

    assert(ncol==nrow);

    integer LDVR = nrow;
    integer LDVL = 1;
    integer NB = 64;
    integer LWORK = (NB+1)*nrow;
    doublecomplex *WORK = new doublecomplex[(NB+1)*nrow];

    if (WORK==NULL) {
        delete[] V;
        delete[] alpha;
        delete[] beta;
        throw db_matrix_error("Error: could not allocate memory!");
    }

    for(int i=0; i< (NB+1)*nrow; ++i)
        WORK[i]=zero;

    doublereal *RWORK = new doublereal[3*nrow];
    if (RWORK==NULL) {
        delete[] V;
        delete[] alpha;
        delete[] beta;
        delete[] WORK;
        throw db_matrix_error("Error: could not allocate memory!");
    }

    for(int i=0; i< 3*nrow; ++i)
        RWORK[i]=0.0;

    doublecomplex *VL=NULL;
    doublecomplex *VR=NULL;

    if(eigval==NULL) {
        delete[] V;
        delete[] alpha;
        delete[] beta;
        delete[] WORK;
        delete[] RWORK;
        throw db_matrix_error("Error: you should give eigval.");
    } else {
        VR= new doublecomplex[nrow*ncol];
        if (VR==NULL) {
            delete[] V;
            delete[] alpha;
            delete[] beta;
            delete[] WORK;
            delete[] RWORK;
            throw db_matrix_error("Error: could not allocate memory!");
        }

        int nelem=nrow*ncol;
        for(int i=0; i< nelem; ++i)
            VR[i]=zero;

        #ifdef MEMORY_TRACK
            total_elements+=nrow*ncol;
            if(total_elements>max_elements) max_elements=total_elements;
        #endif

    }

    integer ilo, ihi;
    double *lscale = new double[N];
    double *rscale = new double[N];

    double *rconde = new double[N];
    double *rcondv = new double[N];

    bool *bwork = new bool[N];
    integer *iwork = new integer[N+2];

    for(int i=0; i< N; ++i) {
        lscale[i]=rscale[i]=0.0;
        rconde[i]=0.0;
        rcondv[i]=0.0;
    }

    double abnrm, bbnrm;

    zggevx("N", "N","V", "N", &N,
        AA, &N, BB, &N, alpha, beta,
        VL, &LDVL, VR, &LDVR, &ilo, &ihi,
        lscale, rscale, &abnrm, &bbnrm,
        rconde, rcondv,
        WORK, &LWORK, RWORK, iwork, bwork, &info);


    delete[] lscale;
    delete[] rscale;

    delete[] rconde;
    delete[] rcondv;

    delete[] iwork;
    delete[] bwork;
    delete[] RWORK;
    delete[] WORK;

    if (eigval!=NULL){
        if(eigval->A!=NULL) {
            delete[] eigval->A;
            eigval->A=NULL;

            #ifdef MEMORY_TRACK
                total_elements-=eigval->nrow*eigval->ncol;
            #endif
        }
        eigval->A = new complex<double>[nrow*ncol];

        // In reality, LAPACK documentation explicitly warns against this
        // naive way of calculating eigenvalues, since beta[i] might be
        // equal to zero quite easily.
        // See http://www.netlib.org/lapack/lug/node35.html

        for(int i=0; i<N; ++i) {
            complex<double> *aleph = (complex<double> *)&alpha[i];
            complex<double> *beth = (complex<double> *)&beta[i];

            if((*beth)!=0.0) {
                eigval->A[i*N+i]=(*aleph)/(*beth);
            } else
                eigval->A[i*N+i]=NAN;
        }
        eigval->ncol = ncol;
        eigval->nrow = nrow;

    } else {
    }

    delete[] V;
    delete[] alpha;
    delete[] beta;
    delete[] A;
    A=(complex<double>*)VR;

    double max=0;

    // No information is given concerning eigenvalues normalization.
    // So we normalize them.

    int i,j;

    for(j=0; j<getNcol(); ++j) {
        max=0;
        for (i=0; i<getNrow();++i) {
            if(abs(operator()(i,j)) > max) {
                max=abs(operator()(i,j));
            }
        }
        for (i=0; i<getNrow();++i) {
            operator()(i,j)/=max;
        }
    }

    // In Fortran, ints are coded on 32 bits. If this code is compiled in
    // a 64 bits system, ints will be 64 bits, thus displaying rubbish.
    info=info & 0xFFFFFFFFl;

    cout << "LAPACK ZGGEVX executed, info = " << info << "\n";

    return *this;
}


/** Execute a fast Fourier transform using FFTW3.
    If type is false, perform a FORWARD FFT.
    If type is true, perform a BACKWARD FFT.
*/
db_matrix db_matrix::fft_a(bool type)
{
    int nx = nrow;
    int ny = ncol;
    db_matrix d(nrow, ncol);

    #ifdef USE_KISSFFT
    // Define dimensions (e.g., 8x8 2D FFT)
    int dims[2] = {nx, ny};
    int ndims = 2;
    
    // Allocate configuration
    kiss_fftnd_cfg cfg = kiss_fftnd_alloc(dims, ndims, type?0:1, NULL, NULL);
    
    // Allocate input and output buffers
    kiss_fft_cpx *in = (kiss_fft_cpx *)malloc(2*nx * ny * sizeof(kiss_fft_cpx));
    kiss_fft_cpx *out= (kiss_fft_cpx *)malloc(2*nx * ny * sizeof(kiss_fft_cpx));
    
    // Populate input data...
    int row;
    complex<double> v;
    for (unsigned int i = 0; i < nx; ++i){
        row=i*ny;
        for (unsigned int j = 0; j < ny; ++j){
            v=operator()(i,j);
            in[row+j].r = v.real();
            in[row+j].i = v.imag();
        }
    }
    
    // Perform FFT
    kiss_fftnd(cfg, in, out);
    // Clean up
    //free(cfg);
    //free(in);
    for (int i = 0; i < nx; ++i){
        row=i*ny;
        for (int j = 0; j < ny; ++j){
            d(i,j)=complex<double>(out[row+j].r,out[row+j].i);
        }
    }
    //free(out);
    

    #elif defined(USE_FFTW)

    fftw_complex *in;
    fftw_complex *out;
    fftw_plan plan;

    sem_wait (mutex_fftw);  /*- sync start -*/
    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*nx*ny);
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*nx*ny);
    sem_post (mutex_fftw);  /*- sync end -*/

    int row;
    complex<double> v;
    for (unsigned int i = 0; i < nx; ++i){
        row=i*ny;
        for (unsigned int j = 0; j < ny; ++j){
            v=operator()(i,j);
            in[row+j][0] = v.real();
            in[row+j][1] = v.imag();
        }
    }

    sem_wait (mutex_fftw);  /*- sync start -*/
    if(type) {
        plan = fftw_plan_dft_2d(nx, ny, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    } else {
        plan = fftw_plan_dft_2d(nx, ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    sem_post (mutex_fftw);  /*- sync end -*/

    fftw_execute(plan);

    sem_wait (mutex_fftw);  /*- sync start -*/
    fftw_free (in);
    fftw_destroy_plan(plan);
    sem_post (mutex_fftw);  /*- sync end -*/

    for (int i = 0; i < nx; ++i){
        row=i*ny;
        for (int j = 0; j < ny; ++j){
            d(i,j)=complex<double>(out[row+j][0],out[row+j][1]);
        }
    }

    sem_wait (mutex_fftw);  /*- sync start -*/
    fftw_free (out);
    sem_post (mutex_fftw);/*- sync end -*/
    #else
        cout <<"FFTW library not included. I can not compute fft2\n";
    #endif
    return d;
}

/** Adapt the size of the actual matrix to the value specified.
    If the wanted size is bigger than the size of the actual matrix,
    a zero pad strategy will be applied. If it is the opposite case, the
    harmonics will be truncated to the wanted size.
*/
db_matrix db_matrix::zero_pad(int si, int sj)
{
    db_matrix grand(si,sj);
    // All elements are set at zero while initizializing

    int limitrow = min(nrow, si);
    int limitcol = min(ncol, sj);

    for (unsigned int i=0; i<limitrow/2+1; ++i) {
        for (unsigned int j=0; j<limitcol/2+1; ++j) {
            grand(i,j)=operator()(i,j);
            if(j<(limitcol-1)/2)
                grand(i,sj - j - 1)=operator()(i,ncol - j - 1);
            if(i<(limitrow-1)/2)
                grand(si - i - 1,j)=operator()(nrow - i - 1,j);
            if(j<(limitcol-1)/2 && i<(limitrow-1)/2)
                grand(si - i - 1, sj - j - 1)=operator()(nrow-i-1,ncol- j - 1);
        }
    }

    return grand;
}

// If called when all matrices are supposed to be eliminated, this gives the
// number of remaining elements.
void db_matrix::leaks(void)
{
    #ifdef MEMORY_TRACK
        cout<<"File "<<__FILE__<<" at line "<<__LINE__<<"\n";
        cout << "Leaked elements in memory: "<<total_elements
                <<"\n";
    #endif
}

//// JEROME:


// The same as Matlab!
db_matrix db_matrix::fftshift()
{

    // Matrix size calculation.
    // The variables dimx and dimy are the number of harmonics taken in x
    // and in y, by taking into account only zero and positive frequencies.

    int size_x=ncol;
    int size_y=nrow;
    int dimx = int(floor(double(size_x)/2.0));
    int dimy = int(floor(double(size_y)/2.0));
    int shiftx = dimx;
    int shifty = dimy;

    db_matrix R(nrow, ncol);

    int i,j;


    // if the size of the matrix is even or odd, the algorithm changed a little
    // bit the term parityx and parityy takes this into account

    int parityx = ncol & 1;
    int parityy = nrow & 1;

    for(j=0; j<dimx; ++j) {
        // Swap the first and the third quadrant
        for(i=0; i<dimy; ++i)
            R(i,j)=operator()(i+shifty+parityy,j+shiftx+parityx);
        // Swap the second and the fourth quadrant
        for(i=0; i<dimy+parityy; ++i)
            R(i+shifty,j)=operator()(i,j+shiftx+parityx);
    }


    for(j=0; j<dimx+parityx; ++j) {
        // Swap the first and the third quadrant
        for(i=0; i<dimy; ++i)
            R(i,j+shiftx)=operator()(i+shifty+parityy,j);
        // Swap the second and the fourth quadrant
        for(i=0; i<dimy+parityy; ++i)
            R(i+shifty,j+shiftx)=operator()(i,j);

    }

    return R;
}

// The same as Matlab!
db_matrix db_matrix::ifftshift()
{

    // Matrix size calculation.
    // The variables dimx and dimy are the number of harmonics taken in x
    // and in y, by taking into account only zero and positive frequencies.

    int size_x=ncol;
    int size_y=nrow;
    int dimx = int(floor(double(size_x)/2.0));
    int dimy = int(floor(double(size_y)/2.0));
    int shiftx = dimx;
    int shifty = dimy;

    db_matrix R(nrow, ncol);

    int i,j;

    // if the size of the matrix is even or odd, the algorithm change a little
    // bit the term parityx and parityy take this into account

    int parityx = ncol & 1;
    int parityy = nrow & 1;

    for(j=0; j<dimx; ++j) {
        // Swap the first and the third quadrant
        for(i=0; i<dimy; ++i)
            R(i+shifty+parityy,j+shiftx+parityx)=operator()(i,j);
        // Swap the second and the fourth quadrant
        for(i=0; i<dimy+parityy; ++i)
            R(i,j+shiftx+parityx)=operator()(i+shifty,j);
    }


    for(j=0; j<dimx+parityx; ++j) {
        // Swap the first and the third quadrant
        for(i=0; i<dimy; ++i)
            R(i+shifty+parityy,j)=operator()(i,j+shiftx);
        // Swap the second and the fourth quadrant
        for(i=0; i<dimy+parityy; ++i)
            R(i,j)=operator()(i+shifty,j+shiftx);

    }


    return R;
}

// copy a smaller matrix M into the larger matrix M1
// M1.copy(M,i,j);
// Where i and j designate the place where to copy
db_matrix db_matrix::copy(const db_matrix &M, int i, int j)
{
    int size_x = M.ncol;
    int size_y = M.nrow;

    int p;
    int q;
    for (p = 0; p < size_y; ++p) {
        for (q = 0; q < size_x; ++q) {
            operator()(i+p,j+q) = M(p,q);
        }
    }
    return *this;
}
// add a part of a matrix M into the matrix M1
// Where mi and mj designate from which we have to copy
// Where ni and nj designate the place where to copy
// Where size_x and size_y are the size of the element to copy
// sign is a coefficient which is multiplied to the elements
db_matrix db_matrix::add(const db_matrix &M, int sign ,int ni, int nj,
            int mi, int mj,int size_x, int size_y)
{

    int p;
    int q;
    complex<double> s=complex<double>(sign,0);

    for (p = 0; p < size_y; ++p) {
        for (q = 0; q < size_x; ++q) {
            operator()(ni+p,nj+q) += s * M(p+mi,q+mj);
        }
    }
    return *this;
}

// return the absolute value of an integer:
inline int abs(int i)
{
    if (i < 0 )
        i = -i;

    return i;
}


// Enroll the vector comming from fft transform depending on the symmetry
// with a given offset called shift
// syntaxe: M_fft.enrollvector(xsymmetric, ysymmetric,);
db_matrix db_matrix::fft2vector(bool xsymmetric, bool ysymmetric, int shift)
{

    // Matrix size calculation.
    // The variables dimx and dimy are the number of harmonics taken in x
    // and in y, by taking into account only zero and positive frequencies

    int size_x=ncol;
    int size_y=nrow;
    int dimx;
    int dimy;

    int shiftx;
    int shifty;

    // the size of the matrix is reduced in the case of symmetry
    if (xsymmetric) {
        dimx = int(floor(double(size_x)/2.0) + 1.0);
        shiftx = int(floor(double(size_x)/2.0));
    }
    else {
        dimx = size_x;
        shiftx =0;
    }
    // the size of the matrix is reduced in the case of symmetry
    if (ysymmetric) {
        dimy = int(floor(double(size_y)/2.0) + 1.0);
        shifty = int(floor(double(size_y)/2.0));;
    }
    else {
        dimy = size_y;
        shifty = 0;
    }

    // order the Fourier transform coefficient
    db_matrix ordered = fftshift();

    //index
    int i,j;

    db_matrix R(2*dimy*dimx,1);

    for(i=0; i<dimy; ++i) {
        for(j=0;j<dimx;++j) {
            R(i*dimx+j+shift,0)=ordered(i+shifty,j+shiftx);
        }
    }

    return R;
}

/** Convert a vector to fft transform matrix depending on the symmetry
    This is an unrolling routine. No FFT calculation is done yet, but the
    terms of a 1D vector are placed in a 2D matrix, at the right place
    for the 2D FFT.

    if applyShift is true, we will output the y component
    syntax: fft_M = Vector.enrollvector(symx, symy, snux, snuy, applyShift,
        calcH);
*/
db_matrix db_matrix::vector2fft(symmetry_e symx, symmetry_e symy,
    int snux, int snuy, bool applyShift, bool calcH, bool calcz)
{

    // Matrix size calculation.
    // The variables dimx and dimy are the number of harmonics taken in x
    // and in y, by taking into account only zero and positive frequencies

    int size_x=snux;
    int size_y=snuy;
    int dimx;
    int dimy;

    int shiftx;
    int shifty;
    int shift;

    int factorX=1;
    int factorY = 1;
    bool xsymmetry;
    bool ysymmetry;

    // the size of the matrix is reduced in the case of symmetry
    if (symx == symmetric || symx == anti_symmetric)    {
        dimx = int(floor(double(size_x)/2.0) + 1.0);
        shiftx = int(floor(double(size_x)/2.0) );
        xsymmetry = true;
        if (calcz) {
            // Ez is odd
            if (symx == anti_symmetric)
                factorX = -1;
        } else {
            // Ex is odd
            if (symx == symmetric && !applyShift)
                factorX = -1;
            // Ey is odd
            if (symx == anti_symmetric && applyShift)
                factorX = -1;
        }
        // H is inverted with respect to E
        if (calcH)
            factorX *= -1;
    } else {
        dimx = size_x;
        shiftx =0;
        xsymmetry = false;
    }
    // the size of the matrix is reduced in the case of symmetry
    if (symy == symmetric || symy == anti_symmetric)    {
        dimy = int(floor(double(size_y)/2.0) + 1.0);
        shifty = int(floor(double(size_y)/2.0));
        ysymmetry = true;
        if (calcz) {
            // Ez is odd
            if (symy == anti_symmetric)
                factorY = -1;
        } else {
            // Ey is odd
            if (symy == symmetric && applyShift)
                factorY = -1;
            // Ex is odd
            if (symy == anti_symmetric && !applyShift)
                factorY = -1;
        }
        // H is inverted with respect to E
        if (calcH)
            factorY *= -1;

    } else {
        dimy = size_y;
        shifty = 0;
        ysymmetry = false;
    }

    // computing the y component
    if (applyShift)
        shift = dimx*dimy;
    else
        shift = 0;


    //index
    int i,j;

    db_matrix R(snuy,snux);

    for(i=0; i<dimy; ++i) {
        for(j=0;j<dimx;++j) {

            // forth quadrant
            R(i+shifty,j+shiftx)=operator()(i*dimx+j+shift,0);

            if (ysymmetry) {
                // first quadrant
                R(shifty-i,shiftx+j)= complex<double>(factorY,0) *
                                operator()(i*dimx+j+shift,0);
            }

            if (xsymmetry && ( shiftx-j >= 0 && j != 0) ) {

                // third quadrant
                R(i+shifty,shiftx-j)= complex<double>(factorX,0) *
                                    operator()(i*dimx+j+shift,0);

                if (xsymmetry && (shifty-i >= 0 && i != 0)) {
                    // second quadrant
                    R(shifty-i,shiftx-j)= complex<double>(factorX*factorY,0) *
                             operator()(i*dimx+j+shift,0);
                }
            }
        }
    }
    return R.ifftshift();
}

// return the normal field at a given px,py coordinate
complex<double> db_matrix::expectedNormalField(double px, double py, double tot_x, double tot_y)
{
    complex<double> n;

        int posx = round((px/tot_x+0.5)*getNcol());
        int posy = round((py/tot_y+0.5)*getNrow());

        if (posx<0 || posy<0 ||
            posx>=getNcol() || posy>=getNrow()) {
            return -1;
        }
        n = operator()(posy, posx);

    return n;
}
