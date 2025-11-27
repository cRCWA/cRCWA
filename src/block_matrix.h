/***************************************************************************
*     CLASS db_matrix                                                      *
*     Davide Bucci, CROMA     2008-present                                 *
*     Jérôme Michallon, CROMA     2012                                     *
*     MINATEC-INPG, 3, parvis Luis Neel                                    *
*     38016, Grenoble CEDEX, France                                        *
*                                                                          *
*     bucci@minatec.inpg.fr                                                *
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

#ifndef BLOCK_MATRIX_H
#define BLOCK_MATRIX_H
// Comment to check for assertion in the code. This should be useful for debug
#define NDEBUG

// Comment to check for assertion in the code. This should be useful for debug
#define NDEBUG

#include <complex>
#include <cassert>
#include <iostream>

// Types of symmetry considered
typedef enum symmetry_e_t{
    anti_symmetric, symmetric, no_symmetry, previous_version
} symmetry_e;

using namespace std;


void init_semaphore_FFTW(void);
void delete_semaphore_FFTW(void);

class db_matrix {
private:
    complex<double>* A;
    int ncol;
    int nrow;

    static double sinc(double alpha) {return (alpha==0)?1:sin(alpha)/alpha;}
    db_matrix fft_a(bool type);

public:
    // Constructors ***********************************************************

    // Standard constructor: fabricate a null matrix in which the A pointer is
    // not defined.
    db_matrix(void);

    // Standard constructor: fabricate a matrix of the given size
    db_matrix(int col, int row);

    // Copy constructor: element by element copy of the given matrix
    db_matrix(const db_matrix& D);

    // Destructor *************************************************************
    ~db_matrix(void);

    // Static members *********************************************************

    // Create a test matrix of the given size
    static db_matrix createTestMatrix(const int nrow, const int ncol);

    // Create a matrix full of ones, of the given size
    static db_matrix createOnesMatrix(const int nrow, const int ncol);

    // Create an unit matrix
    static db_matrix createUnitMatrix(const int nrow, const int ncol);

    // Bloc-merge two matrices on the same row. They must have the same number
    // of rows
    static db_matrix mergeMatrixOnRow(const db_matrix &M1, const db_matrix &M2);

    // Bloc-merge two matrices on the same column. They must have the same
    // number of columns
    static db_matrix mergeMatrixOnColumn(const db_matrix &M1,
        const db_matrix &M2);
    static db_matrix mergeMatrixQuad(const db_matrix &M1,
        const db_matrix &M2,const db_matrix &M3,
        const db_matrix &M4);
    // Create the Hadamard product between two matrices
    static db_matrix hadamard(const db_matrix &M1, const db_matrix &M2);

    // Calculate the Fourier transform of a rectangle and store it in the usual
    // FFT way.
    static db_matrix rectangleFT(double sx,// Total size in the x dir.
        double sy,              // Total size in the y direction
        double periodx,         // x-period (size) of the window
        double periody,         // y-period (size) of the window
        int Nx,                 // number of x harmonics
        int Ny,                 // number of y harmonics
        double centerx,         // center of the rectangle in the window: x
        double centery);        // center of the rectangle in the window: y

    static void leaks(void);

   // Non-static members *****************************************************

    void printSize(void)
    {
        cout<<"Matrix formed by "<<nrow<<" rows and "<<ncol<<" columns."<<endl;
    }
    // Create the Hadamard product between two matrices.
    // Overwrite the actual matrix with the results.
    db_matrix &hadamard(const db_matrix &M2);

    db_matrix toeplitz_ones(double nux,
        double nuy,
        int p1x,
        int p1y,
        int p2x,
        int p2y);

    db_matrix fftshift();
    db_matrix ifftshift();

    long getSize();

     // Modified derivative calculation
    db_matrix toeplitz_deriv(double nux,
        double nuy,
        int p1x,
        int p1y,
        double ksinthetax,
        double ksinthetay);

    db_matrix windowing(void);

    db_matrix &scalemult(const complex<double> alpha, const db_matrix &d);

    // Explicitely eliminate a matrix to save space.
    void kill(void);

    // Print the matrix on the screen.
    void printMatrix(void);
    // Print the matrix on a file.
    void printMatrixOnFile(string filename);


    db_matrix toeplitz_mod(double nux,
        double nuy,
        int p1x,
        int p1y,
        int p2x,
        int p2y);

    //JEROME :
    // Creation of 1D toeplitz matrix
    db_matrix toeplitz1D(int elementsize,
         int colomnVector, bool symmetry, int alpha);
    // for copying a small matrix in a larger one
    db_matrix copy(const db_matrix &M, int i, int j);
    // for adding a small matrix to a larger one
    db_matrix add(const db_matrix &M, int sign ,int ni, int nj,
            int mi, int mj,int size_x, int size_y);
    // compute the toeplitz matrices in case of symmetry
    db_matrix toeplitz_sym( bool xsymmetric, bool ysymmetric,
            int alpha, int beta, bool IsOddx, bool IsOddy);
    // compute the derivative matrices in case of symmetry
    db_matrix toeplitz_deriv_sym(bool xsymmetric, bool ysymmetric,
                double nux, double nuy, int p1x, int p1y);
    // enroll the vector comming from a fft transform
    db_matrix fft2vector(bool xsymmetric, bool ysymmetric, int shift);
    // convert a vector the vector  to fft transform matrix
    db_matrix vector2fft(symmetry_e symx, symmetry_e symy, int snux, int snuy,
        bool applyShift, bool calcH, bool calcz);
    // sum different toeplitz matrices in case of symmetry
    db_matrix &toeplitz1D_sum(db_matrix &Toeplitz1Dp, int elementsize,
                 int gamma, bool IsOdd);

    // return the normal field at a given px,py coordinate
    complex<double> expectedNormalField(double px, double py,
        double tot_x, double tot_y);

    // Obtain the number of rows
    int getNrow(void) {return nrow;}

    // Obtain the number of columns
    int getNcol(void) {return ncol;}

    // operators overloading
    db_matrix& operator=(const db_matrix &d);

    db_matrix& operator*=(const complex<double> &c);
    db_matrix& operator/=(const complex<double> &c);

    db_matrix& operator*=(const db_matrix &d);
    db_matrix& operator+=(const db_matrix &d);
    db_matrix& operator-=(const db_matrix &d);
    db_matrix operator-();

    db_matrix zero_pad(int si, int sj);

    inline complex<double>& operator()(const int i, const int j)
    {
        /*if(i>=nrow || j>=ncol) {
            cerr << i << " " << nrow << "  " << j << " " << ncol << "\n";
        }*/
        assert (i<nrow && j<ncol);
        assert (i>=0 && j>=0);
        assert (A!=NULL);

        // Data are row-based as in Fortran
        return A[i+j*nrow];
    }

    inline complex<double>& operator()(const int i, const int j) const
    {
        /*if(i>=nrow || j>=ncol) {
            cerr << i << " " << nrow << "  " << j << " " << ncol << "\n";
        }*/
        assert (i<nrow && j<ncol);
        assert (i>=0 && j>=0);
        assert (A!=NULL);

        // Data are row-based as in Fortran
        return A[i+j*nrow];
    }

    /* Invert the matrix. Call LAPACK routines ZGETRF and ZGETRI
       The original matrix is overwritten with the new results.
       It returns the inverted matrix.
    */
    db_matrix &invert();

    /* Invert the matrix. Call the LAPACK ZGEEV routine
       The original matrix is overwritten with the eigenvectors.
       It returns the pointer to the eigenvalues vectors.
    */
    db_matrix &eig(db_matrix *eigvect);

    /* Calculate generalized eigenvalues and eigenvectors, as defined by the
       equation A*u=lambda*B*u, where A and B are two square matrix, lambda is
       the associated eigenvalue, u is a column vector representing the
       eigenvector. Calls the LAPACK routine ZGGEVX.
    */
    db_matrix &eigGen(db_matrix matrix2, db_matrix *eigvect);

    /* Calculates the 2D Fast Fourier Transform of the given matrix.
       Employs the FFTW3 library
    */
    inline db_matrix fft2(void) {return fft_a(false);}

    // Transpose the matrix.
    db_matrix transpose();

    /* Calculates the 2D Inverse Fast Fourier Transform of the given matrix.
       Employs the FFTW3 library
    */
    inline db_matrix ifft2(void){return fft_a(true);}

    /* Returns true if the matrix is empty.

    */
    bool isEmpty(void){return A==NULL;}


    friend db_matrix &multsum(const db_matrix &AA, const db_matrix &B,
        db_matrix &C);

    friend db_matrix &multsumscale(const complex<double> aleph,
        const db_matrix &AA,
        const db_matrix &B, const complex<double> beth, db_matrix &C);

protected:

};    // class bloc_matrix


// Non-member operators
inline db_matrix operator*(const complex<double> c, const db_matrix &d)
{ db_matrix a=d; a*=c; return a; }
inline db_matrix operator*(const db_matrix &d, const complex<double> c)
{ db_matrix a=d; a*=c; return a; }

inline db_matrix operator/(const complex<double> c, const db_matrix &d)
{ db_matrix a=d; a/=c; return a; }

inline  const db_matrix operator+(const db_matrix &o,const db_matrix &d) {
    db_matrix r = o;
    r += d;
    return r;
}
inline  const db_matrix operator-(const db_matrix &o,const db_matrix &d)  {
    db_matrix r = o;
    r -= d;
    return r;
}
inline const db_matrix operator*(const db_matrix &o,const db_matrix &d) {
    db_matrix r = o;
    r *= d;
    return r;
}

// Class used to be trhown as an exception

class db_matrix_error {
    string errMess;
public:

    db_matrix_error(){errMess="Error";}
    db_matrix_error(string m){errMess=m;}
    ~db_matrix_error(){}

    string getMess(){return errMess;}

};

#endif
