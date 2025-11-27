#include "block_matrix.h"
/** Toeplitz matrix generation and handling.

*/

// Truncate higher orders of a Toeplitz matrix. Use only when circulant is true
#define forcezero false

// Define this as true to use circulant matrix instead of Toeplitz's
#define circulant false



/** Modified bloc-Toeplitz matrix calculation

*/
db_matrix db_matrix::toeplitz_ones(double nux,
    double nuy,
    int p1x,
    int p1y,
    int p2x,
    int p2y)
{

    // Matrix size calculation.
    // The variables dimx and dimy are the number of harmonics taken in x
    // and in y, by taking into account only zero and positive frequencies.

    int size_x=ncol;
    int size_y=nrow;
    int dimx = int(floor(double(size_x)/2.0) + 1.0);
    int dimy = int(floor(double(size_y)/2.0) + 1.0);
    int dimension = dimx * dimy;
    int shiftx = size_x;
    int shifty = size_y;

    // index

    int i;
    int j;
    int p;
    int q;

    int indexa;
    int indexb;
    int indexa1;
    int indexb1;
    int indexc;
    int indexd;

    double mod;

    db_matrix R(dimension, dimension);
    bool tozero=false;

    for(i=0; i<dimy; ++i) {
        for(j=0; j<dimy; ++j) {
            for(p=0; p<dimx; ++p) {
                for(q=0; q<dimx; ++q) {
                    indexa=i-j;
                    indexb=p-q;
                    indexc=j-dimy/2;
                    indexd=q-dimx/2;

                    if(circulant) {
                        tozero=false;

                        if (indexa>dimy/2) {
                            indexa = indexa - dimy;
                            tozero = true;
                        }

                        if (indexb>dimx/2) {
                            indexb = indexb - dimx;
                            tozero = true;
                        }
                        if (indexa<-dimy/2) {
                            indexa = indexa + dimy;
                            tozero = true;
                        }

                        if (indexb<-dimx/2) {
                            indexb = indexb + dimx;
                            tozero = true;
                        }

                        if(tozero && forcezero) {
                            R(i*dimx+p,j*dimx+q)=0;
                            continue;
                        }
                    }

                    if(indexa<0)
                        indexa1=shifty+indexa;
                    else
                        indexa1=indexa;

                    if(indexb<0)
                        indexb1=shiftx+indexb;
                    else
                        indexb1=indexb;

                    if(indexc > dimy/2)
                        indexc = -(dimy-indexc);
                    if(indexd > dimx/2)
                        indexd = -(dimx-indexd);

                    mod = 1;

                    if (p1x > 0)
                        mod *= nux * indexb;
                    if (p1x > 1)
                        mod *= nux * indexb;
                    if (p1y > 0)
                        mod *= nuy * indexa;
                    if (p1y > 1)
                        mod *= nuy * indexa;

                    if (p2x > 0)
                        mod *= nux * indexd;
                    if (p2x > 1)
                        mod *= nux * indexd;
                    if (p2y > 0)
                        mod *= nuy * indexc;
                    if (p2y > 1)
                        mod *= nuy * indexc;

                    R(i*dimx+p,j*dimx+q)=mod;
                }
            }
        }
    }

    return R;
}

/** Modified bloc-Toeplitz matrix calculation for derivatives.
    NOTE: a multiplication times the imaginary unit is NOT done here.
    so this is a real matrix.
*/
db_matrix db_matrix::toeplitz_deriv(double nux,
    double nuy,
    int p1x,
    int p1y,
    double ksinthetax,
    double ksinthetay)
{

    // Matrix size calculation.
    // The variables dimx and dimy are the number of harmonics taken in x
    // and in y, by taking into account only zero and positive frequencies.

    int size_x=ncol;
    int size_y=nrow;
    int dimx = int(floor(double(size_x)/2.0) + 1.0);
    int dimy = int(floor(double(size_y)/2.0) + 1.0);
    int dimension = dimx * dimy;
    int shiftx = dimx/2;
    int shifty = dimy/2;

    // index

    int i;
    int j;
    //int p;
    //int q;

    db_matrix R(dimension, dimension);

    if(p1y>0) { // Y derivative
        for(i=0; i<dimy; ++i) {
            for(j=0;j<dimx;++j) {
                R(i*dimx+j,i*dimx+j)=nuy*p1y*(i-shifty)+ksinthetay;
            }
        }
    } else {    // X derivative
        for(i=0;i<dimy;++i) {
            for(j=0;j<dimx;++j) {
                    R(i*dimx+j,i*dimx+j)=nux*p1x*(j-shiftx)+ksinthetax;
                }
        }
    }

    return R;
}

/* Modified bloc-Toeplitz matrix calculation
   Useful when harmonics are stored in the FFT way
*/
db_matrix db_matrix::toeplitz_mod(double nux,
    double nuy,
    int p1x,
    int p1y,
    int p2x,
    int p2y)
{

    // Matrix size calculation.
    // The variables dimx and dimy are the number of harmonics taken in x
    // and in y, by taking into account only zero and positive frequencies.

    int size_x=ncol;
    int size_y=nrow;
    int dimx = int(floor(double(size_x)/2.0) + 1.0);
    int dimy = int(floor(double(size_y)/2.0) + 1.0);
    int dimension = dimx * dimy;
    int shiftx = size_x;
    int shifty = size_y;

    // index

    int i;
    int j;
    int p;
    int q;

    int indexa;
    int indexb;
    int indexa1;
    int indexb1;
    int indexc;
    int indexd;

    double mod;

    db_matrix R(dimension, dimension);
    bool tozero = false;

    for(i=0; i<dimy; ++i) {
        for(j=0; j<dimy; ++j) {
            for(p=0; p<dimx; ++p) {
                for(q=0; q<dimx; ++q) {
                    indexa=i-j;
                    indexb=p-q;
                    indexc=j-dimy/2;
                    indexd=q-dimx/2;

                    /* December 13 2008 correction: harmonics are not stored
                       in the same order as does the usual FFT if we are using
                       usual Toeplitz matrices!

                    indexc=j;
                    indexd=q;
                    */

                    if(circulant) {
                        tozero=false;

                        if (indexa>dimy/2) {
                            indexa = indexa - dimy;
                            tozero = true;
                        }

                        if (indexb>dimx/2) {
                            indexb = indexb - dimx;
                            tozero = true;
                        }
                        if (indexa<-dimy/2) {
                            indexa = indexa + dimy;
                            tozero = true;
                        }

                        if (indexb<-dimx/2) {
                            indexb = indexb + dimx;
                            tozero = true;
                        }

                        if(tozero && forcezero) {
                            R(i*dimx+p,j*dimx+q)=0;
                            continue;
                        }
                    }


                    if(indexa<0)
                        indexa1=shifty+indexa;
                    else
                        indexa1=indexa;

                    if(indexb<0)
                        indexb1=shiftx+indexb;
                    else
                        indexb1=indexb;

                    if(indexc > dimy/2)
                        indexc = -(dimy-indexc);
                    if(indexd > dimx/2)
                        indexd = -(dimx-indexd);

                    /* December 13 correction:
                    if(indexc > dimy/2)
                        indexc = -(dimy-indexc);
                    if(indexd > dimx/2)
                        indexd = -(dimx-indexd);
                    */
                    mod = 1;

                    if (p1x > 0)
                        mod *= nux * indexb;
                    if (p1x > 1)
                        mod *= nux * indexb;
                    if (p1y > 0)
                        mod *= nuy * indexa;
                    if (p1y > 1)
                        mod *= nuy * indexa;

                    if (p2x > 0)
                        mod *= nux * indexd;
                    if (p2x > 1)
                        mod *= nux * indexd;
                    if (p2y > 0)
                        mod *= nuy * indexc;
                    if (p2y > 1)
                        mod *= nuy * indexc;

                    R(i*dimx+p,j*dimx+q)=operator()(indexa1,indexb1)*mod;
                }
            }
        }
    }

    return R;
}

/** Create a matrix in 1D case corresponding to:
   sum over cx B(i+alpha*cx) C(cx)
   If symmetry is true, the sum is calculated from 0 to (Sx-1)
   Otherwise, from -(Sx-1) to (Sx-1)

   If alpha is -1, the function will return a Toeplitz matrix
   Otherwise, it is not really a Toeplitz matrix (NOTE: and what the is then?)

   columnVector select a vector in the matrix
   elementSize refert to the size of each element in the matrix Vector
*/

inline int ABS(int x)
{
    return x>0?x:-x;
}
db_matrix db_matrix::toeplitz1D(int elementsize, int columnVector,
    bool symmetry, int alpha)
{

    // Matrix size calculation.
    // The variable dim is the number of harmonics
    // by taking into account only zero and positive frequencies.
    // shifty is the index of the zero frequency

    int size =nrow;
    int dim;
    int nbelements = size/elementsize;
    int shifty =  int(floor(double(nbelements)/2.0) + 1.0) - 1;

    if (symmetry)// the size of the matrix is reduced in the case of symmetry
        dim = int(floor(double(nbelements)/4.0) + 1.0);
    else
        dim = int(floor(double(nbelements)/2.0) + 1.0);

    // if we have alpha = +1 and no symmetry, shift is used to get the correct
    // matrix
    // this is used while calculating the poynting
    int shift = 0;
    if (!symmetry && alpha == +1)
        shift = int(floor(double(nbelements)/4.0) + 1.0) - 1;;

    // index
    int i;
    int j;

    // for copying the element
    int p;
    int q;


    db_matrix R(dim*elementsize, dim*elementsize);

    for(j=0; j<dim; ++j) {
        for(i=0;i<dim; ++i) {

            // sum over cx B(i+alpha*cx) C(cx)
            for (p = 0; p < elementsize; ++p) {
                for (q = 0; q < elementsize; ++q) {

                    // check if the index are out of the matrix, this can append
                    // when alpha = +1 without symmetries
                    if (((i-shift) + alpha*(j-shift) )+shifty < nbelements &&
                            ((i-shift) + alpha*(j-shift) )+shifty >=0 ) {

                        // copy the element from V into R
                        if (symmetry)
                            R(i*elementsize + p,j*elementsize + q) =
                                operator()((ABS((i-shift)+alpha*(j-shift)) +
                                    shifty) * elementsize  + p,
                                            columnVector*elementsize + q);
                        else
                            R(i*elementsize + p,j*elementsize + q) =
                                operator()((((i-shift)+alpha*(j-shift)) +
                                    shifty) * elementsize  + p,
                                            columnVector*elementsize + q);
                    }
                }
            }

        }
    }
    return R;
}



/** Computation of the sum of the 1D toeplitz matrices in case of symmetry
    Syntax: Toeplitz1D.toeplitz1D_sum(Toeplitz1D+,gamma,IsOdd);

   if B is even (IsOdd false)
        if alpha = 1
            Ci = A0*Bi + sum (B|i-ax| + Bi+ax)Aax

        if alpha = -1
            Ci = A0*Bi + sum (B|i-ax| - Bi+ax)Aax

   if B is odd (IsOdd true)
        if alpha = 1
            Ci = A0*Bi + sum (Bi-ax + Bi+ax)Aax     for i-ax > 0
            Ci = A0*Bi + sum (-B|i-ax| + Bi+ax)Aax  for i-ax < 0

        if alpha = -1
            Ci = A0*Bi + sum (Bi-ax - Bi+ax)Aax     for i-ax > 0
            Ci = A0*Bi + sum (-B|i-ax| - Bi+ax)Aax      for i-ax < 0

   elementSize refers to the size of each element of the Toeplitz1D
*/
db_matrix &db_matrix::toeplitz1D_sum(db_matrix &Toeplitz1Dp, int elementsize,
                                int gamma, bool IsOdd)
{
    int dimy = nrow/elementsize;
    int dimx = ncol/elementsize;
    db_matrix element;


    // If B is odd, we need to remove twice the value B|i-ax|
    //  so we end up with -B|i-ax| .
    if (IsOdd) {
        for(int i=0; i<dimy; ++i) {
            for(int j=i+1;j<dimx;++j) {
                add(*this,-2.0,i*elementsize,j*elementsize,i*elementsize,
                        j*elementsize,elementsize, elementsize);
            }
        }

    }

    // for every element, we compute the sum of the toeplitz matrix
    // except for the first colomn
    for(int i=0; i<dimy; ++i) {
        for(int j=1;j<dimx;++j) {
            add(Toeplitz1Dp,gamma,i*elementsize,j*elementsize,i*elementsize,
                    j*elementsize,elementsize, elementsize);
        }
    }

    // if A is odd, the first colomn must be null except for the constant
    // term B0
    if (gamma == -1) {
        // remove the first row:
        for(int i=1; i<dimy; ++i) {
            add(*this,-1.0,i*elementsize,0,
                    i*elementsize,0,elementsize, elementsize);
        }
    }

    /*int dimy = nrow/elementsize;
    int dimx = ncol/elementsize;
    db_matrix element;

    int delta;

    // First colomn
    if (gamma == -1 || IsOdd) {
        // When Aax and Bbx are not even, the first colomn is empty except
        // for the first row and the first colomn
        for(int i=1;i<dimy;++i) {
            add(*this,-1,i*elementsize,0,i*elementsize,0,
                            elementsize, elementsize);
        }
    }
    // first row
    if ((gamma == 1 && !IsOdd)  || (gamma == -1 && IsOdd)){
        delta = 1;
    } else {
        delta =-1;
    }
    for(int j=1;j<dimx;++j) {
        add(Toeplitz1Dp,delta,0,j*elementsize,0,
                j*elementsize,elementsize, elementsize);
    }

    // everything except first row and colomn
    for(int i=1; i<dimy; ++i) {
        for(int j=1;j<dimx;++j) {
            add(Toeplitz1Dp,gamma,i*elementsize,j*elementsize,i*elementsize,
                    j*elementsize,elementsize, elementsize);
        }
    }
*/
    return *this;

}

/** Create a matrix in 2D case:
   in case of symmetry, the sum is calculated from 0 to (Sx-1)/2
   Otherwise, from -(Sx-1)/2 to (Sx-1)/2

   if B is even (Oddx and Oddy false)
        if alpha = 1
            Ci = A0*Bi + sum (B|i-ax| + Bi+ax)Aax

        if alpha = -1
            Ci = A0*Bi + sum (B|i-ax| - Bi+ax)Aax

   if B is odd (Oddx and Oddy)
        if alpha = 1
            Ci = A0*Bi + sum (Bi-ax + Bi+ax)Aax     for i-ax > 0
            Ci = A0*Bi + sum (-Bi-ax + Bi+ax)Aax    for i-ax < 0

        if alpha = -1
            Ci = A0*Bi + sum (Bi-ax - Bi+ax)Aax     for i-ax > 0
            Ci = A0*Bi + sum (-Bi-ax - Bi+ax)Aax    for i-ax < 0
*/
db_matrix db_matrix::toeplitz_sym(bool xsymmetric, bool ysymmetric,
                            int alpha, int beta, bool IsOddx, bool IsOddy)
{

    // Matrix size calculation.
    // The variables dimx and dimy are the number of harmonics taken in x
    // and in y, by taking into account only zero and positive frequencies.

    int size_x=ncol;
    int size_y=nrow;
    int dimx = int(floor(double(size_x)/2.0));
    int dimy = int(floor(double(size_y)/2.0));

    // the size of the matrix is reduced in the case of symmetry
    if (xsymmetric)
        dimx = int(floor(double(size_x)/4.0) + 1.0);
    else
        dimx = int(floor(double(size_x)/2.0) + 1.0);

    // the size of the matrix is reduced in the case of symmetry
    if (ysymmetric)
        dimy = int(floor(double(size_y)/4.0) + 1.0);
    else
        dimy = int(floor(double(size_y)/2.0) + 1.0);


    // index

    int i;

    // sort the coefficient comming from the FFT
    // This matrix consists of vertical vector containing the x
    // coefficient for a given y coefficent (column)
    db_matrix FFT_ordered = fftshift().transpose();

    // creation of a vector containing 1D toeplitz matrices:
    db_matrix Vector(size_y*dimx, dimx);
    db_matrix Toeplitz1Dm(dimx, dimx);
    db_matrix Toeplitz1Dp(dimx, dimx);

    for(i=0;i<size_y; ++i) {

        Toeplitz1Dm = FFT_ordered.toeplitz1D(1, i, xsymmetric, -1);
        if (xsymmetric == true) {
            Toeplitz1Dp = FFT_ordered.toeplitz1D(1, i, xsymmetric, +1);
            Toeplitz1Dm.toeplitz1D_sum(Toeplitz1Dp, 1,alpha,IsOddx);
        }
        /*
        cout << "Toeplitz1D i  = " << i << endl;
        Toeplitz1Dm.printMatrix();*/
        Vector.copy(Toeplitz1Dm, i*dimx,0);
    }

    // free the space that will not be used anymore
    Toeplitz1Dm.kill();
    Toeplitz1Dp.kill();

    // creation of the toeplitz matrix using a vector of 1D toeplitz matrices
    // (for the y components)
    db_matrix Toeplitzm(dimx*dimy, dimx*dimy);
    db_matrix Toeplitzp(dimx*dimy, dimx*dimy);

    Toeplitzm = Vector.toeplitz1D(dimx, 0, ysymmetric, -1);
    if (ysymmetric == true) {
        Toeplitzp = Vector.toeplitz1D(dimx, 0, ysymmetric, +1);
        Toeplitzm.toeplitz1D_sum(Toeplitzp, dimx,beta,IsOddy);
    }

    return Toeplitzm;
}

/** Computation of derivative matrices in case of symmetry
*/
db_matrix db_matrix::toeplitz_deriv_sym(bool xsymmetric, bool ysymmetric,
             double nux, double nuy, int p1x, int p1y)
{

    // Matrix size calculation.
    // The variables dimx and dimy are the number of harmonics taken in x
    // and in y, by taking into account only zero and positive frequencies.

    int size_x=ncol;
    int size_y=nrow;
    int dimx;
    int dimy;

    int shiftx;
    int shifty;

    // the size of the matrix is reduced in the case of symmetry
    if (xsymmetric) {
        dimx = int(floor(double(size_x)/4.0) + 1.0);
        shiftx = 0;
    } else {
        dimx = int(floor(double(size_x)/2.0) + 1.0);
        shiftx = dimx/2 ;
    }
    // the size of the matrix is reduced in the case of symmetry
    if (ysymmetric) {
        dimy = int(floor(double(size_y)/4.0) + 1.0);
        shifty = 0;
    } else {
        dimy = int(floor(double(size_y)/2.0) + 1.0);
        shifty = dimy/2;
    }

    int dimension = dimx * dimy;

    // index

    int i;
    int j;


    db_matrix R(dimension, dimension);

    if(p1y>0) { // Y derivative
        for(i=0; i<dimy; ++i) {
            for(j=0;j<dimx;++j) {
                R(i*dimx+j,i*dimx+j)=nuy*p1y*(i-shifty);
            }
        }
    } else {
        for(i=0;i<dimy;++i) {
            for(j=0;j<dimx;++j) {
                R(i*dimx+j,i*dimx+j)=nux*p1x*(j-shiftx);
            }
        }
    }

    return R;
}
