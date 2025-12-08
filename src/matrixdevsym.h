/***************************************************************************
*     CLASS matrixdevsym                                                   *
*     Davide Bucci, CROMA                                                  *
*     Jérôme Michallon, CROMA     2012                                     *
*     MINATEC-INPG, 3, parvis Luis Neel                                    *
*     38016, Grenoble CEDEX, France                                        *
*                                                                          *
*                                                                          *
***************************************************************************/

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

// Prevent multiple includes
#ifndef MATRIXDEVSYM_H
#define MATRIXDEVSYM_H
#include <iostream>


#include "block_matrix.h"
#include "matrixDev.h"

class fNonDevSymCart : public matrixDev {

    double alpha;

public:
    fNonDevSymCart(void) {alpha = -1.0;}
    bool setAngles(double ksinthetax, double ksinthetay)
    {
        return true;
    }

    void setAlpha(double a) {alpha = a;}
    double getAlpha() {return alpha;}

    db_matrix & createX1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR);

    db_matrix &createX2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X2,
                    db_matrix &M_fft, db_matrix &N_fft, db_matrix &S,
                    db_matrix &F, db_matrix &TTR);

    db_matrix &createX3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X3,
                    db_matrix &O_fft, db_matrix &M_fft, db_matrix &S,
                    db_matrix &G, db_matrix &TTR);

    db_matrix &createX4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR);

    db_matrix & createY1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO);

    virtual db_matrix &createY2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y2,
                    db_matrix &M_fft, db_matrix &Q_fft, db_matrix &Qm1_fft,
                    db_matrix &S,
                    db_matrix &F, db_matrix &TTO);

    db_matrix &createY3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y3,
                    db_matrix &M_fft, db_matrix &P_fft, db_matrix &Pm1_fft,
                    db_matrix &S,
                    db_matrix &G, db_matrix &TTO);

    db_matrix &createY4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO);
};

class fNonDevSym : public matrixDev {

    double alpha;

public:
    bool setAngles(double ksinthetax, double ksinthetay)
    {
        return true;
    }


    fNonDevSym(void) {alpha = -1.0;}

    void setAlpha(double a) {alpha = a;}
    double getAlpha() {return alpha;}

    db_matrix & createX1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR);

    db_matrix &createX2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X2,
                    db_matrix &M_fft, db_matrix &N_toep, db_matrix &S,
                    db_matrix &F, db_matrix &TTR);

    db_matrix &createX3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X3,
                    db_matrix &O_fft, db_matrix &M_toep, db_matrix &S,
                    db_matrix &G, db_matrix &TTR);

    db_matrix &createX4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR);

    db_matrix & createY1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO);

    virtual db_matrix &createY2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y2,
                    db_matrix &M_fft, db_matrix &Q_toep, db_matrix &Qm1_toep,
                    db_matrix &S,
                    db_matrix &F, db_matrix &TTO);

    db_matrix &createY3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y3,
                    db_matrix &M_fft, db_matrix &P_toep, db_matrix &Pm1_toep,
                    db_matrix &S,
                    db_matrix &G, db_matrix &TTO);

    db_matrix &createY4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO);
};
/** Implementation of the normal field vector. See:
    "Normal vector method for convergence improvement using the RCWA for
    crossed gratings", Thomas Schuster, Johannes Ruoff, Norbert Kerwien,
    Stephan Rafler, and Wolfgang Osten  JOSA A, Vol. 24, Issue 9, pp. 2880-2890
    (2007)
    http://dx.doi.org/10.1364/JOSAA.24.002880

*/
class fNFsym : public matrixDev {


public:

    bool setAngles(double ksinthetax, double ksinthetay)
    {
        return true;
    }
    // Some useful matrices which will contain the difference between the
    // Toeplitz matrix built by epsilon and the inverse of the Toeplitz
    // matrix built from 1/epsilon.
    db_matrix Delta_z;
    db_matrix Delta_r;


    // Those matrices will contain the vector field describing the normal
    // to the surfaces.

    db_matrix Nxx;
    db_matrix Nxy;
    db_matrix Nyy;

    db_matrix Nxx_fft;
    db_matrix Nxy_fft;
    db_matrix Nyy_fft;

    // input matrices are also stored:
    db_matrix Nxx_input;
    db_matrix Nxy_input;
    db_matrix Nyy_input;

    db_matrix & createX1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR);

    db_matrix &createX2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X2,
                    db_matrix &M_fft, db_matrix &N_toep, db_matrix &S,
                    db_matrix &F, db_matrix &TTR);

    db_matrix &createX3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X3,
                    db_matrix &O_fft, db_matrix &M_toep, db_matrix &S,
                    db_matrix &G, db_matrix &TTR);

    db_matrix &createX4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR);

    db_matrix & createY1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO);

    virtual db_matrix &createY2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y2,
                    db_matrix &M_fft, db_matrix &Q_toep, db_matrix &Qm1_toep,
                    db_matrix &S,
                    db_matrix &F, db_matrix &TTO);

    db_matrix &createY3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y3,
                    db_matrix &M_fft, db_matrix &P_toep, db_matrix &Pm1_toep,
                    db_matrix &S,
                    db_matrix &G, db_matrix &TTO);

    db_matrix &createY4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO);
};
#endif
