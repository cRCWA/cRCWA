/***************************************************************************
*     CLASS matrixdevsym                                                   *
*     Davide Bucci, CROMA     2008-present                                 *
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

#include <iostream>
#include "matrixdevsym.h"
#include "block_matrix.h"
/** Implementation of:
    matdev nds
    matdev las alpha

*/
db_matrix & fNonDevSymCart::createX1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{
    X1 = Kr;
    X1 *= TTR;
    X1 *= Kz;
    X1 *= -1.0;

    X1 *=(1.0/omega);

    //cout<<"X1 CARTESIEN... ";
    //cout.flush();
    return X1;
}

db_matrix &fNonDevSymCart::createX2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X2,
                    db_matrix &M_fft, db_matrix &N_toep, db_matrix &S,
                    db_matrix &F, db_matrix &TTR)
{

    X2 =Kr;
    X2 *= TTR;
    X2 *= Kr;

    db_matrix toto = N_toep;
    toto *= (-omega*omega);
    X2 += toto;

    X2 *= (1.0/omega);

    return X2;
}

db_matrix &fNonDevSymCart::createX3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X3,
                    db_matrix &O_fft, db_matrix &M_toep, db_matrix &S,
                    db_matrix &G, db_matrix &TTR)
{

    X3 = Kz;
    X3*=TTR;
    X3*= Kz;

    X3 *=-1.0;

    db_matrix toto = M_toep;
    toto *= (omega*omega);
    X3 += toto;


    X3 *=(1.0/omega);

    return X3;
}
db_matrix &fNonDevSymCart::createX4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{
    X4 = Kz;
    X4 *=TTR;
    X4 *= Kr;

    X4 *= (1.0/omega);

    return X4;
}
db_matrix & fNonDevSymCart::createY1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{

    Y1 =  Kr;
    Y1 *= TTO;
    Y1 *= Kz;

    Y1 *=(1.0/omega);

    return Y1;

}
db_matrix &fNonDevSymCart::createY2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y2,
                    db_matrix &M_fft, db_matrix &Q_toep, db_matrix &Qm1_toep,
                    db_matrix &S,
                    db_matrix &F, db_matrix &TTO)
{
    Y2 = Kr;
    Y2 *=TTO;
    Y2 *= Kr;
    Y2 *=-1.0;

    db_matrix toto = Q_toep;

    if(alpha>=0) {
        toto *= alpha;
        db_matrix toto1 = Qm1_toep;
        toto1.invert();
        toto1 *= (1.0-alpha);

        toto += toto1;
        toto1.kill();
    }

    toto *= (omega*omega);
    Y2+=toto;

    Y2 *= (1.0/omega);
    return Y2;
}

db_matrix &fNonDevSymCart::createY3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y3,
                    db_matrix &M_fft, db_matrix &P_toep, db_matrix &Pm1_toep,
                    db_matrix &S,
                    db_matrix &G, db_matrix &TTO)
{

    Y3 = Kz;
    Y3 *= TTO;
    Y3 *=Kz;

    db_matrix toto = P_toep;

    if(alpha>=0) {
        toto *= (1.0-alpha);
        db_matrix toto1 = Pm1_toep;
        toto1.invert();
        toto1 *= alpha;

        toto += toto1;
        toto1.kill();
    }
    toto *= (-omega*omega);
    Y3 += toto;

    Y3 *=(1.0/omega);

    return Y3;
}
db_matrix &fNonDevSymCart::createY4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{
    Y4 = Kr;
    Y4 *= TTO;
    Y4 *= Kz;
    Y4*=-1.0;

    Y4 *= (1.0/omega);
    return Y4;
}
/****************************************************************************/
db_matrix & fNonDevSym::createX1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{
    X1 = F;
    X1 *= Kr;
    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        X1 *= S;
    }
    X1 *= TTR;
    X1 *= G;
    X1 *= Kz;
    X1 *= -1.0;

    X1 *=(1.0/omega);

    return X1;
}

db_matrix &fNonDevSym::createX2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X2,
                    db_matrix &M_fft, db_matrix &N_toep, db_matrix &S,
                    db_matrix &F, db_matrix &TTR)
{

    X2 = F;
    X2 *=Kr;

    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        X2*= S;
    }

    X2 *= TTR;
    X2 *= F;
    X2 *= Kr;

    db_matrix toto = N_toep;
    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        toto *= S;
    }
    toto *= (-omega*omega);
    X2 += toto;

    X2 *= (1.0/omega);

    return X2;
}

db_matrix &fNonDevSym::createX3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X3,
                    db_matrix &O_fft, db_matrix &M_toep, db_matrix &S,
                    db_matrix &G, db_matrix &TTR)
{
    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        X3 = S;
        X3 *= G;
    } else {
        X3 = G;
    }
    X3 *= Kz;
    X3*=TTR;
    X3 *=G;
    X3*= Kz;

    X3 *=-1.0;

    db_matrix toto = M_toep;
    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        toto *= S;
    }
    toto *= (omega*omega);
    X3 += toto;


    X3 *=(1.0/omega);

    return X3;
}
db_matrix &fNonDevSym::createX4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{
    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        X4 = S;
        X4 *= G;
    } else {
        X4 = G;
    }
    X4 *= Kz;

    X4 *=TTR;

    X4 *=F;
    X4 *= Kr;

    X4 *= (1.0/omega);

    return X4;
}
db_matrix & fNonDevSym::createY1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{

    Y1 = F;
    Y1 *=  Kr;
    if(isBent){
        Y1 *= S;
    }
    Y1 *= TTO;
    Y1 *= G;
    Y1 *= Kz;

    Y1 *=(1.0/omega);

    return Y1;

}
db_matrix &fNonDevSym::createY2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y2,
                    db_matrix &M_fft, db_matrix &Q_toep, db_matrix &Qm1_toep,
                    db_matrix &S,
                    db_matrix &F, db_matrix &TTO)
{
    Y2 = F;
    Y2 *= Kr;
    if(isBent){
        Y2 *= S;
    }
    Y2 *=TTO;
    Y2 *=F;
    Y2 *= Kr;
    Y2 *=-1.0;

    db_matrix toto = Q_toep;

    //cout<<"alpha="<<alpha<<"\n";
    if(alpha>=0) {
        toto *= alpha;
        db_matrix toto1 = Qm1_toep;
        toto1.invert();
        toto1 *= (1.0-alpha);

        toto += toto1;
        toto1.kill();
    }

    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        toto *= S;
    }
    toto *= (omega*omega);
    Y2+=toto;

    Y2 *= (1.0/omega);

    return Y2;
}

db_matrix &fNonDevSym::createY3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y3,
                    db_matrix &M_fft, db_matrix &P_toep, db_matrix &Pm1_toep,
                    db_matrix &S,
                    db_matrix &G, db_matrix &TTO)
{
    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        Y3 = S;
        Y3 *= G;
    } else {
        Y3 = G;
    }


    Y3 *= Kz;
    Y3 *= TTO;

    Y3 *= G;
    Y3 *=Kz;

    db_matrix toto = P_toep;

    if(alpha>=0) {
        toto *= (1.0-alpha);
        db_matrix toto1 = Pm1_toep;
        toto1.invert();
        toto1 *= alpha;

        toto += toto1;
        toto1.kill();
    }

    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        toto *= S;
    }
    toto *= (-omega*omega);
    Y3 += toto;

    Y3 *=(1.0/omega);

    return Y3;
}
db_matrix &fNonDevSym::createY4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{
    /*if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        Y4 = S;
        Y4 *= F;
    } else {
        Y4 = F;
    }
    Y4 *= Kr;
    Y4 *= TTO;
    Y4 *= G;
    Y4 *= Kz;
    Y4*=-1.0;*/

//JEROME
    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        Y4 = S;
        Y4 *= G;
    } else {
        Y4 = G;
    }
    Y4 *= Kz;
    Y4 *= TTO;
    Y4 *= F;
    Y4 *= Kr;
    Y4 *=-1.0;

    Y4 *= (1.0/omega);
    return Y4;
}
/** matdev nfs normal_vector_x.f3d normal_vector_y.f3d

    Implementation of the normal field vector. See:
    "Normal vector method for convergence improvement using the RCWA for
    crossed gratings", Thomas Schuster, Johannes Ruoff, Norbert Kerwien,
    Stephan Rafler, and Wolfgang Osten  JOSA A, Vol. 24, Issue 9, pp. 2880-2890
    (2007)
    http://dx.doi.org/10.1364/JOSAA.24.002880

    Order on which the matrices are created:

    X2, X4, X1, X3, Y2, Y1, Y3, Y4.

    In Y2, we create Delta_x used also by Y1 (and then destroyed)
    In Y3, we create Delta_y, then used and destroyed by Y4

*/

db_matrix & fNFsym::createX1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{
    X1 = F;
    X1 *= Kr;
    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        X1 *= S;
    }
    X1 *= TTR;
    X1 *= G;
    X1 *= Kz;
    X1 *= -1.0;

    X1 *=(1.0/omega);
    return X1;
}

db_matrix &fNFsym::createX2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X2,
                    db_matrix &M_fft, db_matrix &N_toep, db_matrix &S,
                    db_matrix &F, db_matrix &TTR)
{

    X2 = F;
    X2 *= Kr;

    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        X2*= S;
    }
    X2 *= TTR;
    X2 *= F;
    X2 *= Kr;

    db_matrix toto = N_toep;
    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        toto *= S;
    }
    toto *= (-omega*omega);
    X2 += toto;

    X2 *= (1.0/omega);

    return X2;
}

db_matrix &fNFsym::createX3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X3,
                    db_matrix &O_fft, db_matrix &M_toep, db_matrix &S,
                    db_matrix &G, db_matrix &TTR)
{
    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        X3 = S;
        X3 *= G;
    } else {
        X3 = G;
    }
    X3 *= Kz;
    X3*=TTR;
    X3 *=G;
    X3*= Kz;

    X3 *=-1.0;

    db_matrix toto = M_toep;
    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        toto *= S;
    }
    toto *= (omega*omega);
    X3 += toto;


    X3 *=(1.0/omega);

    return X3;
}
db_matrix &fNFsym::createX4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{
    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        X4 = S;
        X4 *= G;
    } else {
        X4 = G;
    }
    X4 *= Kz;

    X4 *=TTR;

    X4 *=F;
    X4 *= Kr;

    X4 *= (1.0/omega);

    return X4;
}
db_matrix & fNFsym::createY1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{

    Y1 = F;
    Y1 *=  Kr;
    if(isBent){
        Y1 *= S;
    }
    Y1 *= TTO;
    Y1 *= G;
    Y1 *= Kz;

    db_matrix toto = Delta_z;
/*  toto*=Nr3;
    toto*=Nz1;*/

    toto *= Nxy;

    if(isBent){
        toto*=S;
    }
    toto*= (-omega*omega);

    Y1 += toto;

    Y1 *=(1.0/omega);

//JEROME
    Delta_z.kill();
    return Y1;

}
db_matrix &fNFsym::createY2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y2,
                    db_matrix &M_fft, db_matrix &Q_toep, db_matrix &Qm1_toep,
                    db_matrix &S,
                    db_matrix &F, db_matrix &TTO)
{
    Delta_z = Q_toep - Qm1_toep.invert();



    Y2 = F;
    Y2 *= Kr;
    if(isBent){
        Y2 *= S;
    }
    Y2 *=TTO;
    Y2 *=F;
    Y2 *= Kr;
    Y2 *=-1.0;

    db_matrix toto = Q_toep;

//  toto -= Delta_z*Nz3*Nz2;
    toto -= Delta_z*Nyy;

    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        toto *= S;
    }

    toto *= (omega*omega);
    Y2+=toto;

    Y2 *= (1.0/omega);
    return Y2;
}

db_matrix &fNFsym::createY3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y3,
                    db_matrix &M_fft, db_matrix &P_toep, db_matrix &Pm1_toep,
                    db_matrix &S,
                    db_matrix &G, db_matrix &TTO)
{
    Delta_r = P_toep - Pm1_toep.invert();

    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        Y3 = S;
        Y3 *= G;
    } else {
        Y3 = G;
    }


    Y3 *= Kz;
    Y3 *= TTO;

    Y3 *= G;
    Y3 *= Kz;

    db_matrix toto = P_toep;

//  toto -= Delta_r*Nr2*Nr1;
    toto -= Delta_r*Nxx;

    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        toto *= S;
    }
    toto *= (-omega*omega);
    Y3 += toto;

    Y3 *=(1.0/omega);

    return Y3;
}
db_matrix &fNFsym::createY4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{
    /*if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        Y4 = S;
        Y4 *= F;
    } else {
        Y4 = F;
    }
    Y4 *= Kr;
    Y4 *= TTO;
    Y4 *= G;
    Y4 *= Kz;
    Y4 *=-1.0;*/
//JEROME
    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        Y4 = S;
        Y4 *= G;
    } else {
        Y4 = G;
    }
    Y4 *= Kz;
    Y4 *= TTO;
    Y4 *= F;
    Y4 *= Kr;
    Y4 *=-1.0;

    db_matrix toto = Delta_r;

/*  toto *=Nr2;
    toto *=Nz2;*/
    toto *=Nxy;

    if(isBent) {
        if(S.isEmpty()) {
            cerr << "Warning: The S matrix is not supposed to be empty.";
        }
        toto *= S;
    }
    toto *= (omega*omega);

    Y4 += toto;

    Y4 *= (1.0/omega);
//JEROME
    Delta_r.kill();
    return Y4;
}
