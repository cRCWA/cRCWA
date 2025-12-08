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
#include "matrixDev.h"
#include "block_matrix.h"


/** Developments for matdev be

*/
db_matrix & fPMLbefore::createX1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{

    /*
    X1 = ((isBent?0:1.0)*(complex<double>(0,1))*
          G*db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,0,1))
          )-S*F*G*(db_matrix::hadamard(TTR,
          M_fft.toeplitz_ones(nux,nuy,1,0,0,1))
         +db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,1)));
    */


    X1 = TTR;
    X1.hadamard(M_fft.toeplitz_ones(nux,nuy,1,0,0,1));
    db_matrix pipo = TTR;
    pipo.hadamard(M_fft.toeplitz_ones(nux,nuy,0,0,1,1));
    X1+=pipo;

    if(isBent) {
        pipo = S;
        pipo *=F;
    } else {
        pipo = F;
    }
    pipo *= -1.0;

    pipo *=G;

    pipo *= X1;
    X1=pipo;
    pipo.kill();

    if(isBent) {
        pipo = (complex<double>(0,1))*
          db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,0,1));
        X1= multsum(G,pipo,X1);
        pipo.kill();
    }

    X1 *=(1.0/omega);

    cout<<"X1... ";
    cout.flush();
    return X1;
}

db_matrix &fPMLbefore::createX2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X2,
                    db_matrix &M_fft, db_matrix &N_fft, db_matrix &S,
                    db_matrix &F, db_matrix &TTR)
{
/*
    db_matrix pipo2=db_matrix::hadamard(TTR,
        M_fft.toeplitz_ones(nux,nuy,1,0,1,0))
       + db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,2,0));

    X2 = ((isBent?0:1.0)*(complex<double>(0,-1))*F*
         db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,0)))  +
        S*((complex<double>(-omega*omega,0) *
        N_fft.toeplitz_mod(nux, nuy, 0,0,0,0) +
        F*F*pipo2));
    pipo2.kill();
*/


    db_matrix pipo =complex<double>(-omega*omega,0) *
        N_fft.toeplitz_mod(nux, nuy, 0,0,0,0);

    X2 = db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,1,0,1,0));
    X2 += db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,2,0));

    pipo = multsum(F*F, X2, pipo);

    if(isBent) {
        X2 = S;
        X2 *=pipo;
    } else {
        X2 = pipo;
    }

    pipo.kill();

    if(isBent){
        X2 = multsum(F,(complex<double>(0,-1))*
            db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,0)),X2);
    }

    X2 *= (1.0/omega);

    cout<<"X2... ";
    cout.flush();
    return X2;
}

db_matrix &fPMLbefore::createX3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X3,
                    db_matrix &O_fft, db_matrix &M_fft, db_matrix &S,
                    db_matrix &G, db_matrix &TTR)
{
/*
    db_matrix pipo55=db_matrix::hadamard(TTR,
        M_fft.toeplitz_ones(nux,nuy,0,1,0,1))
        + db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,0,2));

    X3 = S*((omega*omega)*
        M_fft.toeplitz_mod(nux, nuy, 0,0,0,0) -
        G*G*pipo55);

    pipo55.kill();
*/

    db_matrix pipo = TTR;
    pipo.hadamard(M_fft.toeplitz_ones(nux,nuy,0,1,0,1));
    pipo+=db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,0,2));

    if(isBent) {
        X3 = S;
        X3 *= G;
    } else {
        X3 = G;
    }
    X3 *= -1.0;

    X3 *=G;
    X3 *=pipo;
    pipo.kill();

    if(isBent) {
        X3 = multsum(S, (omega*omega)*
            M_fft.toeplitz_mod(nux, nuy, 0,0,0,0), X3);
    } else {
        X3 += (omega*omega)*
            M_fft.toeplitz_mod(nux, nuy, 0,0,0,0);
    }

    X3 *=(1.0/omega);

    return X3;
}
db_matrix &fPMLbefore::createX4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{
/*
    X4 = S*F*G*(db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,1,1,0))+
        db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,1)));
*/
    if(isBent) {
        X4 = S;
        X4 *=F;
    } else {
        X4 = F;
    }

    X4 *=G;
    db_matrix pipo = TTR;
    pipo.hadamard(M_fft.toeplitz_ones(nux,nuy,0,1,1,0));
    pipo+=db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,1));
    X4 *=pipo;
    pipo.kill();

    X4 *= (1.0/omega);

    return X4;
}
db_matrix & fPMLbefore::createY1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{

    Y1 = (((isBent?0:1.0)*(complex<double>(0,-1)))*
        G*db_matrix::hadamard(TTO, M_fft.toeplitz_ones(nux,nuy,0,0,0,1))) +
        S*F*G*(db_matrix::hadamard(TTO, M_fft.toeplitz_ones(nux,nuy,1,0,0,1)) +
        db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,1)));
/*

    if(isBent) {
        Y1 = S;
        Y1 *= F;
    } else {
        Y1 = F;
    }
    Y1 *= G;

    db_matrix pipo =
        db_matrix::hadamard(TTO, M_fft.toeplitz_ones(nux,nuy,1,0,0,1));
    pipo += db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,1));
    Y1 *= pipo;
    pipo.kill();

    if(isBent) {
        Y1 = multsum(complex<double>(0,-1)*G, // ***************
            db_matrix::hadamard(TTO, M_fft.toeplitz_ones(nux,nuy,0,0,0,1)), Y1);
    }   */
    Y1 *=(1.0/omega);
    return Y1;

}
db_matrix &fPMLbefore::createY2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y2,
                    db_matrix &M_fft, db_matrix &Q_fft, db_matrix &Qm1_fft,
                    db_matrix &S,
                    db_matrix &F, db_matrix &TTO)
{
    /*
    db_matrix pipo5=(
        db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,1,0,1,0))
        + db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,2,0)));

    Y2 = (((isBent?0:1.0)*(complex<double>(0,1)))*
        F*db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,0)))+
        S*((omega*omega)* Q_fft.toeplitz_mod(nux, nuy, 0,0,0,0) -
          F*F*pipo5); // ************************* <----
    pipo5.kill();
    */


    Y2 = F;
    Y2 *=F;
    Y2 *=-1.0;
    db_matrix pipo =
        db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,1,0,1,0));
    pipo +=  db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,2,0));
    Y2 *=pipo;
    pipo.kill();
    pipo = Q_fft.toeplitz_mod(nux, nuy, 0,0,0,0);
    pipo *=(omega*omega);
    Y2 += pipo;



    if(isBent) {
        pipo = Y2;
        Y2 = S;
        Y2 *= pipo;
    } else {
        //Y2 = pipo;
    }

    pipo.kill();
    if(isBent) {
        Y2 = multsum(complex<double>(0,-1)*F, // ************
            db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,0)), Y2);
    }


    Y2 *= (1.0/omega);
    cout<<"Y2... ";
    cout.flush();
    return Y2;
}

// Kills TTO, S
db_matrix &fPMLbefore::createY3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y3,
                    db_matrix &M_fft, db_matrix &P_fft,db_matrix &Pm1_fft,
                    db_matrix &S,
                    db_matrix &G, db_matrix &TTO)
{

    Y3 = S*((-omega*omega)* P_fft.toeplitz_mod(nux, nuy, 0,0,0,0)+
        G*G*(db_matrix::hadamard(TTO, M_fft.toeplitz_ones(nux,nuy,0,1,0,1))+
             db_matrix::hadamard(TTO, M_fft.toeplitz_ones(nux,nuy,0,0,0,2))));
    /*
    Y3=TTO;
    Y3.hadamard(M_fft.toeplitz_ones(nux,nuy,0,1,0,1));
    Y3+=db_matrix::hadamard(TTO, M_fft.toeplitz_ones(nux,nuy,0,0,0,2));
    // DO NOT USE MATRIX TTO AS A TEMP!!!
    TTO.kill();
    TTO =G;
    TTO *=G;
    // G.kill(); // why this? It is used in Y4.
    db_matrix pipo=P_fft.toeplitz_mod(nux, nuy, 0,0,0,0);
    pipo*=(-omega*omega);
    pipo = multsum(TTO, Y3, pipo);

    if(isBent) {
        Y3 = S;
        //S.kill();
        Y3 *= pipo;

    } else {
        Y3 = pipo;
    }


    pipo.kill(); */

    Y3 *=(1.0/omega);

    return Y3;
}
db_matrix &fPMLbefore::createY4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{


    Y4 =-S*F*G*(db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,1,1,0)) +
        db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,1)));
/*
    db_matrix pipo =
        db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,1,1,0));
    pipo += db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,1));

    if(isBent) {
        Y4 = S;
        Y4 *= F;
    } else {
        Y4 = F;
    }
    // F is no longer needed.
    F.kill();

    Y4 *= G;
    Y4 *=-1.0;

    Y4 *=pipo;
    pipo.kill();
    */
    Y4 *= (1.0/omega);
    return Y4;
}



/*  Developments for matdev an

*/
db_matrix & fPMLafter::createX1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{




    X1 = (isBent?1.0:0)*(complex<double>(0,1))*
          db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,0,1))*G;

    if(isBent) {
        X1-=S*(db_matrix::hadamard(TTR,
            M_fft.toeplitz_ones(nux,nuy,1,0,0,1))
            +db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,1)))
            *F*G;
    } else {
        X1-=(db_matrix::hadamard(TTR,
            M_fft.toeplitz_ones(nux,nuy,1,0,0,1))
            +db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,1)))
            *F*G;

    }

    X1 *=(1.0/omega);

    return X1;
}

db_matrix &fPMLafter::createX2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X2,
                    db_matrix &M_fft, db_matrix &N_fft, db_matrix &S,
                    db_matrix &F, db_matrix &TTR)
{
    db_matrix pipo=db_matrix::hadamard(TTR,
        M_fft.toeplitz_ones(nux,nuy,1,0,1,0))
       + db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,2,0));


    X2 = (isBent?1.0:0)*(complex<double>(0,-1))*
         db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,0))*F;
    if(isBent) {
        X2+=S*((complex<double>(-omega*omega,0) *
            N_fft.toeplitz_mod(nux, nuy, 0,0,0,0) +
            pipo*(F*F)));
     } else {
        X2+=((complex<double>(-omega*omega,0) *
            N_fft.toeplitz_mod(nux, nuy, 0,0,0,0) +
            pipo*(F*F)));
     }

    pipo.kill();

    X2 *= (1.0/omega);

    return X2;
}

db_matrix &fPMLafter::createX3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X3,
                    db_matrix &O_fft, db_matrix &M_fft, db_matrix &S,
                    db_matrix &G, db_matrix &TTR)
{



    db_matrix pipo=db_matrix::hadamard(TTR,
        M_fft.toeplitz_ones(nux,nuy,0,1,0,1))
        + db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,0,2));
    // PB SIGNES???
    if(isBent) {
        X3 = S*((omega*omega)*
            M_fft.toeplitz_mod(nux, nuy, 0,0,0,0) -
            pipo*(G*G));
    } else {
        X3 = ((omega*omega)*
            M_fft.toeplitz_mod(nux, nuy, 0,0,0,0) -
            pipo*(G*G));
    }

    pipo.kill();

    X3 *=(1.0/omega);

    return X3;
}
db_matrix &fPMLafter::createX4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{

    if(isBent) {
        X4 = S*(db_matrix::hadamard(TTR,
            M_fft.toeplitz_ones(nux,nuy,0,1,1,0))+
            db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,1)))
            *(F*G);
    } else {
        X4 =(db_matrix::hadamard(TTR,
            M_fft.toeplitz_ones(nux,nuy,0,1,1,0))+
            db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,1)))
            *(F*G);
    }

    X4 *= (1.0/omega);

    return X4;
}
db_matrix & fPMLafter::createY1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{

    Y1 = (isBent?1.0:0.0)*(complex<double>(0,-1))*
        db_matrix::hadamard(TTO, M_fft.toeplitz_ones(nux,nuy,0,0,0,1))*G;

    if(isBent) {
        Y1+=S*(db_matrix::hadamard(TTO,
            M_fft.toeplitz_ones(nux,nuy,1,0,0,1)) +
            db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,1)))
            *F*G;
    } else {
        Y1+=(db_matrix::hadamard(TTO,
            M_fft.toeplitz_ones(nux,nuy,1,0,0,1)) +
            db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,1)))
            *F*G;
    }
    Y1 *=(1.0/omega);

    return Y1;

}
db_matrix &fPMLafter::createY2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y2,
                    db_matrix &M_fft, db_matrix &Q_fft,db_matrix &Qm1_fft,
                    db_matrix &S,
                    db_matrix &F, db_matrix &TTO)
{
    db_matrix pipo=db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,1,0,1,0))
        + db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,2,0));

    Y2 = (isBent?1.0:0.0)*(complex<double>(0,1))*
        db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,0))*F;
    //************************* <----
    if(isBent){
        Y2+=S*((omega*omega)* Q_fft.toeplitz_mod(nux, nuy, 0,0,0,0) -
          pipo*F*F);
    } else {
        Y2+=((omega*omega)* Q_fft.toeplitz_mod(nux, nuy, 0,0,0,0) -
          pipo*F*F);
    }

    pipo.kill();



    Y2 *= (1.0/omega);
    cout<<"Y2... ";
    cout.flush();
    return Y2;
}

db_matrix &fPMLafter::createY3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y3,
                    db_matrix &M_fft, db_matrix &P_fft,db_matrix &Pm1_fft,
                    db_matrix &S,
                    db_matrix &G, db_matrix &TTO)
{
    db_matrix pipo1=db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,1,0,1))
        + db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,0,2));


    if(isBent) {
        Y3 = S*((-omega*omega)* P_fft.toeplitz_mod(nux, nuy, 0,0,0,0)+
            pipo1*(G*G));
    } else {
        Y3 = ((-omega*omega)* P_fft.toeplitz_mod(nux, nuy, 0,0,0,0)+
            pipo1*(G*G));
    }
    pipo1.kill();


    Y3 *=(1.0/omega);

    return Y3;
}
db_matrix &fPMLafter::createY4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{

    if(isBent){
        Y4 =-S*(db_matrix::hadamard(TTO,
            M_fft.toeplitz_ones(nux,nuy,0,1,1,0)) +
            db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,1)))
            *(F*G);
    } else {
        Y4 =-1.0*(db_matrix::hadamard(TTO,
            M_fft.toeplitz_ones(nux,nuy,0,1,1,0)) +
            db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,1)))
            *(F*G);
    }

    Y4 *= (1.0/omega);
    return Y4;
}


/*  Developments for matdev af

*/



db_matrix & fPMLafterOPT::createX1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{
    X1 = db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,1,0,0,1));
    X1 += db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,1));
    X1 *= -1.0;
    X1 *= F;
    X1 *= G;

    if(isBent) {
        db_matrix toto = X1;
        X1 = S;
        X1 *= toto;
        toto = db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,0,1));
        toto *= G;
        toto *= complex<double>(0,1);
        X1 +=toto;
        toto.kill();
    }
/*
    X1 = (isBent?0:1.0)*(complex<double>(0,1))*
          db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,0,1))*G;

    if(isBent) {
        X1-=S*(db_matrix::hadamard(TTR,
            M_fft.toeplitz_ones(nux,nuy,1,0,0,1))
            +db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,1)))
            *F*G;
    } else {
        X1-=(db_matrix::hadamard(TTR,
            M_fft.toeplitz_ones(nux,nuy,1,0,0,1))
            +db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,1)))
            *F*G;

    }
    */

    X1 *=(1.0/omega);
    return X1;
}

db_matrix &fPMLafterOPT::createX2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X2,
                    db_matrix &M_fft, db_matrix &N_fft, db_matrix &S,
                    db_matrix &F, db_matrix &TTR)
{

    X2 = db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,1,0,1,0));
    X2 += db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,2,0));
    X2 *= F;
    X2 *= F;
    X2 += complex<double>(-omega*omega,0) *
            N_fft.toeplitz_mod(nux, nuy, 0,0,0,0);

    if(isBent) {
        db_matrix toto = X2;
        X2 = S;
        X2 *= toto;
        toto = db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,0));
        toto *= F;
        toto *= complex<double>(0,-1.0);
        X2 += toto;
    }
    /*
    db_matrix pipo=db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,1,0,1,0))
       + db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,2,0));


    X2 = (isBent?0:1.0)*(complex<double>(0,-1))*
         db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,0))*F;
    if(isBent) {
        X2+=S*((complex<double>(-omega*omega,0) *
            N_fft.toeplitz_mod(nux, nuy, 0,0,0,0) +
            pipo*(F*F)));
     } else {
        X2+=((complex<double>(-omega*omega,0) *
            N_fft.toeplitz_mod(nux, nuy, 0,0,0,0) +
            pipo*(F*F)));
     }

    pipo.kill();
    */

    X2 *= (1.0/omega);

    return X2;
}

db_matrix &fPMLafterOPT::createX3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X3,
                    db_matrix &O_fft, db_matrix &M_fft, db_matrix &S,
                    db_matrix &G, db_matrix &TTR)
{
    X3 = db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,1,0,1));
    X3 += db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,0,2));
    X3 *= G;
    X3 *= G;
    X3 *= -1.0;
    X3 += (omega*omega)*
            M_fft.toeplitz_mod(nux, nuy, 0,0,0,0);
    if(isBent) {
        db_matrix toto = X3;
        X3 = S;
        X3 *= toto;
    }
    /*
    db_matrix pipo=db_matrix::hadamard(TTR,
        M_fft.toeplitz_ones(nux,nuy,0,1,0,1))
        + db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,0,2));

    if(isBent) {
        X3 = S*((omega*omega)*
            M_fft.toeplitz_mod(nux, nuy, 0,0,0,0) -
            pipo*(G*G));
    } else {
        X3 = ((omega*omega)*
            M_fft.toeplitz_mod(nux, nuy, 0,0,0,0) -
            pipo*(G*G));
    }

    pipo.kill();
    */
    X3 *=(1.0/omega);

    return X3;
}
db_matrix &fPMLafterOPT::createX4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{

    if(isBent) {
        X4 = S*(db_matrix::hadamard(TTR,
            M_fft.toeplitz_ones(nux,nuy,0,1,1,0))+
            db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,1)))
            *(F*G);
    } else {
        X4 =(db_matrix::hadamard(TTR,
            M_fft.toeplitz_ones(nux,nuy,0,1,1,0))+
            db_matrix::hadamard(TTR, M_fft.toeplitz_ones(nux,nuy,0,0,1,1)))
            *(F*G);
    }

    X4 *= (1.0/omega);

    return X4;
}
db_matrix & fPMLafterOPT::createY1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{

    Y1 = db_matrix::hadamard(TTO, M_fft.toeplitz_ones(nux,nuy,1,0,0,1));
    Y1 += db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,1));
    Y1 *=F;
    Y1 *=G;
    if(isBent) {
        db_matrix toto = Y1;
        Y1 = S;
        Y1 *=toto;
        toto = db_matrix::hadamard(TTO, M_fft.toeplitz_ones(nux,nuy,0,0,0,1));
        toto *= complex<double>(0,-1.0);
        toto *= G;
        Y1 +=toto;
    }

/*
    Y1 = (isBent?0:1.0)*(complex<double>(0,-1))*
        db_matrix::hadamard(TTO, M_fft.toeplitz_ones(nux,nuy,0,0,0,1))*G;

    if(isBent) {
        Y1+=S*(db_matrix::hadamard(TTO,
            M_fft.toeplitz_ones(nux,nuy,1,0,0,1)) +
            db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,1)))
            *F*G;
    } else {
        Y1+=(db_matrix::hadamard(TTO,
            M_fft.toeplitz_ones(nux,nuy,1,0,0,1)) +
            db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,1)))
            *F*G;
    }
*/

    Y1 *=(1.0/omega);

    return Y1;

}
db_matrix &fPMLafterOPT::createY2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y2,
                    db_matrix &M_fft, db_matrix &Q_fft, db_matrix &Qm1_fft,
                    db_matrix &S,
                    db_matrix &F, db_matrix &TTO)
{

    Y2 = db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,1,0,1,0));
    Y2+= db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,2,0));
    Y2*=F;
    Y2*=F;
    Y2*=-1.0;
    Y2+=(omega*omega)* Q_fft.toeplitz_mod(nux, nuy, 0,0,0,0);
    if(isBent){
        db_matrix toto = Y2;
        Y2 = S;
        Y2 *=toto;
        toto = db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,0));
        toto *= F;
        toto *= complex<double>(0,1);

        Y2 += toto;
    }


/*
    db_matrix pipo=db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,1,0,1,0))
        + db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,2,0));

    Y2 = (isBent?0:1.0)*(complex<double>(0,1))*
        db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,0))*F;
    // ************************* <----
    if(isBent){
        Y2+=S*((omega*omega)* Q_fft.toeplitz_mod(nux, nuy, 0,0,0,0) -
          pipo*F*F);
    } else {
        Y2+=((omega*omega)* Q_fft.toeplitz_mod(nux, nuy, 0,0,0,0) -
          pipo*F*F);
    }

    pipo.kill();
*/


    Y2 *= (1.0/omega);
    return Y2;
}

db_matrix &fPMLafterOPT::createY3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y3,
                    db_matrix &M_fft, db_matrix &P_fft, db_matrix &Pm1_fft,
                    db_matrix &S,
                    db_matrix &G, db_matrix &TTO)
{
    Y3=db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,1,0,1));
    Y3+= db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,0,2));
    Y3*=G;
    Y3*=G;
    Y3+= (-omega*omega)* P_fft.toeplitz_mod(nux, nuy, 0,0,0,0);


    if(isBent) {
        db_matrix pipo1 = Y3;
        Y3 = S;
        Y3*= pipo1;
    }

    /*
        db_matrix pipo1=db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,1,0,1))
        + db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,0,2));


    if(isBent) {
        Y3 = S*((-omega*omega)* P_fft.toeplitz_mod(nux, nuy, 0,0,0,0)+
            pipo1*(G*G));
    } else {
        Y3 = ((-omega*omega)* P_fft.toeplitz_mod(nux, nuy, 0,0,0,0)+
            pipo1*(G*G));
    }
    pipo1.kill();
*/


    Y3 *=(1.0/omega);

    return Y3;
}
db_matrix &fPMLafterOPT::createY4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{
    Y4 = db_matrix::hadamard(TTO,
            M_fft.toeplitz_ones(nux,nuy,0,1,1,0));
    Y4 += db_matrix::hadamard(TTO,M_fft.toeplitz_ones(nux,nuy,0,0,1,1));

    Y4*=F;
    Y4*=G;
    Y4*= -1.0;


    if(isBent){
        db_matrix toto=Y4;
        Y4 = S;
        Y4 *= toto;
    }


    Y4 *= (1.0/omega);

    return Y4;
}

/** Implementation of :
    matdev nd
    matdev la alpha

*/

db_matrix & fNonDev::createX1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{
    X1 = F;
    X1 *= M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);
    if(isBent) {
        X1 *= S;
    }
    X1 *= TTR;
    X1 *= G;
    X1 *= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);
    X1 *= -1.0;

    X1 *=(1.0/omega);

    return X1;
}

db_matrix &fNonDev::createX2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X2,
                    db_matrix &M_fft, db_matrix &N_fft, db_matrix &S,
                    db_matrix &F, db_matrix &TTR)
{

    X2 = F;
    X2 *= M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);

    if(isBent) {
        X2*= S;
    }
    X2 *= TTR;
    X2 *= F;
    X2 *= M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);

    db_matrix toto = N_fft.toeplitz_mod(nux, nuy, 0,0,0,0);
    if(isBent) {
        toto *= S;
    }
    toto *= (-omega*omega);
    X2 += toto;

    X2 *= (1.0/omega);

    return X2;
}

db_matrix &fNonDev::createX3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X3,
                    db_matrix &O_fft, db_matrix &M_fft, db_matrix &S,
                    db_matrix &G, db_matrix &TTR)
{
    if(isBent) {
        X3 = S;
        X3 *= G;
    } else {
        X3 = G;
    }
    X3 *= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);
    X3*=TTR;
    X3 *=G;
    X3*= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);

    X3 *=-1.0;

    db_matrix toto = M_fft.toeplitz_mod(nux, nuy, 0,0,0,0);
    if(isBent) {
        toto *= S;
    }
    toto *= (omega*omega);
    X3 += toto;


    X3 *=(1.0/omega);

    return X3;
}
db_matrix &fNonDev::createX4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{
    if(isBent) {
        X4 = S;
        X4 *= G;
    } else {
        X4 = G;
    }
    X4 *= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);

    X4 *=TTR;

    X4 *=F;
    X4 *= M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);

    X4 *= (1.0/omega);

    return X4;
}
db_matrix & fNonDev::createY1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{

    Y1 = F;
    Y1 *=  M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);
    if(isBent){
        Y1 *= S;
    }
    Y1 *= TTO;
    Y1 *= G;
    Y1 *= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);

    Y1 *=(1.0/omega);

    return Y1;

}
db_matrix &fNonDev::createY2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y2,
                    db_matrix &M_fft, db_matrix &Q_fft, db_matrix &Qm1_fft,
                    db_matrix &S,
                    db_matrix &F, db_matrix &TTO)
{
    Y2 = F;
    Y2 *= M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);
    if(isBent){
        Y2 *= S;
    }
    Y2 *=TTO;
    Y2 *=F;
    Y2 *= M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);
    Y2 *=-1.0;

    db_matrix toto = Q_fft.toeplitz_mod(nux, nuy,0,0,0,0);
    //cout<<"alpha="<<alpha<<"\n";
    if(alpha>=0) {
        toto *= alpha;
        db_matrix toto1 = Qm1_fft.toeplitz_mod(nux, nuy, 0,0,0,0);
        toto1.invert();
        toto1 *= (1.0-alpha);

        toto += toto1;
        toto1.kill();
    }
    if(isBent) {
        toto *= S;
    }
    toto *= (omega*omega);
    Y2+=toto;
    Y2 *= (1.0/omega);
    return Y2;
}

db_matrix &fNonDev::createY3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y3,
                    db_matrix &M_fft, db_matrix &P_fft, db_matrix &Pm1_fft,
                    db_matrix &S,
                    db_matrix &G, db_matrix &TTO)
{
    if(isBent) {
        Y3 = S;
        Y3 *= G;
    } else {
        Y3 = G;
    }


    Y3 *= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);
    Y3 *= TTO;

    Y3 *= G;
    Y3 *= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);

    db_matrix toto = P_fft.toeplitz_mod(nux, nuy, 0,0,0,0);

    if(alpha>=0) {
        toto *= (1.0-alpha);
        db_matrix toto1 = Pm1_fft.toeplitz_mod(nux, nuy, 0,0,0,0);
        toto1.invert();
        toto1 *= alpha;

        toto += toto1;
        toto1.kill();
    }

    if(isBent) {
        toto *= S;
    }
    toto *= (-omega*omega);
    Y3 += toto;

    Y3 *=(1.0/omega);

    return Y3;
}
db_matrix &fNonDev::createY4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{
    if(isBent) {
        Y4 = S;
        Y4 *= F;
    } else {
        Y4 = F;
    }
    Y4 *= M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);
    Y4 *= TTO;
    Y4 *= G;
    Y4 *= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);
    Y4*=-1.0;

    Y4 *= (1.0/omega);
    return Y4;
}



/** matdev nf normal_vector_x.f3d normal_vector_y.f3d

    Implementation of the normal field vector. See:
    "Normal vector method for convergence improvement using the RCWA for
    crossed gratings", Thomas Schuster, Johannes Ruoff, Norbert Kerwien,
    Stephan Rafler, and Wolfgang Osten  JOSA A, Vol. 24, Issue 9, pp. 2880-2890
    (2007)
    http://dx.doi.org/10.1364/JOSAA.24.002880

    Order on which the matrices are created:

    X2, X4, X1, X3, Y2, Y1, Y3, Y4.

    In Y2, we create Delta_x used also by Y1 (and then destroyed)
    In Y4, we create Delta_y, then used and destroyed by Y3

*/

db_matrix & fNF::createX1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{
    X1 = F;
    X1 *= M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);
    if(isBent) {
        X1 *= S;
    }
    X1 *= TTR;
    X1 *= G;
    X1 *= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);
    X1 *= -1.0;

    X1 *=(1.0/omega);

    return X1;
}

db_matrix &fNF::createX2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X2,
                    db_matrix &M_fft, db_matrix &N_fft, db_matrix &S,
                    db_matrix &F, db_matrix &TTR)
{

    X2 = F;
    X2 *= M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);

    if(isBent) {
        X2*= S;
    }
    X2 *= TTR;
    X2 *= F;
    X2 *= M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);

    db_matrix toto = N_fft.toeplitz_mod(nux, nuy, 0,0,0,0);
    if(isBent) {
        toto *= S;
    }
    toto *= (-omega*omega);
    X2 += toto;

    X2 *= (1.0/omega);

    return X2;
}

db_matrix &fNF::createX3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X3,
                    db_matrix &O_fft, db_matrix &M_fft, db_matrix &S,
                    db_matrix &G, db_matrix &TTR)
{
    if(isBent) {
        X3 = S;
        X3 *= G;
    } else {
        X3 = G;
    }
    X3 *= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);
    X3*=TTR;
    X3 *=G;
    X3*= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);

    X3 *=-1.0;

    db_matrix toto = M_fft.toeplitz_mod(nux, nuy, 0,0,0,0);
    if(isBent) {
        toto *= S;
    }
    toto *= (omega*omega);
    X3 += toto;


    X3 *=(1.0/omega);

    return X3;
}
db_matrix &fNF::createX4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &X4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTR)
{
    if(isBent) {
        X4 = S;
        X4 *= G;
    } else {
        X4 = G;
    }
    X4 *= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);

    X4 *=TTR;

    X4 *=F;
    X4 *= M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);

    X4 *= (1.0/omega);

    return X4;
}
db_matrix & fNF::createY1(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y1,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{

    Y1 = F;
    Y1 *=  M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);
    if(isBent){
        Y1 *= S;
    }
    Y1 *= TTO;
    Y1 *= G;
    Y1 *= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);

    db_matrix toto = Delta_z;
    toto*=Nr;
    toto*=Nz;
    if(isBent){
        toto*=S;
    }
    toto*= (-omega*omega);

    Y1 += toto;

    Y1 *=(1.0/omega);

    return Y1;

}
db_matrix &fNF::createY2(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y2,
                    db_matrix &M_fft, db_matrix &Q_fft, db_matrix &Qm1_fft,
                    db_matrix &S,
                    db_matrix &F, db_matrix &TTO)
{
    Delta_z = Q_fft.toeplitz_mod(nux, nuy,0,0,0,0)-
        (Qm1_fft.toeplitz_mod(nux, nuy,0,0,0,0)).invert();


    Nz = Nz_fft.toeplitz_mod(nux, nuy,0,0,0,0);
    Nr = Nr_fft.toeplitz_mod(nux, nuy,0,0,0,0);


    Y2 = F;
    Y2 *= M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);
    if(isBent){
        Y2 *= S;
    }
    Y2 *=TTO;
    Y2 *=F;
    Y2 *= M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);
    Y2 *=-1.0;

    db_matrix toto = Q_fft.toeplitz_mod(nux, nuy,0,0,0,0);

    toto -= Delta_z*Nz*Nz;

    if(isBent) {
        toto *= S;
    }

    toto *= (omega*omega);
    Y2+=toto;

    Y2 *= (1.0/omega);

    return Y2;
}

db_matrix &fNF::createY3(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y3,
                    db_matrix &M_fft, db_matrix &P_fft, db_matrix &Pm1_fft,
                    db_matrix &S,
                    db_matrix &G, db_matrix &TTO)
{
    Delta_r = P_fft.toeplitz_mod(nux, nuy,0,0,0,0)-
        (Pm1_fft.toeplitz_mod(nux, nuy,0,0,0,0)).invert();

    if(isBent) {
        Y3 = S;
        Y3 *= G;
    } else {
        Y3 = G;
    }


    Y3 *= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);
    Y3 *= TTO;

    Y3 *= G;
    Y3 *= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);

    db_matrix toto = P_fft.toeplitz_mod(nux, nuy, 0,0,0,0);

    toto -= Delta_r*Nr*Nr;

    if(isBent) {
        toto *= S;
    }
    toto *= (-omega*omega);
    Y3 += toto;

    Y3 *=(1.0/omega);

    return Y3;
}
db_matrix &fNF::createY4(double nux, double nuy, double omega,
                    bool isBent,
                    db_matrix &Y4,
                    db_matrix &M_fft, db_matrix &S,
                    db_matrix &F,
                    db_matrix &G, db_matrix &TTO)
{
//JEROME
/*  if(isBent) {
        Y4 = S;
        Y4 *= F;
    } else {
        Y4 = F;
    }
    Y4 *= M_fft.toeplitz_deriv(nux,nuy,1,0);
    Y4 *= TTO;
    Y4 *= G;
    Y4 *= M_fft.toeplitz_deriv(nux,nuy,0,1);
    Y4 *=-1.0;
*/
    if(isBent) {
        Y4 = S;
        Y4 *= F;
    } else {
        Y4 = F;
    }
    Y4 *= M_fft.toeplitz_deriv(nux,nuy,0,1,ksinthetax,ksinthetay);
    Y4 *= TTO;
    Y4 *= G;
    Y4 *= M_fft.toeplitz_deriv(nux,nuy,1,0,ksinthetax,ksinthetay);
    Y4 *=-1.0;


    db_matrix toto = Delta_r;

    toto *=Nr;
    toto *=Nz;


    if(isBent) {
        toto *= S;
    }
    toto *= (omega*omega);

    Y4 += toto;

    Y4 *= (1.0/omega);
    return Y4;
}
