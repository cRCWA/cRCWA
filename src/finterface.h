/***************************************************************************
*     CLASS finterface "Field interface"                                   *
*     Davide Bucci, CROMA     april 2010                                   *
*     MINATEC-INPG, 3, parvis Luis Neel                                    *
*     38016, Grenoble CEDEX, France                                        *
*                                                                          *
*     bucci@minatec.inpg.fr                                                *
*                                                                          *
*     Version: 1.0.0                                                       *
****************************************************************************/

// Prevent multiple includes
#ifndef FINTERFACE_H
#define FINTERFACE_H

#include "section.h"

#include "block_matrix.h"

class finterface {
private:
public:

    // The S matrix of the current section.
    db_matrix Tpp;  // S21
    db_matrix Rpm;  // S22
    db_matrix Rmp;  // S11
    db_matrix Tmm;  // S12

    // The global S matrix (from the beginning of the structure to this
    // particular section.
    db_matrix TRpm; // S22
    db_matrix TTpp; // S21

    //finterface &createInterface(class section &, class section &);

    finterface &createInterface(db_matrix &Wt_m1,
    db_matrix &Wtp1, db_matrix &Vt_m1, db_matrix &Vtp1, db_matrix &Wtp1_m1,
    db_matrix &Wt, db_matrix &Vtp1_m1, db_matrix &Vt, db_matrix &Pt,
    db_matrix &Ptp1);


    static finterface createSmatrix(finterface *s,
        int numOfSections);
    static void inpoutp(finterface *s, int numOfInterfaces,
        class section *p);

};

#endif
