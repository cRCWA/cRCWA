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
