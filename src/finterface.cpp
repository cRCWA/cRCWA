/***************************************************************************
*     CLASS finterface "Field interface"                                   *
*     Davide Bucci, CROMA                                                  *
*     MINATEC-INPG, 3, parvis Luis Neel                                    *
*     38016, Grenoble CEDEX, France                                        *
*                                                                          *
*     bucci@minatec.inpg.fr                                                *
*                                                                          *
****************************************************************************/

#include "finterface.h"
#include "section.h"

/**  From the eigenvector matrices of the electric and magnetic field
    of two sections, calculates the S matrix of the interface between
    the two. For each section, the eigenvalue and eigenvector matrices
    should have been calculated before calling to createInterface.

    @param t the first section
    @param tp1 the second section

*/

finterface& finterface::createInterface(db_matrix &Wt_m1,
    db_matrix &Wtp1, db_matrix &Vt_m1, db_matrix &Vtp1, db_matrix &Wtp1_m1,
    db_matrix &Wt, db_matrix &Vtp1_m1, db_matrix &Vt, db_matrix &Pt,
    db_matrix &Ptp1)
{

    // Tpp=(Wt_m1*Wtp1+Vt_m1*Vtp1);

    Tpp = Wt_m1;
    Tpp *= Wtp1;
    Rpm = Tpp;
    multsum(Vt_m1, Vtp1, Tpp);

    Tpp.invert();

    // Rpm=(-Wt_m1*Wtp1+Vt_m1*Vtp1);

    multsumscale(complex<double>(1.0,0), Vt_m1, Vtp1,
        complex<double>(-1.0,0), Rpm);

    Rpm=Tpp*Rpm;
    Rpm*=Ptp1;

    //Tpp*=2.0;
    //Tpp*=Pt;
    complex<double> deux=2.0;

    Tpp.scalemult(deux, Pt);

    // Tmm=(Wtp1_m1*Wt+Vtp1_m1*Vt);

    Tmm = Wtp1_m1;
    Tmm *= Wt;
    Rmp = Tmm;
    multsum(Vtp1_m1, Vt, Tmm);

    Tmm.invert();

    // Rmp=(-Wtp1_m1*Wt+Vtp1_m1*Vt);

    multsumscale(complex<double>(1.0,0), Vtp1_m1, Vt,
        complex<double>(-1.0,0), Rmp);


    Rmp=Tmm*Rmp;
    Rmp*=Pt;

    //Tmm*=2.0;
    //Tmm*=Ptp1;

    Tmm.scalemult(deux, Ptp1);

    return *this;
}

/** From an array of interfaces, calculates the S matrix of the global structure
    For each interface, the elements of the S matrix should have been created.
    This routine processes every interface in order to calculate the global S
    matrix.

    @param s the array of interfaces to be used for calculations
    @param numOfInterface the number of interfaces used in the array


*/
finterface finterface::createSmatrix(finterface *s, int numOfInterfaces)
{
    finterface F;
    db_matrix I=db_matrix::createUnitMatrix(s[0].Tpp.getNcol(),
        s[0].Tpp.getNcol());

    F.Tpp = s[0].Tpp;
    F.Rpm = s[0].Rpm;
    F.Rmp = s[0].Rmp;
    F.Tmm = s[0].Tmm;

    db_matrix TppN;
    db_matrix RpmN;
    db_matrix RmpN;
    db_matrix TmmN;

    db_matrix pipo;

    s[0].TRpm=s[0].Rpm;
    s[0].TTpp=s[0].Tpp;

    for(int i=1; i<numOfInterfaces-1; ++i) {

        //pipo=I-F.Rpm*s[i].Rmp;

        pipo=I;
        multsumscale(complex<double>(-1.0, 0), F.Rpm,
             s[i].Rmp, complex<double>(1.0, 0), pipo);
        pipo.invert();

        TppN = s[i].Tpp;
        TppN *=pipo;
        RpmN = TppN;
        TppN*= F.Tpp;

        RpmN*=F.Rpm;

        RpmN*=s[i].Tmm;
        RpmN += s[i].Rpm;

        RmpN = F.Tmm;

        //pipo=I-s[i].Rmp*F.Rpm;
        pipo=I;
        multsumscale(complex<double>(-1.0, 0), s[i].Rmp,
            F.Rpm, complex<double>(1.0, 0), pipo);

        //multsum (-s[i].Rmp,F.Rpm, pipo);
        /*multsumscale (complex<double>(-1.0,0.0), s[i].Rmp,F.Rpm,
            complex<double>(1.0,0.0), pipo);*/

        pipo.invert();
        RmpN *=pipo;
        RmpN *=s[i].Rmp;
        RmpN *=F.Tpp;
        RmpN +=F.Rmp;

        TmmN = F.Tmm;
        TmmN*=pipo;
        TmmN*=s[i].Tmm;

        F.Tpp = TppN;   // S21
        F.Rpm = RpmN;   // S22
        F.Rmp = RmpN;   // S11
        F.Tmm = TmmN;   // S12

        s[i].TRpm=RpmN;
        s[i].TTpp=TppN;
    }

    // At the end of the procedure, F will contain the global S matrix.

    return F;
}


/**  Once the S matrices of each interface have been calculated, as well
    as the total calculation window, calculates the field excitation of
    each section. Note that the propagative component of the excitation
    of the first interface must been defined as well as the couterpropagative
    component of the excitation of the last interface. In other words, one
    should set:

    p[0].sWp
    p[numOfSections-1].sWm

    This last expression can seem to be an error. Of course, if we have N
    interfaces between N+1 sections, we will have... N+1 sections!

    @param s the array of interfaces to be used for calculations
    @param numOfSections the number of interfaces used in the array
    @param p the pointer at the section array

*/
void finterface::inpoutp(finterface *s, int numOfSections, section *p)
{

    db_matrix I=db_matrix::createUnitMatrix(s[0].Tpp.getNcol(),
        s[0].Tpp.getNcol());

    /*  Here we start creating counterpropagative fields. We start from
        the right and propagate the counterpropagative excitation back to
        the beginning of the structure.

    */
    db_matrix toto;
    for(int i=numOfSections-2; i>=0; --i) {
        if(i>0) {
            //p[i].sWm = I-s[i].Rmp*s[i-1].TRpm;

            p[i].sWm=I;
            toto = s[i-1].TRpm;
            //toto *=-1.0;
            //multsum(s[i].Rmp, toto, p[i].sWm);

            multsumscale(complex<double>(-1.0,0.0),s[i].Rmp, toto,
                complex<double>(1.0,0.0), p[i].sWm);

            toto.kill();

            p[i].sWm.invert();

            // p[i].sWm *= (s[i].Rmp*s[i-1].TTpp*p[0].sWp+s[i].Tmm*p[i+1].sWm);
            toto = s[i].Rmp;
            toto *= s[i-1].TTpp;
            toto *= p[0].sWp;

            multsum(s[i].Tmm, p[i+1].sWm, toto);
            p[i].sWm *= toto;
            toto.kill();
        } else {
            //p[0].sWm = s[0].Rmp*p[0].sWp+s[0].Tmm*p[1].sWm;
            p[0].sWm = s[0].Rmp;
            p[0].sWm*= p[0].sWp;
            multsum(s[0].Tmm, p[1].sWm, p[0].sWm);
        }
    }

    /* Here we propagate the propagative fields. We start from the left and
       we propagate fields towards the right of the structure.
    */
    for(int i=1; i<numOfSections; ++i) {
        //p[i].sWp = s[i-1].Tpp*p[i-1].sWp+s[i-1].Rpm*p[i].sWm;
        p[i].sWp = s[i-1].Tpp;
        p[i].sWp *= p[i-1].sWp;
        multsum(s[i-1].Rpm, p[i].sWm, p[i].sWp);
    }
}
