/***************************************************************************
*     CLASS section                                                        *
*     Davide Bucci, CROMA     march 2008 - april 2010                      *
*     Jérôme Michallon, CROMA     2012                                     *
*     MINATEC-INPG, 3, parvis Luis Neel                                    *
*     38016, Grenoble CEDEX, France                                        *
*                                                                          *
*     bucci@minatec.inpg.fr                                                *
*                                                                          *
*     Version: 1.2   (bent waveguides and propagation)                     *
****************************************************************************
*    C++ class containing informations about the waveguide structure.      *
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

#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <limits>

#include "block_matrix.h"
#include "parsefile.h"
#include "structure.h"
#include "matrixDev.h"
#include "finterface.h"
#include "parallelize.h"
#include "phys_constants.h"
#include "commands.h"

using namespace std;

/* Version history

    1.2   Feb, 9, 2010, propagation
    1.1.2 Nov, 22 2009, added the Spectrum command
    1.1.1 May, 26 2009, bent waveguide version
    1.1.0 Apr, 1 2009, added coordinates transform PML, see for example [1]
    1.0.3 Feb, 26 2009, corrected a bug in the output of the rectangle command.


 [1] J.P. Hugonin et al. "Fourier modal methods for modelling optical dielectric
     waveguides", Optical and quantum electronics (2005) 37:107-119

 [2] J.P. Hugonin, P. Lalanne "Perfectly matched layers as nonlinear coordinate
     transforms: a generalized formalization", J. Opt. Soc. Am. A, vol. 22, n. 9
     september 2005

 [3] B. Martin "Etude et réalisation d'un spectromètre compact en optique
     intégrée sur verre", PhD report, INP Grenoble, 23 janvier 2009

*/


section::section()
{
    reset();
}

section::section(section &m)
{
    section::operator= (m);
}

section& section::operator= (const section &m)
{
    // Some ideas from there:
    // http://stackoverflow.com/questions/12423058

    // Prevent self-assignment
    if(this==&m)
        return *this;

    section_data::operator= (m);

    // DB: Not a very elegant solution, in my opinion.
    // maybe crtr should point to something stored OUTSIDE class section?

    if (m.crtr == &m.PMLafterOPT)
        crtr = &PMLafterOPT;
    else if (m.crtr == &m.PMLafter)
        crtr = &PMLafter;
    else if (m.crtr == &m.PMLbefore)
        crtr = &PMLbefore;
    else if (m.crtr == &m.NonDev)
        crtr = &NonDev;
    else if (m.crtr == &m.NormalField)
        crtr = &NormalField;
    else if (m.crtr == &m.NonDevSym)
        crtr = &NonDevSym;
    else if (m.crtr == &m.NormalFieldSym)
        crtr = &NormalFieldSym;
    else {
        cerr<<"ensureSize: error cannot do the link for the matrix"
             " development.\n";
    }
    return *this;
}

section::~section()
{

}

void section::reset(void)
{
    father = NULL;
    qx=qy=0;
    gamma=1;
    high_imag=.1;
    low_imag=-.1;
    isBent = false;
    isSubstrateSet=false;
    tot_z=0;
    useAzimuthalOrder=false;
    isBent = false;
    radius=0;

    // Set the default matrix creation strategy for PMLs
    crtr = &PMLafterOPT;
}

/** Calculate the integral of the Poynting vector (power flow).
    @param p The structure
    @param EE matrix associated to magnetic field in the given quota
    @param HH matrix associated to electric field in the given quota
    @param column column of the previous matrices to be considered

*/
double section::integralPoynting(structure *p, db_matrix &EE, db_matrix &HH,
    int column)
{
    // The number of rows of EE is the double of the one of matrix associated
    // to the x and y components of the fields, because in EE the two
    // components are arranged in a single column.
    int nlines = EE.getNrow()/2;

    // Calculation of the mode normalization. First of all, we obtain
    // the Fourier coefficients of all components involved in the
    // z component of the Poynting vector.
    db_matrix Hxc(nlines,1);    // Complex conjugate of magnetic fields.
    db_matrix Hyc(nlines,1);

    db_matrix Ex(nlines,1);     // Electric fields.
    db_matrix Ey(nlines,1);

    // Here we extract the components from the V and W matrices,
    // respectively for the magnetic and electric fields.

    int shift = 0;
    for(int j=0; j<nlines; ++j) {
        Hxc(j,0) = conj(HH(j,column));          // Calculate the conjugate
        Hyc(j,0) = conj(HH(j+nlines,column));   // at the same time, for H.
        Ex(j,0) = EE(j,column);
        Ey(j,0) = EE(j+nlines,column);
    }

    // Here we calculate the Fourier coefficients of the Poynting
    // vector, which for the z component is given by:
    // Sz = Ex conj(Hy) - Ey conj(Hx)
    // where the multiplication becomes a convolution of the Fourier
    // coefficients (thus involving the multiplication by a Toeplitz
    // matrix). The loops represent the calculation of the 0 term of
    // that convolution.

    complex<double> Poynting_0 = 0;

    // In case of symmetry, the computation of the Poynting vector is 2 or 4*
    // times the value of the Poynting vector, except for the constant term.

    complex<double> factorX = complex<double>(1,0);
    complex<double> factorY = complex<double>(1,0);
    // dimx and dimy refer to the number of harmonics considered for x and y
    double dimx;
    double dimy;
    dimx = p->dimx;
    dimy = p->dimy;

    // Correction of the sizes when a symmetry is taken into account.
    if (p->symx == anti_symmetric || p->symx== symmetric) {
        factorX *= complex<double>(2,0);
        dimx = int(floor(double(dimx)/4.0) + 1.0);
    } else {
        dimx = int(floor(double(dimx)/2.0) + 1.0);
    }
    if (p->symy == anti_symmetric || p->symy== symmetric) {
        factorY *= complex<double>(2,0);
        dimy = int(floor(double(dimy)/4.0) + 1.0);
    } else {
        dimy = int(floor(double(dimy)/2.0) + 1.0);
    }

    // Compute the term independent of x and y:
    Poynting_0 += Ex(0,0)*Hyc(0,0)-Ey(0,0)*Hxc(0,0);
    // Compute the Poynting vector for term independent of y

    for(int kk=1; kk<dimx; ++kk) {
        Poynting_0 += factorX*(Ex(kk,0)*Hyc(kk,0)-Ey(kk,0)*Hxc(kk,0));
    }

    for(int ll=1; ll<dimy; ++ll) {
        // Compute the term independent of x:
        Poynting_0 += factorY*(Ex(ll*dimx,0)*Hyc(ll*dimx,0)-Ey(ll*dimx,0)*
            Hxc(ll*dimx,0));
        // Compute the Poynting vector for y and y dependent term
        for(int kk=1; kk<dimx; ++kk) {
            Poynting_0 += factorY*factorX*(Ex(kk+ll*dimx,0)*Hyc(kk+ll*dimx,0)-
                Ey(kk+ll*dimx,0)*Hxc(kk+ll*dimx,0));
        }
    }

    return -0.5*Poynting_0.real();
}

/** Calculate the integral of the Poynting vector (power flow) inside
    a rectangle.
    @param p The structure
    @param EE matrix associated to magnetic field in the given quota
    @param HH matrix associated to electric field in the given quota
    @param column column of the previous matrices to be considered
    @param wx Width of the calculation rectangle
    @param wy Height of the calculation rectangle
    @param px Position of the center of the rectangle (in x)
    @param py Position of the center of the rectangle (in y)
*/
double section::integralPoynting_rectangle(structure *p, db_matrix &EE,
    db_matrix &HH, int column, double wx, double wy, double px, double py)
{
    // dimx and dimy refer to the number of harmonics considered for x and y
    // shiftx, shifty point at the 0 frequency
    int dimx;
    int dimy;
    dimx = p->dimx;
    dimy = p->dimy;

    dimx = int(floor(double(dimx)/2.0) + 1.0);
    dimy = int(floor(double(dimy)/2.0) + 1.0);

    bool xsymmetric=false, ysymmetric=false;
    int shiftx, shifty;

    // if the number of harmonics is not 4*p+1, it can give bad values, this
    // comes from the fact that we need 4 time more harmonics to create a
    // vector of size p+1
    if (((int)((dimx-1)/4)*4+1 != dimx && dimx != 1) ||
         ((int)((dimy-1)/4)*4+1 != dimy && dimy != 1))
        cout << "ERROR: you must use a number of harmonics that follow the"
            " equation" << endl << " harmonics = 4*p+1 with p integer to get"
            " a good value with this function." << endl;

    if (p->symx == anti_symmetric || p->symx== symmetric)
        xsymmetric = true;
    if (p->symy == anti_symmetric || p->symy== symmetric)
        ysymmetric = true;

    shiftx= dimx/4;
    shifty= dimy/4;

    // V contains both the x and y fields one on the top of the other.
    // nlines is therefore the size of each component.

    int nlines = HH.getNrow()/2;

    // Calculation of the mode normalization. First of all, we allocate
    // memory for the Fourier coefficients of all components involved in the
    // z component of the Poynting vector.

    db_matrix Hxc(nlines,1);    // Conjugate of magnetic fields
    db_matrix Hyc(nlines,1);

    db_matrix Ex(nlines,1);     // Electric fields
    db_matrix Ey(nlines,1);

    // Here we extract the components from the V and W matrices,
    // respectively for the magnetic and electric fields.

    for(int j=0; j<nlines; ++j) {
        Hxc(j,0) = conj(HH(j,column));
        Hyc(j,0) = conj(HH(j+nlines,column));
        Ex(j,0) = EE(j,column);
        Ey(j,0) = EE(j+nlines,column);
    }


    // Extend the vector in case of symmetry so that we end up with the same
    // case as without symmetries:
    db_matrix Hxc_ext = Hxc.vector2fft(p->symx,p->symy,dimx,dimy,false,
        false,false).fft2vector(false,false,0);
    db_matrix Hyc_ext = Hyc.vector2fft(p->symx,p->symy,dimx,dimy,false,
        false,false).fft2vector(false,false,0);


    // Since the Toeplitz matrices need twice more harmonics, here we reduce
    // the size of the vector by a factor 2.

    dimx = dimx/2+1;

    db_matrix Hxc2(dimx*(dimy/2+1),1);
    db_matrix Hyc2(dimx*(dimy/2+1),1);

    for(int j=0; j<dimy/2+1; ++j) {
        for(int i=0; i<dimx; ++i) {
            Hxc2(j*dimx+i,0) = Hxc_ext((j+shifty)*(dimx*2-1)+(i+shiftx),0);
            Hyc2(j*dimx+i,0) = Hyc_ext((j+shifty)*(dimx*2-1)+(i+shiftx),0);
        }
    }
    //free the memory
    Hxc.kill();
    Hyc.kill();

    Hxc_ext.kill();
    Hyc_ext.kill();

    // Here we calculate the Fourier coefficients of the Poynting
    // vector, which for the z component is given by:
    // Sz = Ex conj(Hy) - Ey conj(Hx),
    // where the multiplication becomes a convolution of the Fourier
    // coefficients (thus involving the multiplication by a Toeplitz
    // matrix).

    // In case of symmetry, the Poynting vector is always even because Ex and
    // Hy, Ey and Hx have the same parity

    // Compute the convolution of a field times a conjugate of a second field
    // Hence, the operation is the following (a sign '+' instead of a sign '-'
    // for a classical Toeplitz matrix):
    // Ci = sum (Bi+ax) Aax

    // Convert the vector Ex into usual FFT matrix and then order it.
    db_matrix FFT_ordered = Ex.vector2fft(p->symx,p->symy,dimx*2-1,dimy,false,
                         false,false).fftshift().transpose();

    // Creation of a vector containing 1D Toeplitz matrices:
    db_matrix Vector(dimy*dimx, dimx);
    db_matrix Toeplitz1Dp;

    for(int i=0;i< dimy; ++i) {
        Toeplitz1Dp = FFT_ordered.toeplitz1D(1, i, false, +1);
        Vector.copy(Toeplitz1Dp, i*dimx,0);
    }
    // Creation of the Toeplitz matrix using a vector of 1D Toeplitz matrices
    // (for the y components)
    db_matrix Ex_toep = Vector.toeplitz1D(dimx, 0, false, +1);
    // Do the same for Ey:

    // Convert the vector Ey into usual FFT matrix and then order it.
    FFT_ordered = Ey.vector2fft(p->symx,p->symy,dimx*2-1,dimy,false,
                         false,false).fftshift().transpose();

    for(int i=0;i< dimy; ++i) {
        Toeplitz1Dp = FFT_ordered.toeplitz1D(1, i, false, +1);
        Vector.copy(Toeplitz1Dp, i*dimx,0);
    }

    // Creation of the Toeplitz matrix using a vector of 1D Toeplitz matrices
    // (for the y components)
    db_matrix Ey_toep = Vector.toeplitz1D(dimx, 0, false, +1);

    // Free memory:
    FFT_ordered.kill();
    Toeplitz1Dp.kill();
    Vector.kill();

    // Compute the Poynting vector. The Toeplitz matrices allow for taking
    // into account that this is done in the Fourier space and it is thus
    // a real convolution which is done.
    db_matrix Sz = Ex_toep * Hyc2 - Ey_toep * Hxc2;

    // The following loop is the calculation of
    // the integral of the Poynting vector in the given rectangle.

    int sx,sy;
    double nux=2.0*M_PI/p->tot_x;
    double nuy=2.0*M_PI/p->tot_y;
    complex<double> Poynting_0 = 0;
    complex<double> contrib_y = 0;
    complex<double> contrib_x = 0;

    px=px+p->tot_x/2;
    py=py+p->tot_y/2;


    dimy = dimy/2+1;

    // The integral is done in the Fourier domain, by considering that
    // integrating the Poynting vector in a rectangular domain is there
    // represented by the integration of the harmonics of the Poynting
    // vector, convoluted with cardinal sinuses which are the Fourier
    // transform of the rectangular integration window.

    // For more details, see Davide's handbook AFMM 2, page 64 and 80-81

    for(int ll=0; ll<dimy; ++ll) {
        // Compute the term independent of x at first.

        // Harmonic index for y.
        sy = ll - shifty;
        if (sy !=0)
            contrib_y = sin((double)sy*nuy*wy/2.0)/((double)sy*nuy*wy/2.0);
        else
            contrib_y = 1;

        contrib_y *= exp(-complex<double>(0,1)*(double)sy*nuy*py);

        // compute the Poynting vector for y and y dependent term
        for(int kk=0; kk<dimx; ++kk) {
            // Harmonic index for x.
            sx = kk - shiftx;
            if (sx !=0)
                contrib_x=Sz(kk+ll*dimx,0)*
                    sin((double)sx*nux*wx/2.0)/((double)sx*nux*wx/2.0);
            else
                contrib_x = 1;

            contrib_x *= exp(-complex<double>(0,1)*(double)sx*nux*px);

            Poynting_0 +=  Sz(kk+ll*dimx,0)*contrib_y*contrib_x;
        }
    }

    return -0.5*Poynting_0.real() *wx*wy;
}

/** Specify that in this section, an anisotropic PML has to be used.
    @param wx the width of the PML region
    @param wy the height of the PML region
    @param alpha the complex coefficient.
*/
void section::set_pml(double wx, double wy, complex <double> alpha)
{
    if(!isSubstrateSet) {
        throw parsefile_commandError("pml: The substrate should be set"
            " before trying to define structures.");
    }

    // Add anisotropic PMLs
    PML_definition(wx, wy, alpha);
    cout << "Added PML layers alpha: "<<alpha<<", x thickness: "<<wx<<
            " m, y thickness: "<<wy<<" m.\n";
}

/** Specify that for the current section, a coordinate-transform PML strategy
    has to be applied. See J. Hugonin, P. Lalanne, I. Villar, and I. Matias,
    "Fourier modal methods for modeling optical dielectric waveguides,"
    Optical and quantum electronics, vol. 37, no. 1, pp. 107-119, 2005.

    @param qdx the TOTAL width of the PML region
    @param qdy the TOTAL height of the PML region
    @param g the complex PML factor
*/
void section::set_pml_transf(double qdx, double qdy, complex<double> g)
{
        if(!isSubstrateSet) {
            throw parsefile_commandError("pml_transf: The substrate should"
                " be set before trying to define structures.");
        }

        qx=qdx;
        qy=qdy;
        gamma=g;

        cout <<"Size of the coordinates transform region: "<<qdx
             <<" m x "<<qdy<<" m, complex factor: "<<gamma<<"\n";

}

/** Set the lowest bounds for the selection of "interesting" modes by solve
    on the basis of the effective refractive index.
    @param n_g the minimum real and imaginary part for a mode to be considered
        "interesting" after the mode search.
*/
void section::set_lowindex(complex<double>n_l)
{
    n_0=n_l.real();
    low_imag=n_l.imag();
    cout << "Lowest index real part set to: "<<n_l.real()<<" \n";
    cout << "Lowest index imaginary part set to: "<<low_imag<<"\n";
}

/** Set the highest bounds for the selection of "interesting" modes by solve
    on the basis of the effective refractive index.
    @param n_g the highersr real and imaginary part for a mode to be considered
        "interesting" after the mode search.
*/
void section::set_highindex(complex<double>n_h)
{
    n_g=complex<double>(n_h.real());
    high_imag=n_h.imag();
    cout << "Highest index real part set to: "<<n_h.real()<<" \n";
    cout << "Highest index imaginary part set to: "<<high_imag<<" \n";
}

/** Add a rectangular core in the current section.
    @param n_g the refractive index of the core
    @param wx the width of the rectangle
    @param wy the height of the rectangle
    @param px the horisontal position of the rectangle
    @param py the vertical position of the rectangle
*/
void section::add_rectangle(complex<double>n_g, double wx, double wy,
    double px, double py)
{
    if(!isSubstrateSet) {
        throw parsefile_commandError("rectangle: The substrate should be"
            " set before trying to define structures.");
    }

    // Save the rectangle in the list
    rect_wg t;

    t.cx=px;
    t.cy=py;
    t.h=wy;
    t.w=wx;
    t.n=n_g;

    rects.push_back(t);
    structure *p = father;

    P_fft+=(n_g*n_g-n_0*n_0)*EPS_0*
        db_matrix::rectangleFT(wx, wy, p->tot_x, p->tot_y, p->dimx,
        p->dimy, p->tot_x/2+px, p->tot_y/2+py);
    Q_fft+=(n_g*n_g-n_0*n_0)*EPS_0*
        db_matrix::rectangleFT(wx, wy, p->tot_x, p->tot_y, p->dimx,
        p->dimy, p->tot_x/2+px, p->tot_y/2+py);
    R_fft+=(n_g*n_g-n_0*n_0)*EPS_0*
        db_matrix::rectangleFT(wx, wy, p->tot_x, p->tot_y, p->dimx,
        p->dimy, p->tot_x/2.0+px, p->tot_y/2.0+py);

    Pm1_fft+=(1.0/(n_g*n_g)-1.0/(n_0*n_0))/EPS_0*
        db_matrix::rectangleFT(wx, wy, p->tot_x, p->tot_y, p->dimx,
        p->dimy, p->tot_x/2+px, p->tot_y/2+py);

    Qm1_fft+=(1.0/(n_g*n_g)-1.0/(n_0*n_0))/EPS_0*
        db_matrix::rectangleFT(wx, wy, p->tot_x, p->tot_y, p->dimx,
        p->dimy, p->tot_x/2+px, p->tot_y/2+py);
    Rm1_fft+=(1.0/(n_g*n_g)-1.0/(n_0*n_0))/EPS_0*
        db_matrix::rectangleFT(wx, wy, p->tot_x, p->tot_y, p->dimx,
        p->dimy, p->tot_x/2.0+px, p->tot_y/2.0+py);

    // Save the index having maximum value of the real part.
    if ((n_g.real()) > (this->n_g.real())) this->n_g=n_g;

    cout << "Generated rectangular section waveguide with core refractive"
        " index: "<<n_g<< "\nSize "<<wx<<" m x "<<wy<<" m at location "
        <<px<<" m x "<<py<<" m\n";
}

/* Set a value for the substrate of this section.
    @param subs the complex value of the substrate.
*/
void section::set_substrate(complex<double>subs)
{
    structure *p = getFather();
    if (p->tot_x==0 || p->tot_y==0)
        throw parsefile_commandError("substrate: calculation window size"
            " not defined (use command size).");

    if (p->dimx==0 || p->dimy==0)
        throw parsefile_commandError("substrate: number of harmonics not"
            " defined (use command harmonics).");
    
    n_0=subs;
    cout << "Substrate refractive index set to: "<<n_0<<"\n";

    cout.flush();
    // p contains the current object and should therefore be a structure
    M_fft=MU_0 * db_matrix::rectangleFT(p->tot_x, p->tot_y,
        p->tot_x, p->tot_y, p->dimx, p->dimy, p->tot_x/2, p->tot_y/2);
    Mm1_fft=1.0/MU_0 * db_matrix::rectangleFT(p->tot_x, p->tot_y,
        p->tot_x, p->tot_y, p->dimx, p->dimy, p->tot_x/2, p->tot_y/2);

    N_fft=MU_0* db_matrix::rectangleFT(p->tot_x, p->tot_y,
        p->tot_x, p->tot_y, p->dimx, p->dimy, p->tot_x/2, p->tot_y/2);
    Nm1_fft=1.0/MU_0* db_matrix::rectangleFT(p->tot_x, p->tot_y,
        p->tot_x, p->tot_y, p->dimx, p->dimy, p->tot_x/2, p->tot_y/2);

    O_fft=MU_0* db_matrix::rectangleFT(p->tot_x, p->tot_y,
        p->tot_x, p->tot_y, p->dimx, p->dimy, p->tot_x/2, p->tot_y/2);
    Om1_fft=1.0/MU_0* db_matrix::rectangleFT(p->tot_x, p->tot_y,
        p->tot_x, p->tot_y, p->dimx, p->dimy, p->tot_x/2, p->tot_y/2);

    P_fft=n_0*n_0*EPS_0*db_matrix::rectangleFT(
        p->tot_x, p->tot_y, p->tot_x, p->tot_y, p->dimx, p->dimy,
        p->tot_x/2, p->tot_y/2);
    Pm1_fft=1.0/(n_0*n_0*EPS_0)*
        db_matrix::rectangleFT(p->tot_x, p->tot_y, p->tot_x, p->tot_y,
        p->dimx, p->dimy, p->tot_x/2, p->tot_y/2);

    Q_fft=n_0*n_0*EPS_0* db_matrix::rectangleFT(
        p->tot_x, p->tot_y, p->tot_x, p->tot_y, p->dimx, p->dimy,
        p->tot_x/2, p->tot_y/2);
    Qm1_fft=1.0/(n_0*n_0*EPS_0)*
        db_matrix::rectangleFT(p->tot_x, p->tot_y, p->tot_x, p->tot_y,
        p->dimx, p->dimy, p->tot_x/2, p->tot_y/2);

    R_fft=n_0*n_0*EPS_0* db_matrix::rectangleFT(
        p->tot_x, p->tot_y, p->tot_x, p->tot_y, p->dimx, p->dimy,
        p->tot_x/2, p->tot_y/2);
    Rm1_fft=1.0/(n_0*n_0*EPS_0)*
        db_matrix::rectangleFT(p->tot_x, p->tot_y, p->tot_x, p->tot_y,
        p->dimx, p->dimy, p->tot_x/2, p->tot_y/2);

    isSubstrateSet=true;
}

/** Set the bending radius of this section.
    @param r the bending radius.
*/
void section::set_bend(double r)
{
    radius=r;
    isBent = (r!=0);

    cout << "Bending radius set to: "<<radius<<  " m" <<
             ((r==0)?" (straight waveguide).\n":".\n");
}

/** Execute the order command.
    @param min the minimum order for the retained interesting modes.
    @param max the maximum order for the retained interesting modes.
*/
void section::do_order(double min, double max)
{
    useAzimuthalOrder = true;
    minAzimuthalOrder = min;
    maxAzimuthalOrder = max;

    cout << "I will retain modes on the basis of their azimuthal order"
        " value.\n";
    cout << "min: "<< min<< " max: "<<max<<"\n";
}

/** Execute the outgmodes command.
    @param sj number of points in the x direction
    @param si number of points in the y direction
    @param ftype transverse field component that must be drawn.
    @param erimco type of the output. If the type is S (storage), the file
        name is ignored and a list of matrices containing the modes is
        returned.
    @param fname file name. Ignored if erimco is S.
    @return a list of matrices containing the mode field representation if the
        erimco is equal to S (store). It is an empty list in the other cases.
*/
list<db_matrix> section::do_outgmodes(int sj, int si, char *ftype, 
    rimco_e rimc, const char *fname)
{
    list<db_matrix> storage;
    double f=CELERITY/father->lambda;
    double omega=2.0*M_PI*f;
    double k_0=omega/CELERITY;

    bool shiftToY = false;
    bool calcH = false;
    bool calcD = false;
    bool calcz = false;
    bool improve_representation=false;

    field_e field;
    
    if(!isSubstrateSet) {
        throw parsefile_commandError("outgmodes: The substrate should be"
            " set before trying to write modal eigenfields.");
    }

    int snux=father->dimx/2+1;
    int snuy=father->dimy/2+1;

    //if (si<snuy) si=snuy;
    //if (sj<snux) sj=snux;

    cout << "Fields represented with: "<<sj<<" x "<<si<<" points.\n";

    // TODO: more or less the same code in propagation routines. Should
    // be kept together, if possible.
    if(strcmp(ftype,"Ex")==0) {
        shiftToY = false;
        calcH = false;
        calcz = false;
    } else if(strcmp(ftype, "Ey")==0) {
        shiftToY = true;
        calcH = false;
        calcz = false;
    } else if(strcmp(ftype, "Ez")==0) {
        shiftToY = false;
        calcH = false;
        calcz = true;
    } else if(strcmp(ftype,"Hx")==0) {
        shiftToY = false;
        calcH = true;
        calcz = false;
    } else if(strcmp(ftype, "Hy")==0) {
        shiftToY = true;
        calcH = true;
        calcz = false;
    } else if(strcmp(ftype, "Hz")==0) {
        shiftToY = false;
        calcH = true;
        calcz = true;
    } else if(strcmp(ftype, "Dy")==0) {
        shiftToY = true;
        calcH = false;
        calcD = true;
        calcz = false;
    } else if(strcmp(ftype, "Dx")==0) {
        shiftToY = false;
        calcH = false;
        calcD = true;
        calcz = false;
    }else if(strcmp(ftype, "Dz")==0) {
        shiftToY = false;
        calcH = false;
        calcD = true;
        calcz = true;
    }else if(strcmp(ftype, "Ex2")==0) {
        shiftToY = false;
        calcH = false;
        calcD = false;
        calcz = false;
        improve_representation = true;
    }else if(strcmp(ftype, "Ey2")==0) {
        shiftToY = true;
        calcH = false;
        calcD = false;
        calcz = false;
        improve_representation = true;
    }else if(strcmp(ftype, "Ez2")==0) {
        shiftToY = false;
        calcH = false;
        calcD = false;
        calcz = true;
        improve_representation = true;
    }else {
        throw parsefile_commandError("outgmodes: unrecognized field type."
            " Should be {Ex|Ey|Ez|Ex2|Ey2|Ez2|Hx|Hy|Hz|Dx|Dy|Dz}.");
    }


    // Check if the C matrix is empty
    if (W.isEmpty()) {
        throw parsefile_commandError(
            "outgmodes: the problem must be set and modes calculated with"
            " solve before trying to describe them on a file.");
    }
    double dx=father->tot_x / sj;
    double dy=father->tot_y / si;

    db_matrix out(si, sj);
    int mi=0;
    int mm, ll;

    bool interesting=false;
    complex<double>e;
    complex<double>g;
    double order;
    double quality;
    double r=0;


    // compute epsilon and mu constants if needed
    db_matrix epsilonxy;
    db_matrix epsz;
    db_matrix muz;

    commands::getEpsilonAndMu(*this, calcD, calcz, calcH,improve_representation,
        epsilonxy, epsz, muz);

    // The V matrix will be used if we want to compute Ez,Ez2,Dz, Hx and Hy
    // if the V matrix does not exist, it is created
    if ((calcH && !calcz) || (calcz && !calcH)) {
        if (father->HxField && father->HyField) {
            if (V.getNcol() == 0 && V.getNrow() == 0) {
                V = assembleVmatrix();
            }
        } else {
            throw parsefile_commandError("outgmodes: You need to use 'wants "
                "both' to output the required mode");
        }
    }
    // Find guided modes and get their profile
    for (int jj=0; jj<B.getNrow();++jj) {
        e=B(jj,jj);
        g = father->getEffectiveIndex(e);
        order = father->getOrder(e, radius);
        quality = father->getQuality(e);

        /* We need to decide wheter a given mode is interesting enough to
            retain its characteristics. */
        if(useAzimuthalOrder) {
            /* The first possibility (useful for bent waveguides) is to
                take the azimuthal order and see if it is inside the
                defined window */

            interesting= order >= minAzimuthalOrder &&
                order <= maxAzimuthalOrder;
        } else {
            /* If not, we will do the standard way, by taking the effective
                index and comparing it to the window defined by the minimum
                and maximum refractive indices */

            interesting=
                g.real()> n_0.real() &&
                g.real()< n_g.real() &&
                g.imag()< high_imag &&
                g.imag()> low_imag;
            r = g.real();
        }
        // Search for modes in the interval of interest
        if(interesting) {
            db_matrix select(W.getNrow(), 1);

            // Select the jj-order mode
            // so that W times jj selects the mode
            select(jj, 0)=1;

            // We compute the field:
            db_matrix unused;

            out = commands::getField(*this, select,
                shiftToY, calcH, snux, snuy, sj, si,
                calcD, epsilonxy, epsz, muz,calcz,improve_representation,
                false, unused);

            // We can then write the results on the output file
            if (rimc == O) {
                // Optiwave format
                char s[256];
                sprintf(s, "%s_%s_%s_%d.f3d", fname, ftype, "o", mi++);
                FILE *f=fopen(s,"w");
                father->fileoutputOW(out, f, dx*1e6, dy*1e6);
                cout << "Mode file written (Optiwave format): "<< s<<"\n";
            } else if(rimc == S){
                // Store the results in a set of matrices
                storage.push_back(out);
            } else {
                // Gnuplot format. Put a small header.
                char s[256];
                const char *erimco[]={"r","i","m","c","o","s"};
                sprintf(s, "%s_%s_%s_%d.mode", fname, ftype, erimco[rimc],
                    mi++);
                FILE *f=fopen(s,"w");

                if(useAzimuthalOrder) {
                    fprintf(f, "# order=%lf, Q= %lf\n", order, quality);
                } else {
                    fprintf(f, "# neff=(%lf, %e)\n",g.real(), g.imag());
                }
                father->fileoutputGP(out, f, rimc, dx, dy);
                cout << "Mode file written (Gnuplot format): "<< s<<"\n";
            }
        }
    }
    return storage;
}

/**
    Obtain a representation of the refractive index distribution in this
    section of the waveguide
    @param sj number of columns of the matrix to be obtained (x size)
    @param si number of rows of the matrix to be obtained (y size)
    @param type the type of the output matrix.
    @return a matrix containing the structure.
    
*/
db_matrix section::do_inpstruct(int sj, int si, char* otype)
{
    db_matrix out(si, sj);
    structure *p=father;

    cout << "Structure represented with: "<<sj<<" x "<<si<<" points\n";

    if(strcmp(otype,"i")==0) {
        out=P_fft.zero_pad(si,sj).fft2();
        out/=EPS_0;
        for(int i=0; i<out.getNrow(); ++i) {
            for(int j=0; j<out.getNcol(); ++j) {
                out(i,j)=sqrt(out(i,j));
            }
        }
    } else if(strcmp(otype,"ex")==0) {
        out=P_fft.zero_pad(si,sj).fft2();
    } else if(strcmp(otype,"ey")==0) {
        out=Q_fft.zero_pad(si,sj).fft2();
    } else if(strcmp(otype,"ez")==0) {
        out=R_fft.zero_pad(si,sj).fft2();
    } else if(strcmp(otype,"mux")==0) {
        out=M_fft.zero_pad(si,sj).fft2();
    } else if(strcmp(otype,"muy")==0) {
        out=N_fft.zero_pad(si,sj).fft2();
    } else if(strcmp(otype,"muz")==0) {
        out=O_fft.zero_pad(si,sj).fft2();
    } else if(strcmp(otype,"im")==0 || strcmp(otype,"imo")==0) {
        db_matrix i(si, sj);
        out=i;
        double px, py;
        complex<double> n;
        for (int l=0; l<out.getNcol(); ++l) {
            for(int m=0; m<out.getNrow(); ++m) {
                px=(double(l)*p->tot_x)/out.getNcol()-p->tot_x/2.0;
                py=(double(m)*p->tot_y)/out.getNcol()-p->tot_y/2.0;

                n=expectedRefractiveIndex(px, py);

                out(m,l) = n;
            }
        }
    }  else {
        throw parsefile_commandError("inpstruct: unrecognized output type."
            " Should be {i|im|ex|ey|ez|mux|muy|muz}.");
    }
    return out;
}


/** Store a refractive index distribution
    @param rd the matrix to be used.
*/
void section::store_refractive_index(db_matrix rd)
{
    // Store the background in the index database
    background=rd;

    // Here we know the size of the problem. Let's try to allocate
    // a matrix. If there is an error, the user should be warned

    unsigned int nx = rd.getNcol();
    unsigned int ny = rd.getNrow();

    db_matrix index2(ny,nx);
    db_matrix index2m1(ny,nx);

    // The SUBSTRATE command does allocate the correct amount of
    // space for the harmonics describing the structure.

    const char *vv[]={"substrate","1","0"};
    commands::c_substrate(father, 3, const_cast<char **>(vv));

    for(unsigned int i=0; i<ny; ++i) {
        for(unsigned int j=0; j<nx; ++j) {
            index2(i,j)=rd(i,j)*rd(i,j)*EPS_0;
            index2m1(i,j)=1.0/index2(i,j);
        }
    }
    index2 *= 1.0/nx/ny;
    index2m1 *= 1.0/nx/ny;

    /* Here the index matrix should contain the index distribution
       as read from the input file. We will take a 2D Fast Fourier
       Transform and we will take a suffient number of harmonics.
       In a lot of practical cases, the number of retained harmonics
       will be lower than the points calculated from this FFT.
       In such a case, the series will be truncated by retaining
       only the number of harmonics specified for the calculations
       (see command HARMONICS). This is equivalent to a low-pass
       filter on the spatial frequencies of the structure.
    */

    db_matrix index2_fft = index2.ifft2();
    db_matrix index2m1_fft = index2m1.ifft2();

    int shx=(father->dimx-1)/2+1;
    int shy=(father->dimy-1)/2+1;

    /* Test whether the file has a sufficient number of points.
       Probably, it could be possible to implement a zero padding
       strategy in order to let the user to load a file having a
       number of points less than the number of harmonics currently
       employed.
    */
    if(rd.getNcol()<shx || rd.getNrow()<shy) {
        throw parsefile_commandError("indfile: index matrix should have more "
            "points than the current number of harmonics.\n");
    }

    // positive x and y frequencies
    for(unsigned int i=0; i<(father->dimy+1)/2; ++i) {
        for(unsigned int j=0; j<(father->dimx+1)/2; ++j) {

            P_fft(i,j) = index2_fft(i,j);
            Q_fft(i,j) = index2_fft(i,j);
            R_fft(i,j) = index2_fft(i,j);

            Pm1_fft(i,j) = index2m1_fft(i,j);
            Qm1_fft(i,j) = index2m1_fft(i,j);
            Rm1_fft(i,j) = index2m1_fft(i,j);
        }
    }

    // positive x, negative y
    for(unsigned int i=0; i<(father->dimy-1)/2; ++i) {
        for(unsigned int j=0; j<(father->dimx+1)/2; ++j) {

            P_fft(i+shy,j) = index2_fft(i+ny-shy+1,j);
            Q_fft(i+shy,j) = index2_fft(i+ny-shy+1,j);
            R_fft(i+shy,j) = index2_fft(i+ny-shy+1,j);

            Pm1_fft(i+shy,j) = index2m1_fft(i+ny-shy+1,j);
            Qm1_fft(i+shy,j) = index2m1_fft(i+ny-shy+1,j);
            Rm1_fft(i+shy,j) = index2m1_fft(i+ny-shy+1,j);
        }
    }

    // negative x, positive y
    for(unsigned int i=0; i<(father->dimy+1)/2; ++i) {
        for(unsigned int j=0; j<(father->dimx-1)/2; ++j) {

            P_fft(i,j+shx) = index2_fft(i,j+nx-shx+1);
            Q_fft(i,j+shx) = index2_fft(i,j+nx-shx+1);
            R_fft(i,j+shx) = index2_fft(i,j+nx-shx+1);

            Pm1_fft(i,j+shx) = index2m1_fft(i,j+nx-shx+1);
            Qm1_fft(i,j+shx) = index2m1_fft(i,j+nx-shx+1);
            Rm1_fft(i,j+shx) = index2m1_fft(i,j+nx-shx+1);
        }
    }

    // negative x, negative y
    for(unsigned int i=0; i<(father->dimy-1)/2; ++i) {
        for(unsigned int j=0; j<(father->dimx-1)/2; ++j) {

            P_fft(i+shy,j+shx) = index2_fft(i+ny-shy+1,
                j+nx-shx+1);
            Q_fft(i+shy,j+shx) = index2_fft(i+ny-shy+1,
                j+nx-shx+1);
            R_fft(i+shy,j+shx) = index2_fft(i+ny-shy+1,
                j+nx-shx+1);

            Pm1_fft(i+shy,j+shx) =
                index2m1_fft(i+ny-shy+1,j+nx-shx+1);
            Qm1_fft(i+shy,j+shx)=
                index2m1_fft(i+ny-shy+1,j+nx-shx+1);
            Rm1_fft(i+shy,j+shx) =
                index2m1_fft(i+ny-shy+1,j+nx-shx+1);
        }
    }
}

/** Assemble the W matrix which represent the propagation operator in the
   Fourier space for the electric field.
*/
db_matrix &section::calcWmatrix(double omega)
{
    // Spatial frequency calculation
    double nux=2.0*M_PI/father->tot_x;
    double nuy=2.0*M_PI/father->tot_y;

    if(isBent && radius == 0 && father->calcPropagation) {
        throw parsefile_commandError(
            "solve: I can not propagate with zero radius.\n");
    }


    /* To obtain a modified bloc-Toeplitz matrix, we need to multiply element
       by element (Hadamard product) the original bloc-Toeplitz matrix times
       the modification matrices.

       For example, if A contains the Fourier components of a vector,
       the following instruction will calculate the bloc-Toeplitz matrix:

       A.toeplitz_mod(nux, nuy, 0,0,0,0);

       The modification matrices are calculated the very same way from a "ones"
       matrix. This allows a greather flexibility in defining modified matrices
       and allows to apply modificators even when the original matrix is not
       a bloc-Toeplitz matrix (for example, this is done when applying the
       inverse rule, since the inverse of a bloc-Toeplitz matrix is not a bloc-
       Toeplitz matrix).

       If we consider the first order problem, we use the following A-matrix
       definition:

      [Ex]     [Ex]
      [Ey]     [Ey]
      [Hx] = A [Hx]
      [Hy]     [Hy]

      Where the A matrix is the following bloc-matrix:

        [  0    0      X1     X2   ]
        [  0    0      X3     X4   ]
      A=[  Y1   Y2     0      0    ]
        [  Y3   Y4     0      0    ]

    */

    // Types of symmetry available (J.M.)
    typedef enum symmetry_e_t2{X,Y,none,both} symmetry_e2;

    // Symmetry along x,y, both or none
    symmetry_e2 symmetry = none;
    bool xsymmetry;
    bool ysymmetry;
    int factor;
    int MfactorX, MfactorY;
    int NfactorX, NfactorY;
    int OfactorX, OfactorY;
    int PfactorX, PfactorY;
    int QfactorX, QfactorY;
    int RfactorX, RfactorY;
    symmetry_e symx;
    symmetry_e symy;
    bool NrIsOdd, NzIsOdd;

    db_matrix S;
    db_matrix TTO;
    db_matrix TTR;

    symx = father->symx;
    symy = father->symy;

    // depending on the symmetry, the construction of the matrix will be
    // different
    if (symx == previous_version && symy == previous_version) {
        TTR = R_fft.toeplitz_mod(nux, nuy, 0,0,0,0).invert();

        // To save memory space, we create the S matrix only if needed.
        if (isBent) {
            S = bending(radius).toeplitz_mod(nux, nuy, 0,0,0,0);
        } else {
            S = db_matrix::createUnitMatrix(TTR.getNrow(), TTR.getNcol());
        }

        db_matrix F=transfPML(transf_x).toeplitz_mod(nux, nuy, 0,0,0,0);
        waitSemaphoreIO();
        cout << "Created ";
        postSemaphoreIO();

        X2 = crtr->createX2(nux, nuy, omega, isBent, X2, M_fft, N_fft,
            S, F, TTR);
        db_matrix G=transfPML(transf_y).toeplitz_mod(nux, nuy, 0,0,0,0);
        
        //G.windowing();
        X4 = crtr->createX4(nux, nuy, omega, isBent, X4, M_fft, S, F, G, TTR);
        X1 = crtr->createX1(nux, nuy, omega, isBent, X1, M_fft, S, F, G, TTR);
        X3 = crtr->createX3(nux, nuy, omega, isBent, X3, M_fft, M_fft,
            S, G, TTR);
        TTR.kill();

        TTO=O_fft.toeplitz_mod(nux,nuy,0,0,0,0).invert();
        Y2 = crtr->createY2(nux, nuy, omega, isBent, Y2, M_fft, Q_fft,
            Qm1_fft, S, F, TTO);
        Y1 = crtr->createY1(nux, nuy, omega, isBent, Y1, M_fft, S, F, G, TTO);
        Y3 = crtr->createY3(nux, nuy, omega, isBent, Y3, M_fft, P_fft,
            Pm1_fft, S, G, TTO);
        Y4 = crtr->createY4(nux, nuy, omega, isBent, Y4, M_fft, S, F, G, TTO);

        F.kill();
    } else {
        // Depending on the symmetry, the Toeplitz developements will be
        // different to store the information about symmetry in variables

        // The difference between symmetry and anti-symmetry is taken into
        // account by a sign - taken into account in factor

        factor = 1;
        symmetry = none;
        xsymmetry = false;
        ysymmetry = false;
        if (symy == no_symmetry) {
            if (symx == symmetric || symx == anti_symmetric ) {
                symmetry = X;
                xsymmetry =true;
                if (symx == anti_symmetric)
                    factor = -1;
            }
        }
        if (symx == no_symmetry) {
            if (symy == symmetric || symy == anti_symmetric ) {
                ysymmetry =true;
                symmetry = Y;
                if (symy ==anti_symmetric)
                    factor = -1;
            }
        }
        if (symx == symmetric && symy == anti_symmetric){
            ysymmetry =true;
            xsymmetry =true;
            symmetry = both;
            factor = 1;
        }
        if (symx == anti_symmetric && symy == symmetric){
            ysymmetry =true;
            xsymmetry =true;
            symmetry = both;
            factor = -1;
        }
        // in other cases, no-symmetry will be considered

        switch (symmetry){
        case none:
            MfactorX = 0;
            MfactorY = 0;

            NfactorX = 0;
            NfactorY = 0;

            OfactorX = 0;
            OfactorY = 0;

            PfactorX = 0;
            PfactorY = 0;

            QfactorX = 0;
            QfactorY = 0;

            RfactorX = 0;
            RfactorY = 0;

            break;
        case X:
            MfactorX = -factor;
            MfactorY = 0;

            NfactorX = -factor;
            NfactorY = 0;

            OfactorX = -factor;
            OfactorY = 0;

            PfactorX = -factor;
            PfactorY = 0;

            QfactorX = factor;
            QfactorY = 0;

            RfactorX = factor;
            RfactorY = 0;

            break;
        case Y:
            MfactorX = 0;
            MfactorY = factor;

            NfactorX = 0;
            NfactorY = factor;

            OfactorX = 0;
            OfactorY = -factor;

            PfactorX = 0;
            PfactorY = factor;

            QfactorX = 0;
            QfactorY = -factor;

            RfactorX = 0;
            RfactorY = factor;

            break;
        case both:
            MfactorX = factor;
            MfactorY = factor;

            NfactorX = -factor;
            NfactorY = -factor;

            OfactorX = -factor;
            OfactorY = factor;

            PfactorX = -factor;
            PfactorY = -factor;

            QfactorX = factor;
            QfactorY = factor;

            RfactorX = factor;
            RfactorY = -factor;

            break;
        }

        // TODO In case of symmetry, the input of the matrix developpement is a
        // little bit different instead of having the Fourier
        //  transform as input, we input the toeplitz matrix

        // I think that this way of computing must be used for all matrices
        //  developpement for more compatibility

        // computation of the toeplitz matrices comming from the
        // fourier transform of the structure
        db_matrix N_toeplitz;
        db_matrix M_toeplitz;
        db_matrix P_toeplitz;
        db_matrix Pm1_toeplitz;
        db_matrix Q_toeplitz;
        db_matrix Qm1_toeplitz;

        // compute the derivative matrices that should be included as input
        crtr->Kr = M_fft.toeplitz_deriv_sym(xsymmetry,ysymmetry,nux,nuy,1,0);
        crtr->Kz = M_fft.toeplitz_deriv_sym(xsymmetry,ysymmetry,nux,nuy,0,1);

        //Compute X2


        //TTR = R_fft.toeplitz_mod(nux, nuy, 0,0,0,0).invert();
        TTR = R_fft.toeplitz_sym(xsymmetry, ysymmetry, RfactorX,RfactorY,
            false,false).invert();

        N_toeplitz = N_fft.toeplitz_sym(xsymmetry, ysymmetry, NfactorX,
            NfactorY, false,false);
        M_toeplitz = M_fft.toeplitz_sym(xsymmetry, ysymmetry, MfactorX,
            MfactorY, false,false);


        // In case of symmetry, S,F,G has not be calculated correctly yet!
        // To save memory space, we create the S matrix only if needed.
        if (isBent) {
            //S = bending(radius).toeplitz_mod(nux, nuy, 0,0,0,0);
            S = bending(radius).toeplitz_sym(xsymmetry, ysymmetry, -1,-1,
                false,false);
        } else {
            S = db_matrix::createUnitMatrix(TTR.getNrow(), TTR.getNcol());
        }

        //db_matrix F=transfPML(transf_x).toeplitz_mod(nux, nuy, 0,0,0,0);
        db_matrix F=transfPML(transf_x).toeplitz_sym(xsymmetry, ysymmetry,
             +1,+1, false,false);

        waitSemaphoreIO();
        cout << "Created ";
        postSemaphoreIO();


        X2 = crtr->createX2(nux, nuy, omega, isBent, X2, M_fft,
                    N_toeplitz, S, F, TTR);
        N_toeplitz.kill();

        //Compute X1 and X4
        //db_matrix G=transfPML(transf_y).toeplitz_mod(nux, nuy, 0,0,0,0);
        db_matrix G=transfPML(transf_y).toeplitz_sym(xsymmetry, ysymmetry,
            +1,+1, false,false);
        X4 = crtr->createX4(nux, nuy, omega, isBent, X4, M_fft, S, F, G, TTR);
        X1 = crtr->createX1(nux, nuy, omega, isBent, X1, M_fft, S, F, G, TTR);

        //Compute X3
        X3 = crtr->createX3(nux, nuy, omega, isBent, X3, M_fft,
            M_toeplitz, S, G, TTR);
        M_toeplitz.kill();
        TTR.kill();


        //Compute Y2
        Q_toeplitz = Q_fft.toeplitz_sym(xsymmetry, ysymmetry, QfactorX,
            QfactorY, false,false);
        Qm1_toeplitz = Qm1_fft.toeplitz_sym(xsymmetry, ysymmetry,
            QfactorX, QfactorY, false,false);
        TTO=O_fft.toeplitz_sym(xsymmetry, ysymmetry, OfactorX, OfactorY,
            false,false).invert();
        // In the case of normal field developpement, normal field toeplitz
        // matrices are computed. There are two types of matrices, minus or
        // plus toeplitz matrices.
        if (crtr == &NormalFieldSym ) {

            NormalFieldSym.Nyy =
            NormalFieldSym.Nyy_fft.toeplitz_sym(xsymmetry, ysymmetry,
                QfactorX,QfactorY, false,false);

            NormalFieldSym.Nxy =
            NormalFieldSym.Nxy_fft.toeplitz_sym(xsymmetry, ysymmetry,
                -QfactorX,-QfactorY, true,true);

        }
        Y2 = crtr->createY2(nux, nuy, omega, isBent, Y2, M_fft,
            Q_toeplitz, Qm1_toeplitz, S,F, TTO);

        Q_toeplitz.kill();
        Qm1_toeplitz.kill();

        //Compute Y1
        Y1 = crtr->createY1(nux, nuy, omega, isBent, Y1, M_fft, S, F, G, TTO);

        //Compute Y3
        if (crtr == &NormalFieldSym ) {

            NormalFieldSym.Nxx =
            NormalFieldSym.Nxx_fft.toeplitz_sym(xsymmetry, ysymmetry,
                PfactorX,PfactorX, false,false);

            NormalFieldSym.Nxy =
            NormalFieldSym.Nxy_fft.toeplitz_sym(xsymmetry, ysymmetry,
                -PfactorX,-PfactorX, true,true);
        }

        P_toeplitz = P_fft.toeplitz_sym(xsymmetry, ysymmetry, PfactorX,
            PfactorY, false, false);
        Pm1_toeplitz = Pm1_fft.toeplitz_sym(xsymmetry, ysymmetry,
            PfactorX, PfactorY, false,false);

        Y3 = crtr->createY3(nux, nuy, omega, isBent, Y3, M_fft,
                P_toeplitz, Pm1_toeplitz, S, G, TTO);

        P_toeplitz.kill();
        Pm1_toeplitz.kill();

        //Compute Y4
        Y4 = crtr->createY4(nux, nuy, omega, isBent, Y4, M_fft, S, F, G, TTO);

        // Free space:
        if (crtr == &NormalFieldSym ) {
            NormalFieldSym.Nxx.kill();
            NormalFieldSym.Nxy.kill();
            NormalFieldSym.Nyy.kill();
        }

        crtr->Kr.kill();
        crtr->Kz.kill();
        F.kill();
    }
//FIN

    TTO.kill();
    if(isBent && radius!=0) {
        X1 *= (1.0/radius);
        X2 *= (1.0/radius);
        X3 *= (1.0/radius);
        X4 *= (1.0/radius);
        Y1 *= (1.0/radius);
        Y2 *= (1.0/radius);
        Y3 *= (1.0/radius);
        Y4 *= (1.0/radius);
    }
    /* We assemble the blocks we just calculated, in order to obtain a second
       order problem in which the matrix involved in calculations is 4 times
       smaller than the original A matrix. Only (for example) the electric
       field is calculated.
       This is allowed by the presence of blocks full of zeros in the A matrix.

        B = [X1*Y1+X2*Y3    X1*Y2+X2*Y4] = [DD  EE]
            [X3*Y1+X4*Y3    X3*Y2+X4*Y4]   [FF  GG]
    */


    db_matrix DD =X2;
    DD*=Y3;
    DD=multsum(X1,Y1,DD);
    db_matrix EE =X2;
    EE*=Y4;
    EE=multsum(X1,Y2,EE);
    // We can suppress the unused matrices in order to free space if they are
    // not needed at a later moment

    if (!father->calcPropagation) X1.kill();
    if (!father->calcPropagation) X2.kill();
    db_matrix FF=X4;
    FF*=Y3;
    FF=multsum(X3,Y1,FF);

    if (!father->HxField) Y1.kill();
    if (!father->HyField) Y3.kill();

    db_matrix GG =X4;
    GG*=Y4;
    GG=multsum(X3,Y2,GG);
    if (!father->calcPropagation) X3.kill();
    if (!father->calcPropagation) X4.kill();

    if (!father->HxField) Y2.kill();
    if (!father->HyField) Y4.kill();

    waitSemaphoreIO();
    cout<<"DD, EE, FF, GG \n";
    postSemaphoreIO();
    cout.flush();

    B=db_matrix::mergeMatrixQuad(DD,EE,FF,GG);
    DD.kill();
    EE.kill();
    FF.kill();
    GG.kill();
    return B;
}

// Calculate the Fourier coefficients of the radius matrix
// in order to take into account the bending
db_matrix section::bending(double r)
{
    int ii, jj;
    int i, j;

    int Nx;                 // Size of the resulting matrix
    int Ny;

    Nx=father->dimx;
    Ny=father->dimy;

    db_matrix A(Ny, Nx);

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

            // i and j represent here the y and x harmonic number

            A(ii,jj)=(i==0?1.0:0.0)*
                (j==0?r:complex<double>(0,1)*father->tot_x*exp(complex<double>
                    (0,0*M_PI*j))/
                (2.0*M_PI*j));

       }
    }
    // 3.0/5.0*
    // //exp(complex<double>(0,2.0*M_PI*j*r/father->tot_x))
    return A;
}


/** Calculate the Fourier coefficients of the derivative of the coordinate
    transform function, as described in [2]
*/
db_matrix section::transfPML(
    transf_types type               // x or y type
    )
{
    int ii, jj;
    int i, j;

    double d;               // period (size) of the window
    double q;               // size of the mapped region
    int Nx;                 // Size of the resulting matrix
    int Ny;

    int k;
    bool tozero=false;

    Nx=father->dimx;
    Ny=father->dimy;

    if(type==transf_x) {
        d=father->tot_x;
        q=qx;
    } else {
        d=father->tot_y;
        q=qy;
    }

    db_matrix A(Ny, Nx);

    for(ii=0;ii<Ny;++ii) {
        for(jj=0;jj<Nx;++jj) {


            // Harmonics are stored in the standard FFT way inside of the
            // matrix being created.
            if ((ii+1)<Ny/2.0+1)
               i=ii;
            else
               i=-(Ny-ii);


            if ((jj+1)<Nx/2.0+1)
               j=jj;
            else
               j=-(Nx-jj);

            // i and j represent here the y and x harmonic number

            if (type==transf_x){  // x case
                k=j;
                tozero=i!=0;
            } else {
                k=i;
                tozero=j!=0;
            }

            // Here k represents the x or y harmonic number, depending on
            // which case we are. The variable tozero is true if the harmonic
            // number is above 1 for the axis which is NOT the one we
            // are fabricating the PML

            if(tozero) {
                A(ii,jj)=0;
            } else {
                A(ii,jj)=(k==0?1.0:0.0)-q/2.0/d*
                    ((1.0+gamma/4.0)*sinc(M_PI*(k*q/d))+
                    0.5*sinc(M_PI*(k*q/d-1.0))+
                    0.5*sinc(M_PI*(k*q/d+1.0))-
                    gamma/8.0*sinc(M_PI*(k*q/d-2.0))-
                    gamma/8.0*sinc(M_PI*(k*q/d+2.0)));
            }
       }
    }

    return A;
}


/** Add anisotropic PML at the borders of the calculation window. Type 1
*/
void section::PML_definition_rect1(double ep_x,
    double ep_y,
    double pos_x,
    double pos_y,
    complex<double> alpha)
{
    int dimx=father->dimx;
    int dimy=father->dimy;


    M_fft+=(MU_0/alpha - MU_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // x
    Mm1_fft+=(1.0/(MU_0/alpha)-1.0/MU_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // x

    N_fft+=(MU_0*alpha - MU_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // y
    Nm1_fft+=(1.0/(MU_0*alpha)-1.0/MU_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // x

    O_fft+=(MU_0*alpha - MU_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // z
    Om1_fft+=(1.0/(MU_0*alpha)-1.0/MU_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // x


    P_fft+=(EPS_0/alpha-EPS_0*n_0*n_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // x
    Pm1_fft+=(1.0/(EPS_0/alpha)-1.0/(EPS_0*n_0*n_0))* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // x

    Q_fft+=(EPS_0*alpha-EPS_0*n_0*n_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // y
    Qm1_fft+=(1.0/(EPS_0*alpha)-1.0/(EPS_0*n_0*n_0))* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // y

    R_fft+=(EPS_0*alpha-EPS_0*n_0*n_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // z
    Rm1_fft+=(1.0/(EPS_0*n_0*n_0*alpha)-1.0/(EPS_0*n_0*n_0))* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // z

}

/** Add anisotropic PML at the borders of the calculation window. Type 2
*/
void section::PML_definition_rect2(double ep_x,
    double ep_y,
    double pos_x,
    double pos_y,
    complex<double> alpha)
{
    int dimx=father->dimx;
    int dimy=father->dimy;

    M_fft+=(MU_0*alpha - MU_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // x
    Mm1_fft+=(1.0/(MU_0*alpha) - 1.0/MU_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // x

    N_fft+=(MU_0/alpha - MU_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // y
    Nm1_fft+=(1.0/(MU_0/alpha) - 1.0/MU_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // y

    O_fft+=(MU_0*alpha - MU_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // z
    Om1_fft+=(1.0/(MU_0*alpha)-1.0/MU_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // x


    P_fft+=(EPS_0*alpha-EPS_0*n_0*n_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // x
    Pm1_fft+=(1.0/(EPS_0*alpha)-1.0/(EPS_0*n_0*n_0))* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // x

    Q_fft+=(EPS_0/alpha-EPS_0*n_0*n_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // y
    Qm1_fft+=(1.0/(EPS_0/alpha)-1.0/(EPS_0*n_0*n_0))* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // y

    R_fft+=(EPS_0*alpha-EPS_0*n_0*n_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // z
    Rm1_fft+=(1.0/(EPS_0*alpha)-1.0/(EPS_0*n_0*n_0))* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // z

}

// Add anisotropic PML at the borders of the calculation window. Corners
void section::PML_definition_rectC(double ep_x,
    double ep_y,
    double pos_x,
    double pos_y,
    complex<double> alpha)
{
    int dimx=father->dimx;
    int dimy=father->dimy;
    O_fft+=(MU_0*alpha*alpha- MU_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // z
    Om1_fft+=(1.0/MU_0/alpha/alpha-1.0/MU_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // z

    R_fft+=(EPS_0*alpha*alpha-EPS_0*n_0*n_0)* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // z
    Rm1_fft+=(1.0/(EPS_0*alpha*alpha)-1.0/(EPS_0*n_0*n_0))* db_matrix::rectangleFT(ep_x, ep_y, father->tot_x, father->tot_y, dimx, dimy, pos_x, pos_y); // z

}


// Add anisotropic PML at the borders of the calculation window
void section::PML_definition(double pmldx,
    double pmldy,
    complex<double> alpha)
{
    // PML DEFINITION versus x
    PML_definition_rect1(pmldx, father->tot_y-pmldy*2, pmldx/2, father->tot_y/2,alpha);

    PML_definition_rect1(pmldx, father->tot_y-pmldy*2, father->tot_x-pmldx/2, father->tot_y/2, alpha);

    // PML DEFINITION versus y

    PML_definition_rect2(father->tot_x-pmldx*2, pmldy,  father->tot_x/2, pmldy/2, alpha);

    PML_definition_rect2(father->tot_x-pmldx*2, pmldy,  father->tot_x/2, father->tot_y-pmldy/2, alpha);

    // Corners
    PML_definition_rectC(pmldx, pmldy, pmldx/2, pmldy/2, alpha);

    PML_definition_rectC(pmldx, pmldy, pmldx/2, father->tot_y-pmldy/2, alpha);

    PML_definition_rectC(pmldx, pmldy, father->tot_x-pmldx/2, pmldy/2, alpha);

    PML_definition_rectC(pmldx, pmldy, father->tot_x-pmldx/2, father->tot_y-pmldy/2, alpha);

}

/* Assemble the magnetic field eigenvector matrix.
   Matrices X1, X2, X3 and X4 are needed, as well as W.
   The magnetic field eigenvector matrix is needed when calculating the
   S-matrix terms for each interface.
*/
db_matrix section::assembleVmatrix_alt(void)
{
    db_matrix V=db_matrix::mergeMatrixQuad(X1,X2,X3,X4).invert();
    V*=W;
    complex<double> e;
    int i,j;
    db_matrix Gamma(W.getNrow(), W.getNcol());

    for(i=0; i<W.getNcol(); ++i) {
        e=B(i,i);
        Gamma(i, i)= e;
    }
    V *=Gamma;

    return V;
}

/* Assemble the magnetic field eigenvector matrix.
   Matrices Y1, Y2, Y3 and Y4 are needed, as well as W.
   The magnetic field eigenvector matrix is needed when calculating the
   S-matrix terms for each interface.
   This is a faster algorithm.
*/
db_matrix section::assembleVmatrix(void)
{

    if((X1.getNcol() == 0 && X1.getNrow() == 0) ||
       (X2.getNcol() == 0 && X2.getNrow() == 0) ||
       (X3.getNcol() == 0 && X3.getNrow() == 0) ||
       (X4.getNcol() == 0 && X4.getNrow() == 0)) {
        throw parsefile_commandError("section::assembleVmatrix: At least one "
            "of the X matrices is undefined!");
    }

    db_matrix V=db_matrix::mergeMatrixQuad(X1,X2,X3,X4);

    V.invert();
    V*=W;
    complex<double> e;
    int i,j;
    /* Instead of creating a diagonal matrix of the propagation constants
       and multiplying it times XX, we do the calculations on the fly:
       we use much less memory and it is faster.
    */

    for(i=0; i<W.getNcol(); ++i) {
        e=B(i,i);
        for(j=0; j<W.getNrow(); ++j) {
            V(j, i) *= e;
        }
    }

    return V;
}

// Calculate the propagation matrix for a given section
db_matrix section::propagationMatrix(void)
{
    db_matrix P(B.getNrow(), B.getNcol());
    complex<double>e;

    for (int jj=0; jj<B.getNrow();++jj) {

        // At first, we obtain the propagation constants associated
        // to each eigenmode.

        e=B(jj,jj);

        P(jj,jj)=exp(complex<double>(0,-1)*e*tot_z);
    }

    return P;
}

/** Create the excitation field which can possibly be used for the current
    section.

    @parameter type a constant which determines which kind of excitation should
        be used.
    @paraleter real_neff if the excitation is of type SELECT_REAL_INDEX, this
        value is used for the effective index search: select the mode which has
        the real part of the refractive index which is the closest possible
        to the given value of real_neff.
        If the excitation is of type SELECT_MODE_POSITION, then this parameter
        is converted to integer and is interpreted as the modal position in the
        matrices.

    @parameter coeff the coefficient to be assigned to the excitation of the
        mode being selected.

    @parameter index_yz If the excitation is CONSTANT_FIELD_X or
        CONSTANT_FIELD_Y, this parameter specifies the index for the angle
        with whom the field is propagating with the normal to the interfaces.

    @parameter index_xz If the excitation is CONSTANT_FIELD_X or
        CONSTANT_FIELD_Y, this parameter specifies the index for the angle
        with whom the field is propagating with the normal to the interfaces.

    @return a matrix containing the modal excitation (in the modal base).

*/
db_matrix section::create_excitation(int type, double real_neff,
    complex<double> coeff, int index_yz, int index_xz)
{
    db_matrix excitation(W.getNrow(),1);

    cout << "Creating the excitation field.\n";
    // Select the type of field excitation to be used.

    int i;
    int i_min=-1;
    complex<double> e;
    complex<double> g;

    double im_max=-numeric_limits<double>::max();
    int i_max=-1;
    int shift=0;
    int pos=0;

    switch(type) {
        case SELECT_REAL_INDEX:
            // In this case, we basically search the eigenmode which has the
            // real part of the effective index which is the closest possible
            // to the given value.

            cout <<"Search for the mode with the real part of n_eff as close"
                " as possible to n = "<< real_neff<<".\n";

            // We search for the eigenmode having the highest real part of
            // the effective index, among those which are declared as
            // "interesting" from the previous configuration of the
            // system.

            i_min = select_real(real_neff);
            cout <<"Mode position: "<<i_min<<", coefficient "<<coeff<<endl;
            // Here the mode has been selected.
            if(i_min>0) {
                excitation(i_min,0)=coeff;
                cout <<"Excitation mode n_eff="<<
                    father->getEffectiveIndex(B(i_min,i_min))<<"\n";
            }
            break;

        case SELECT_MODE_POSITION:
            pos=(int)real_neff;
            cout <<"Selected the mode position: "<<pos<<", coefficient "
                <<coeff<<endl;
            if(pos>=0)
                excitation(pos,0)=coeff;
            else
                parsefile_commandError("Position out of boundaries.");
            break;

        case SELECT_FUNDAMENTAL:
            // Search for the fundamental mode.

            // We search for the eigenmode having the highest real part of
            // the effective index, among those which are declared as
            // "interesting" from the previous configuration of the
            // system.

            for(i=0;i<B.getNrow();++i) {
                e=B(i,i);

                if(isInteresting(e)) {
                    g = father->getEffectiveIndex(e);
                    if (g.real() > im_max){
                        im_max = father->getEffectiveIndex(e).real();
                        i_max=i;
                    }
                }
            }
            if(i_max>0) {
                excitation(i_max,0)=coeff;
                cout <<"Excitation mode n_eff="<<
                    father->getEffectiveIndex(B(i_max,i_max))<<"\n";
            }
            break;
        default:
            throw parsefile_commandError(
                "Ups... section::create_excitation is to be checked!\n");
    }

    return excitation;
}

db_matrix section::create_excitation_const(bool isY,
    double real_neff,
    complex<double> coeff , int index_yz, int index_xz)
{
    db_matrix excitation(W.getNrow(),1);

    cout << "Creating the excitation field.\n";
    // Select the type of field excitation to be used.

    int i;
    int i_min=-1;
    complex<double> e;
    complex<double> g;

    double im_max=-numeric_limits<double>::max();
    int i_max=-1;
    int shift=0;
    if(isY)
        shift = W.getNrow()/2;

    int snux=(father->dimx-1)/2+1;
    int snuy=(father->dimy-1)/2+1;

    int mm, ll;

    db_matrix rd_fft(snuy, snux);

    if(index_yz==0 && index_xz==0) {
        cout <<"Constant excitation "<< coeff<<".\n";
        rd_fft(0,0) = coeff;
    } else {
        double angleyz = asin(index_yz*father->lambda/
                (n_0.real()*father->tot_y));
        double anglexz = asin(index_xz*father->lambda/
                (n_0.real()*father->tot_x));
        double results[2];
        results[0]=anglexz/M_PI*180.0;
        results[1]=angleyz/M_PI*180.0;
        father->nP->insertArray("ans", 2, results);
        cout <<"Excitation angles "<< anglexz/M_PI*180.0<<", ";
        cout << angleyz/M_PI*180.0  <<".\n";
        cout <<"Coefficient "<< coeff<<".\n";

        // Prevent divide by zero.
        double pitchy=0;
        if(index_yz!=0)
            pitchy=father->tot_y/index_yz;

        double pitchx=0;
        if(index_xz!=0)
            pitchx=father->tot_x/index_xz;

        cout << "pitchx="<<pitchx<<"  pitchy="<<pitchy<<endl;
        double x,y;
        db_matrix rd(snuy, snux);
        for(int l=0;l<snuy; ++l) {
            y=(l-(snuy)/2.0)*father->tot_y/snuy;
            for(int m=0; m<snux; ++m) {
                x=(m-(snux)/2.0)*father->tot_x/snux;
                rd(l,m)=coeff;

                if(pitchx!=0)
                    rd(l,m)*=exp(complex<double>(0,-1)
                        *(2.0*M_PI*x/pitchx));
                if(pitchy!=0)
                    rd(l,m)*=exp(complex<double>(0,-1)
                        *(2.0*M_PI*y/pitchy));

                rd(l,m)/=(double)(snux*snuy);
            }
        }

        rd_fft = rd.ifft2();
    }
    //db_matrix N(W.getNrow(),1);

    // The read matrix has been constructed according to standard
    // rules followed by FFT. We need to unroll it.
/*
    for(int l=0;l<snuy; ++l) {
        for(int m=0; m<snux; ++m) {
            mm=m+(snux/2);
            if (mm>=snux)
                mm=mm-snux;

            ll=l+(snuy/2);
            if (ll>=snuy)
                ll=ll-snuy;

            N(mm+ll*snux+shift,0) = rd_fft(l,m);

        }
    }*/

    // Enroll the vector depending on the symmetry
    bool xsym = false;
    bool ysym = false;
    if (father->symx == anti_symmetric || father->symx == symmetric) {
        xsym = true;
        snux = snux/2 +1;
    }
    if (father->symy == anti_symmetric || father->symy == symmetric) {
        ysym = true;
        snuy = snuy/2 +1;
    }
    if (isY)
        shift = snux*snuy;

    db_matrix N = rd_fft.fft2vector(xsym,ysym,shift);

    db_matrix Wm1 = W;
    Wm1.invert();

    excitation = Wm1*N;

    return excitation;
}

/** Create the excitation field which can possibly be used for the current
    section using a input file.

    @input is the matrix read from file
    @return a matrix containing the modal excitation (in the modal base).
*/
db_matrix section::create_excitation_from_file(bool excitation_fy,
    db_matrix &rd )
{
    db_matrix excitation(W.getNrow(),1);

    // rd contains the excitation given by the file. We should
    // take its Fourier transform and then switch to the modal
    // representation.

    int shift = 0;

    int nx = rd.getNcol();
    int ny = rd.getNrow();

    rd *= 1.0/nx/ny;

    int snux=(father->dimx-1)/2+1;
    int snuy=(father->dimy-1)/2+1;

    cout << "The wanted size is " << snux << " x " << snuy <<
        " points.\n";
    cout.flush();
    db_matrix rd_fft = rd.ifft2().zero_pad(snuy,snux);



    /*if (excitation_fy)
        shift = snux*snuy;

    db_matrix N( W.getNrow(),1);

    int mm, ll;

    // The read matrix has been constructed according to standard
    // rules followed by FFT. We need to unroll it.

    for(int l=0;l<snuy; ++l) {
        for(int m=0; m<snux; ++m) {
            mm=m+(snux/2);
            if (mm>=snux)
                mm=mm-snux;

            ll=l+(snuy/2);
            if (ll>=snuy)
                ll=ll-snuy;

            N(mm+ll*snux+shift,0) = rd_fft(l,m);

        }
    }*/

//JEROME

    // Enroll the vector depending on the symmetry
    bool xsym = false;
    bool ysym = false;
    if (father->symx == anti_symmetric || father->symx == symmetric) {
        xsym = true;
        snux = snux/2 +1;
    }
    if (father->symy == anti_symmetric || father->symy == symmetric) {
        ysym = true;
        snuy = snuy/2 +1;
    }
    if (excitation_fy)
        shift = snux*snuy;

    db_matrix N = rd_fft.fft2vector(xsym,ysym,shift);

    // N now contains the excitation vector represented in the Fourier
    // domain, with a number of harmonics equal to the one defined at
    // the beginning of the calculation. We need to represent it in
    // the modal space in order to be used for calculation.

    db_matrix Wm1 = W;
    Wm1.invert();

    excitation = Wm1*N;
    return excitation;

}
/** Selects the mode having the real part of the effective index closest to a
    given value

    @param real_neff the wanted real part of the refractive index
    @return the index of the mode

*/
int section::select_real(double real_neff)
{
    double d_min=numeric_limits<double>::max();
    int i_min = -1;
    complex<double> g;

    for(int i=0;i<B.getNrow();++i) {
        g = father->getEffectiveIndex(B(i,i));

        if (abs(g.real()-real_neff) < d_min){
            d_min = abs(g.real()-real_neff);
            i_min=i;
        }
    }

    return i_min;
}

/** Calculates all modes of the interesting section.
*/
void section::modes_calculation()
{
    // Create the W matrix from the structure we defined
    waitSemaphoreIO();
    cout << "Creating the electric field propagation matrix\n";
    postSemaphoreIO();
    calcWmatrix(father->getOmega());

    // Calculate eigenvectors
    waitSemaphoreIO();
    cout<<"Calculating eigenvalues and eigenvectors\n";
    postSemaphoreIO();
    cout.flush();
    B.eig(&W);

    // The calculation of eigenvalues and eigenvectors in this form gives us
    // eigenvalues which are the square of the propagation constants.
    // We memorise in the principal diagonal of B their square roots (thus the
    // propagation constants), by changing appropriately the sign of their
    // imaginary part if we want to force convergent results.

    complex<double>e, g;

    for (int j=0; j<B.getNrow();++j) {
        e=B(j,j);
        g = sqrt(e);

        if(father->ensureConvergence && g.imag()>0) {
            // here the calculated square root can have both the real part
            // as well as the imaginary part negative.
            // If the imaginary part is positive, this means that the solution
            // might be associated with an exploding exponential where
            // the propagation will take place.
            if(g.real()>0) {
                g = complex<double>(g.real(), -g.imag());
            } else {
                g = complex<double>(-g.real(), -g.imag());
            }
        }
        e=B(j,j)= g;
    }
}

/* Find interesting modes and print their effective index.
    To be called after section::mode_calculation().
    The result will be a vector containing all effective indices, with
    alternating real part/imaginary part for all of them.
    The size of the vector will therefore be twice the number of
    "interesting" modes.
*/
vector<double> section::find_interesting()
{
    vector<double> results;
    vector<int> idx;

    find_interesting_g(results, idx);

    return results;
}

/** Calculates both the values of interesting indices
    and their indices.

*/
void section::find_interesting_g(vector<double> &results, vector<int> &idx)
{
    bool interesting=false;
    complex<double>e;
    complex<double>g;
    double order;
    double quality;

    int k=0;
    // Find interesting modes and print their effective index

    if(useAzimuthalOrder) {
        cout<<"Searching for interesting modes: " <<
            minAzimuthalOrder << " <= order <= " <<
            maxAzimuthalOrder<<"\n";
    } else {
        cout<<"Searching for interesting modes: "<<n_0.real()<<
                " < Re{n_eff} < "
                << n_g.real()<<"\n";
        cout<<"                            "<<low_imag<<
                " < Im{n_eff} < "
                << high_imag<<"\n";
    }
    cout.flush();

    for (int j=0; j<B.getNrow();++j) {
        e=B(j,j);

        interesting = isInteresting(e);
        order = father->getOrder(e, radius);
        quality = father->getQuality(e);
        g = father->getEffectiveIndex(e);
        if(interesting) {
            results.insert(results.end(),g.real());
            results.insert(results.end(),g.imag());

            idx.insert(idx.end(),j);

            cout <<"Interesting mode found "<<k++<<", effective index: "
                <<g;
            if(isBent)
                cout<<" azimuthal order: "
                    << order
                    <<" Q = "<< quality
                    <<" ("<<father->getLosses_dB_cm(e)<<" dB/cm)";
            cout<<"\n";
            cout.flush();
        }
    }
    return;
}


/** Determine wether a given mode is "interesting" or not.

*/
bool section::isInteresting(complex<double> e)
{
    complex<double> g;
    double order, quality;
    bool interesting=false;


    g = father->getEffectiveIndex(e);
    order = father->getOrder(e, radius);
    quality = father->getQuality(e);


    /* We need to decide wheter a given mode is interesting enough to
            retain its characteristics. */
    if(useAzimuthalOrder) {
        /* The first possibility (useful for bent waveguides) is to
            take the azimuthal order and see if it is inside the defined
            window */

        interesting= order >= minAzimuthalOrder &&
                     order <= maxAzimuthalOrder;
    } else {
        /* If not, we will do the standard way, by taking the effective
            index and comparing it to the window defined by the minimum
            and maximum refractive indices */

        interesting=
                g.real()> n_0.real() &&
                g.real()< n_g.real() &&
                g.imag()< high_imag &&
                g.imag()> low_imag;
        }
    return interesting;
}

complex<double> section::expectedRefractiveIndex(double px, double py)
{
    complex<double> n=n_0;

    // At first, we check if there is a background index distribution
    // loaded in memory. This happens when the user adopts the INDFILE
    // command.

    if(background.getNrow()>0) {
        int posx = round((px/father->tot_x+0.5)*background.getNcol());
        int posy = round((py/father->tot_y+0.5)*background.getNrow());

        if (posx<0 || posy<0 ||
            posx>=background.getNcol() || posy>=background.getNrow()) {
            return -1;
        }
        n = background(posy, posx);
    }

    // Then, we superpose the rectangles eventually present in the section.
    list<rect_wg>::iterator i;

    for(i=rects.begin(); i != rects.end(); ++i) {
        if (i->isInside(px, py)) {
            n += i->n - n_0;
        }
    }

    return n;
}
