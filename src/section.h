/***************************************************************************
*     CLASS section                                                        *
*     Davide Bucci, CROMA     march 2010                                   *
*     Jérôme Michallon, CROMA     2012                                     *
*     MINATEC-INPG, 3, parvis Luis Neel                                    *
*     38016, Grenoble CEDEX, France                                        *
*                                                                          *
*     bucci@minatec.inpg.fr                                                *
*                                                                          *
*     Version: 1.2                                                         *
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

// Prevent multiple includes
#ifndef SECTION_H
#define SECTION_H

#include <list>
#include <vector>

#include "block_matrix.h"
#include "finterface.h"

#include "matrixDev.h"
#include "matrixdevsym.h"

typedef enum rimco_e_t{R,I,M,C,O,S} rimco_e;


#define ERROR 1

// This is needed since section contains a pointer to the structure class
// We do not need to define the class here, but just provide a placeholder
class structure;


class rect_wg
{
public:
    double cx;
    double cy;
    double w;
    double h;
    complex<double> n;

    // Check if the given point lies inside the rectangle
    bool isInside(double px, double py) {
        return (abs(px-cx)<=(w/2.0) && abs(py-cy)<=(h/2.0))?true:false;
    }
};

class section_data {
public:
    // Index database. Here we store the background from INDFILE command
    db_matrix background;
    // Here we store the data of the rectangles
    list<rect_wg> rects;

    // The matrix developements are dependent on each section.
    fPMLafterOPT PMLafterOPT;
    fPMLafter PMLafter;
    fPMLbefore PMLbefore;
    fNonDev NonDev;
    fNF NormalField;

    fNonDevSym NonDevSym;
    fNFsym NormalFieldSym;

    matrixDev *crtr;

    // Mu_x
    db_matrix M_fft;
    db_matrix Mm1_fft;

    // My_y
    db_matrix N_fft;
    db_matrix Nm1_fft;

    // Mu_z
    db_matrix O_fft;
    db_matrix Om1_fft;

    // Epsilon_x
    db_matrix P_fft;
    db_matrix Pm1_fft;

    // Epsilon_y
    db_matrix Q_fft;
    db_matrix Qm1_fft;

    // Epsilon_z
    db_matrix R_fft;
    db_matrix Rm1_fft;

    // Pointer to the structure containing this section
    structure *father;

    // Propagation constant matrix

    // We memorise in the principal diagonal of B the square roots of
    // the eigenvalues calculated (thus the
    // propagation constants), by changing appropriately the sign of their
    // imaginary part if we want to force convergent results.
    db_matrix B;

    // Electrical field eigenvectors matrix
    db_matrix W;

    // Magnetic field eigenvector matrix
    db_matrix V;

    // Counterpropagative excitation (at the right of this section)
    db_matrix sWm;

    // Propagative excitation (at the left of this section)
    db_matrix sWp;

    // These matrices will be used only if the H field is requested. See
    // the description of the WANTS command.

    db_matrix X1;
    db_matrix X2;
    db_matrix X3;
    db_matrix X4;

    db_matrix Y1;
    db_matrix Y2;
    db_matrix Y3;
    db_matrix Y4;

    // Total length of the section
    double tot_z;

    // Radiant frequency
    // double omega;

    // PML sizes
    double pmldx;
    double pmldy;

    // Substrate refractive index
    complex<double> n_0;

    // Alpha factor, for the PMLs
    complex<double> alpha;

    // Bend radius of the center of the waveguide
    double radius;

    bool isBent;

    // Index of the structure having the biggest real part of the r.i.
    complex<double> n_g;

    // Size of the mapped region when using PMLs described in [1]
    double qx;
    double qy;
    // Complex PML factor as described in [2]
    complex<double> gamma;

    double low_imag;
    double high_imag;

    // Search using the azimuthal order
    bool useAzimuthalOrder;
    double minAzimuthalOrder;
    double maxAzimuthalOrder;

    // transform type
    enum transf_types {transf_x, transf_y};

    // Is true if the substrate has been set
    bool isSubstrateSet;
};

class section:public section_data
{
private:
    // On the basis of the provided information, calculate the refractive
    // index of the given point, without passing by the Fourier transform.
    complex<double> expectedRefractiveIndex(double px, double py);


    void PML_definition_rect1(double ep_x,
        double ep_y,
        double pos_x,
        double pos_y,
        complex<double> alpha);

    void PML_definition_rect2(double ep_x,
        double ep_y,
        double pos_x,
        double pos_y,
        complex<double> alpha);
    void PML_definition_rectC(double ep_x,
        double ep_y,
        double pos_x,
        double pos_y,
        complex<double> alpha);

    db_matrix transfPML(transf_types type);

    db_matrix bending(double r);

    static double sinc(double alpha) {return (alpha==0)?1.0:sin(alpha)/alpha;}

public:

    // Standard constructor
    section();
    // Copy constructor
    section(section &m);
    // Copy operator
    section& operator= (const section &m);

    // Destructor
    ~section();

    // Resets the current section
    void reset(void);


    static double integralPoynting_rectangle(structure *p, db_matrix &V,
        db_matrix &W, int column, double wx, double wy, double px, double py);

    static double integralPoynting(structure *p, db_matrix &V, db_matrix &W,
        int column);

    void set_substrate(complex<double>subs);
    void add_rectangle(complex<double>n_g, double wx, double wy, double px,
        double py);
    void set_lowindex(complex<double>n_g);
    void set_highindex(complex<double>n_g);
    void set_pml(double wx, double wy, complex<double> alpha);
    void set_pml_transf(double qdx, double qdy, complex<double> g);
    void set_bend(double r);
    db_matrix do_inpstruct(int sj, int sx, char* otype);
    list<db_matrix> do_outgmodes(int sj, int si, char *ftype, 
        rimco_e rimc, const char *fname);
    void store_refractive_index(db_matrix rd);
    void do_order(double min, double max);
    db_matrix &getW(void){return W;}
    db_matrix &getB(void){return B;}


    bool isInteresting(complex<double> e);

    void modes_calculation();
    vector<double> find_interesting();
    void find_interesting_g(vector<double> &results, vector<int> &idx);


    void setFather(structure *f)  {father=f;}
    structure *getFather(void) {return father;}

    // Add anisotropic PML at the borders of the calculation window
    void PML_definition(double pmldx,
        double pmldy,
        complex<double> alpha);

    /* Assemble the W matrix which represent the propagation operator in the
        Fourier space.
    */
    db_matrix &calcWmatrix(double omega);

    db_matrix assembleVmatrix(void);
    db_matrix assembleVmatrix_alt(void);
    db_matrix propagationMatrix(void);

    const static int SELECT_FUNDAMENTAL = 0;
    const static int SELECT_REAL_INDEX = 1;
    const static int SELECT_MODE_POSITION = 2;

    db_matrix create_excitation(int type, double real_neff,
        complex<double> coeff, int index_yz, int index_xz);

    db_matrix create_excitation_const(bool isY,
        double real_neff, complex<double> coeff , int index_yz, int index_xz);

    db_matrix create_excitation_from_file(bool excitation_fy, db_matrix &rd );

    int select_real(double real_neff);

    friend class structure;
    friend class finterface;
    friend class commands;
};

#endif
