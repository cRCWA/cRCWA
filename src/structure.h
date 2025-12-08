/***************************************************************************
*     CLASS structure                                                      *
*     Davide Bucci, CROMA                                                  *
*     Jérôme Michallon, CROMA                                              *
*     MINATEC-INPG, 3, parvis Luis Neel                                    *
*     38016, Grenoble CEDEX, France                                        *
*                                                                          *
*     bucci@minatec.inpg.fr                                                *
*                                                                          *
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
#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <cstdlib>
#include "parsefile.h"
#include "block_matrix.h"
#include "section.h"
#include "commands.h"

#define ERROR 1
#define EPS_0 8.8541878176e-12
#define MU_0 (4e-7*M_PI)
#define C_VACUUM 299792458.0
#define INIT_SIZE 10

typedef enum field_e_t{Ex,Ey, Ez, Hx, Hy, Hz, Dx, Dy, Dz} field_e;

/** Class describing the additional output data to be calculated during the
    propagation

*/
class outputdata {
public:
    // if type t ==i
    bool should_record_integral;
    field_e window_field_type;
    double window_centerx;
    double window_centery;
    double window_width;
    double window_height;
    string window_file_name;

    outputdata ()
    {
        should_record_integral = false;
        should_record_generation_rate = false;
        window_field_type = Ex;
        window_centerx = 0;
        window_centery = 0;
        window_width = 0;
        window_height = 0;
        generation_radial_step_size = 0;
        generation_from_z0 = 0;
        generation_to_z1 = 0;
        generation_to_r1 = 0;
        generation_absorptance_power_interpolated = 0;
        generation_real_from_z0 = 0.0;
        generation_real_to_z1 = 0.0;
        generation_nz = 0;
    }

    // type t ==g

    bool should_record_generation_rate;
    double generation_radial_step_size;
    double generation_from_z0;
    double generation_to_z1;
    double generation_real_from_z0;
    double generation_real_to_z1;
    double generation_to_r1;
    double generation_absorptance_power_interpolated;
    int generation_nz;
    string generation_file_name;
};

class structure : public parsefile
{
private:
    // The size of the arrays containing the interfaces and the sections
    int actual_size;

    // The global S matrix
    finterface globalS;

    // The total number of sections contained in the structure
    int number_of_sections;

    // The pointer at the current section
    section *cur;
    // The pointer at the first section
    section *sec_list;
    // List of interfaces between adjacent sections
    finterface *interf_list;

    // Size of the FFT matrices
    /* NOTE. This takes into account the total number of positive and negative
            frequency contributions to be obtained from the FT. In other words,
            if the number of harmonics is n, the size of the matrix which
            should contain them is 2*n+1.
    */
    int dimx;
    int dimy;

    // Symmetry type:
    symmetry_e symx;
    symmetry_e symy;

    // Size of the calculation window
    double tot_x;
    double tot_y;

    // Additional output data to be recorded during the propagation
    outputdata additional_output_data;

    // Matrix containing the eigenvectors for Bloch modes
    db_matrix bloch_eigvect;
    // Matrix (column vector) containing eigenvalues for Bloch modes
    db_matrix bloch_eigval;

    // Specifies that the excitation vectors of each interface have already
    // been calculated and the propagation can be done immediately.
    bool excitationCalculated;

    // Wavelength
    double lambda;

    // Specifies which component of the H field will be calculated
    bool HxField;
    bool HyField;
    bool calcPropagation;

    // Force all modes to have a negative imaginary part of the
    // effective index
    bool ensureConvergence;

    static void fileoutputGP(db_matrix out, FILE *f, rimco_e rimc, double dx,
        double dy);

    static void fileoutputOW(db_matrix out, FILE *f, double dx, double dy);

    double getOrder(complex<double> e, double radius);
    double getQuality(complex<double> e);
    double getLosses_dB_cm(complex<double> e);


    // terms for taking in account the angles in the derivative
    // calculations and in the field output.
    double ksinthetax;
    double ksinthetay;

public:
    complex<double> getEffectiveIndex(complex<double> e);
    static db_matrix getPropagation(db_matrix &B, double tp,
        double tm,db_matrix &excitation_p, db_matrix &excitation_m,
        bool calcH);

    // Low-level implementation of the commands
    void do_size(double sx, double sy);
    void do_harmonics(int nx, int ny);
    void set_wavelength(double l);
    double get_wavelength(void);
    vector<double> do_solve(void);
    void do_assemble(void);
    void set_ensureConvergence(bool s);
    void do_wants(const char *code);
    void do_section(const double tl);
    void do_select(unsigned int number);
    double do_powerz(const double z);
    double do_monitor(const double z, const double wx, const double wy,
        const double px, const double py);

    vector<double> do_bloch(void);





    section *getCur(void) {return cur;}

    // Number of threads allowed during numerical-intensive,
    // possibly parallelizable procedures.
    int number_of_threads_allowed;

    void ensureSize(int newsize);
    double getOmega(void);

    bool symmetryNotTakenIntoAccount(void)
    {
        return (symx==previous_version)&&(symx==previous_version);
    }

    vector<double> calculateBlochModes();
    void outBlochFile(double pp, char *name);
    void outBlochMaxImag(double pp, double pa);
    complex<double> outBlochClosest(double pp, complex<double> pa);
    int indexBlochClosest(double pp, double pa);
    void setAngles(double n0, double thetax, double thetay);
    double getKsinthetax(void);
    double getKsinthetay(void);

    // Parallelized structures:
    static void *calc_eigenvalues(void *threadarg);



    double getTotalLength(void);

    // Standard constructor
    structure(numParser *p);
    // Destructor
    ~structure();

    // Reset the structure
    void reset(void);

    friend class commands;
    friend class section;

    // Parallel calculations callbacks
    static void *calcV(void *threadarg);

};

#endif
