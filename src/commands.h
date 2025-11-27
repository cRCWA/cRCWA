/***************************************************************************
*     CLASS commands                                                       *
*     Davide Bucci, CROMA     march 2008- 2014                             *
*     Jérôme Michallon, CROMA     2012                                     *
*     MINATEC-INPG, 3, parvis Luis Neel                                    *
*     38016, Grenoble CEDEX, France                                        *
*                                                                          *
*     bucci@minatec.inpg.fr                                                *
*                                                                          *
*     Version: 1.2                                                         *
****************************************************************************
*    C++ class containing all the commands recognized by the parser.       *
****************************************************************************/

// Prevent multiple includes
#ifndef COMMANDS_H
#define COMMANDS_H


#include "parsefile.h"
#include "block_matrix.h"
#include "section.h"
#include "structure.h"

#define ERROR 1

db_matrix fft2cent(int dimx, int dimy, db_matrix dd);

// Structure with data for the propagation calculation.
struct propagation_data{
    double dz;  // discretization step in the propagation axis.
    int sizex;  // number of points in x
    int sizey;  // number of points in y
    bool recordAdditionalData;  // Should record additional output data?
    enum rimco_e_t rimc;    // rimc type of output (real part,
                            // imaginary, module, complex)
    char *component;    // type of component to be calculated:
                        // {Ex|Ey|Ez|Hx|Hy|Hz|Dx|Dy|Dz|Ex2|Ey2|Ez2}
    char *fname;    // the file to be written
};

class commands:public parsefile
{
private:

static void process_section(double dx, double dy, double dz,
        section &c_section, double &alpha_t, double &z0, double &z,
        FILE *f, enum rimco_e_t rimc, double &x_t, double &z_t, int &nz,
        db_matrix &excitation_p, db_matrix &excitation_m,
        bool applyShift, bool calcH, int snux, int snuy, int sj, int si,
        FILE *f1, bool calcD, db_matrix &epsilonxy, db_matrix &epsz,
        db_matrix &muz, bool calcz,bool improve_representation, FILE *f2);

protected:
    void init(structure &s);

public:

    // Standard constructor
    commands(structure &s, numParser *p);

    /* Here follows the file reading commands */
    static int c_diesis(parsefile *obj, int argc,char *argv[]);
    static int c_substrate(parsefile *obj, int argc,char *argv[]);
    static int c_size(parsefile *obj, int argc,char *argv[]);
    static int c_harmonics(parsefile *obj, int argc,char *argv[]);
    static int c_wavelength(parsefile *obj, int argc,char *argv[]);
    static int c_rectangle(parsefile *obj, int argc,char *argv[]);
    static int c_pml(parsefile *obj, int argc,char *argv[]);
    static int c_solve(parsefile *obj, int argc,char *argv[]);
    static int c_inpstruct(parsefile *obj, int argc,char *argv[]);
    static int c_outgmodes(parsefile *obj, int argc,char *argv[]);
    static int c_highindex(parsefile *obj, int argc,char *argv[]);
    static int c_lowindex(parsefile *obj, int argc,char *argv[]);
    static int c_pml_transf(parsefile *obj, int argc,char *argv[]);
    static int c_indfile(parsefile *obj, int argc,char *argv[]);
    static int c_bend(parsefile *obj, int argc,char *argv[]);
    static int c_wants(parsefile *obj, int argc,char *argv[]);
    static int c_spectrum(parsefile *obj, int argc,char *argv[]);
    static int c_propagation(parsefile *obj, int argc,char *argv[]);
    static int c_section(parsefile *obj, int argc,char *argv[]);
    static int c_order(parsefile *obj, int argc,char *argv[]);
    static int c_assemble(parsefile *obj, int argc,char *argv[]);
    static int c_excitation(parsefile *obj, int argc,char *argv[]);
    static int c_quit(parsefile *obj, int argc,char *argv[]);
    static int c_help(parsefile *obj, int argc,char *argv[]);
    static int c_coefficient(parsefile *obj, int argc,char *argv[]);
    static int c_select(parsefile *obj, int argc,char *argv[]);
    static int c_carpet(parsefile *obj, int argc,char *argv[]);
    static int c_matdev(parsefile *obj, int argc,char *argv[]);
    static int c_clear(parsefile *obj, int argc,char *argv[]);
    static int c_norm(parsefile *obj, int argc,char *argv[]);
    static int c_outdata(parsefile *obj, int argc,char *argv[]);
    static int c_eigenvc(parsefile *obj, int argc,char *argv[]);
    static int c_eigenam(parsefile *obj, int argc,char *argv[]);
    static int c_memocc(parsefile *obj, int argc,char *argv[]);
    static int c_power(parsefile *obj, int argc,char *argv[]);
//JEROME
    static int c_powerZ(parsefile *obj, int argc,char *argv[]);
    static int c_monitor(parsefile *obj, int argc,char *argv[]);
    static int c_symmetry(parsefile *obj, int argc,char *argv[]);
// DB
    static int c_let(parsefile *obj, int argc,char *argv[]);
    static int c_bloch(parsefile *obj, int argc,char *argv[]);
    static int c_outbloch(parsefile *obj, int argc,char *argv[]);
    static int c_parallel(parsefile *obj, int argc,char *argv[]);
    static int c_selmodes(parsefile *obj, int argc,char *argv[]);
    static int c_modepos(parsefile *obj, int argc,char *argv[]);
    static int c_angles(parsefile *obj, int argc,char *argv[]);

    // Important utility methods used by the commands. (maybe make private?)

    static void getEpsilonAndMu(section &c_section, bool calcD, bool calcz,
        bool calcH, bool improve_representation, db_matrix &epsilonxy,
        db_matrix &epsz, db_matrix &muz);
    static db_matrix getFourierField(section &c_section,db_matrix &fields,
        bool calcH, bool calcD, db_matrix &epsilonxy,
        db_matrix &epsz,db_matrix &muz,bool calcz);

    static db_matrix getField(section &c_section,db_matrix &fields,
        bool applyShift, bool calcH, int snux, int snuy, int sj, int si,
        bool calcD, db_matrix &epsilonxy,
        db_matrix &epsz,db_matrix &muz,bool calcz,bool improve_representation,
        bool additionnalEy, db_matrix &outEy);
    static db_matrix getGeneration_rate(section &c_section,double tp, double tm,
        db_matrix &excitation_p, db_matrix &excitation_m,
        int snux, int snuy,double r1,double dr,
         db_matrix &epsilonxy,
        db_matrix &epsz,db_matrix &muz);

    static db_matrix read_file(char* filename, double &tot_x, double &tot_y,
        double mult);
    static structure * wantedParameters(char *cname,
        int argc, int w, parsefile *obj);

    static void outputfield(db_matrix &out, double dx, double dy,
    section &c_section, double &alpha, double &alpha_t, double &z0, double &z,
    FILE *f, enum rimco_e_t rimc, double &x_t, double &z_t, double &x_c,
    double &z_c, FILE *f1);



    // Parallelized structures
    static void *propagation_structure(void *threadarg);
};

#endif
