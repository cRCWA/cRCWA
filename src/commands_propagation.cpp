/***************************************************************************
*     CLASS commands                                                       *
*     Davide Bucci, CROMA     2008-present                                 *
*     Jérôme Michallon, CROMA     2013                                     *
*     MINATEC-Grenoble INP, 3, parvis Luis Neel                            *
*     38016, Grenoble CEDEX, France                                        *
*                                                                          *
*     bucci@minatec.grenoble-inp.fr                                        *
*                                                                          *
****************************************************************************
*    User commands available for the parsing.                              *
****************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstring>
#include <stdlib.h>
#include <stdio.h>

#include <unistd.h>

#include "block_matrix.h"
#include "structure.h"
#include "section.h"
#include "commands.h"
#include "parallelize.h"
#include "phys_constants.h"

/*  PROPAGATION: once all the propagation matrices have been assembled,
    propagate the excitation fields defined and write on a file the results.
    In order to calculate the propagation of the field, the user should enter
    the commands 'wants propagation' before 'solve' and then 'assemble'.

    Usage:
    propagation e rimc dz nx ny filename

    Parameters:
    e       must be {Ex|Ey|Ez|Ex2|Ey2|Ez2|Hx|Hy|Hz|Dx|Dy|Dz} . THe number 2
            refers as the improved representation of the field. Multiple
            components can be given at the same time, and the calculation will
            be parallelized if the 'parallel' command is employed in the script
            by separating them with a comma (no spaces).
    rimc    must be {r|i|m|c} and indicates if we want to output the real part
            (r), the imaginary part (i), the magnitude (m) or the complex
            value (c) of the fields.

    dz propagation step to be used to output results

    nx number of x points in the output file

    ny number of y points in the output file

    filename the name of the output file which will contain the field
*/
int commands::c_propagation(parsefile *obj, int argc,char *argv[])
{
    structure *p;
    enum rimco_e_t rimc;

    cout<<endl;

    if(argc==7){
        if((p=dynamic_cast<structure *>(obj))==NULL)
            throw parsefile_commandError("Incorrect dynamic cast: programming"
                " error:-(");

        // Check if the W matrix is empty
        if (p->cur->W.isEmpty()) {
            throw parsefile_commandError(
                "propagation: you should first launch solve.");
        }

        if(strcmp(argv[2],"r")==0) {
            rimc=R;
        } else if(strcmp(argv[2],"i")==0) {
            rimc=I;
        } else if(strcmp(argv[2],"m")==0) {
            rimc=M;
        } else if(strcmp(argv[2],"c")==0) {
            rimc=C;
        } else {
            throw parsefile_commandError("propagation: unrecognized rimc."
                " Should be {r|i|m|c}.");
        }

        // At first, we represent the excitation fields in the eigenspace
        int snux=p->dimx/2+1;
        int snuy=p->dimy/2+1;

        cout << "Matrix size: "<<snuy<<" x "<<snux<<"\n";
        cout<<"Calculating eigenmode representation of field excitations.\n";
        cout.flush();

        db_matrix N(p->sec_list[0].W.getNrow(),1);

        if(p->sec_list[p->number_of_sections-1].sWm.isEmpty()) {
            cout << "The backward excitation has not been specified. I will"
                " consider it as absent.\n";
            p->sec_list[p->number_of_sections-1].sWm = N;
        }

        if(p->sec_list[0].sWp.isEmpty()) {
            cout << "The forward excitation has not been specified. I will"
                " consider it as absent.\n";
            p->sec_list[0].sWp = N;
        }

        // p->excitationCalculated=false;
        if(!p->excitationCalculated) {
            cout<<"Calculating excitation vectors for each interface.\n";

            finterface::inpoutp(p->interf_list, p->number_of_sections,
                p->sec_list);
            p->excitationCalculated=true;
        }

        cout.flush();

        double dz=0.0;
        if(sscanf(argv[3], "%lf", &dz)!=1)
            throw parsefile_commandError("propagation: Can not read the size"
                " of the propagation output step.");

        if (dz<0)
            throw parsefile_commandError("propagation: Negative propagation"
                " output step.");

        int sizex=0;
        int sizey=0;

        cout << "Propagation step: dz = "<< dz << " m.\n";

        // Read the number of points to be used for the ouptut file

        if(sscanf(argv[4], "%20d", &sizex)!=1){
            cerr<<"propagation: error while reading parameter nx.\n";
            return ERROR;
        }

        if(sscanf(argv[5], "%20d", &sizey)!=1){
            cerr<<"propagation: error while reading parameter ny.\n";
            return ERROR;
        }

        // Propagate the excitation field
        cout<<"Starting field propagation.\n";
        cout.flush();

        char *comp=argv[1];
        char *fname=argv[6];
        int k=0;
        int l=0;
        int instance=0;
        int used_threads=1;

        struct propagation_data spd[10];
        struct thread_data td[10];

        // Parallelize the execution of propagation on multiple components
        // if necessary.
        do {
            do {
                ++k;
            } while(argv[1][k]!='\0'&&argv[1][k]!=',');
            if(argv[1][k]==',')
                argv[1][k++]='\0';

            if(argv[6][l]=='\0') {
                throw parsefile_commandError("propagation: if you use multiple"
                " components of the field, you should use multiple file"
                " names. How am I supposed to know where to write, otherwise?");
            }

            do {
                ++l;
            } while(argv[6][l]!='\0'&&argv[6][l]!=',');
            if(argv[6][l]==',')
                argv[6][l++]='\0';

            // Data about calculations to be done
            spd[instance].dz= dz;
            spd[instance].sizex= sizex;
            spd[instance].sizey= sizey;
            // Calculates additional output data only in the first instance
            // of the propagation.
            spd[instance].recordAdditionalData = (instance==0);
            spd[instance].rimc= rimc;
            spd[instance].component= comp;
            spd[instance].fname= fname;

            // Generic data about parallelization
            td[instance].p=p;
            td[instance].instance_number=instance;
            td[instance].thread_counter=&used_threads;
            td[instance].payload = &spd[instance];

            parallelize(&td[instance], propagation_structure,
                instance, p->number_of_threads_allowed,
                (argv[1][k]=='\0')?instance+1:10);

            comp=&argv[1][k];
            fname=&argv[6][l];

            instance++;
        } while(argv[1][k]!='\0');
        // Now instance contains the number of launched instances.
        wait_for_threads(td,instance);

        return 0;
    } else {
        cerr<<"propagation: invalid number of parameters.\n";
        return ERROR;
    }
}

/* PROPAGATION low-level routines */

/** This routine writes on a file the result of a propagation calculation in
    a complete structure. It should be called when the excitation vectors
    have already been calculated.
    @param p the structure on which the calculation should be done.
*/
void *commands::propagation_structure(void *threadarg)
{
    struct thread_data *my_data;
    my_data = (struct thread_data *) threadarg;
    structure *p=my_data->p;

    struct propagation_data *pspd;
    pspd = (struct propagation_data*) my_data->payload;

    double dz = pspd->dz;
    int sizex = pspd->sizex;
    int sizey= pspd->sizey;
    enum rimco_e_t rimc = pspd->rimc;
    char *component=pspd->component;
    char *filename = pspd->fname;

    double z0=0.0;
    double z=0;
    double x_t = 0;
    double z_t = 0;
    double alpha_t=0;
    int nz = 0;
    double dx=p->tot_x/sizex, dy=p->tot_y/sizey;
    int snux=p->dimx/2+1;
    int snuy=p->dimy/2+1;

    bool calcD = false; // True if the D vector should be calculated
    bool calcH = false; // True if the H vector should be calculated
    bool calcz = false; // True if the z component should be calculated
    bool shiftToY = false; // True if the y component should be calculated
    bool improve_representation = false; // True if Lalanne-Jurek repr.

    waitSemaphoreIO();

    if(strcmp(component,"Ex")==0) {
        shiftToY = false;
        calcH = false;
        calcz = false;
        cout << "Propagation of the Ex field.\n";
    } else if(strcmp(component, "Ey")==0) {
        shiftToY = true;
        calcH = false;
        calcz = false;
        cout << "Propagation of the Ey field.\n";
    } else if(strcmp(component, "Ez")==0) {
        shiftToY = false;
        calcH = false;
        calcz = true;
        cout << "Propagation of the Ez field.\n";
    } else if(strcmp(component,"Hx")==0) {
        shiftToY = false;
        calcH = true;
        calcz = false;
        cout << "Propagation of the Hx field.\n";
    } else if(strcmp(component, "Hy")==0) {
        shiftToY = true;
        calcH = true;
        calcz = false;
        cout << "Propagation of the Hy field.\n";
    } else if(strcmp(component, "Hz")==0) {
        shiftToY = false;
        calcH = true;
        calcz = true;
        cout << "Propagation of the Hz field.\n";
    } else if(strcmp(component, "Dy")==0) {
        shiftToY = true;
        calcH = false;
        calcD = true;
        calcz = false;
        cout << "Propagation of the Dy field.\n";
    } else if(strcmp(component, "Dx")==0) {
        shiftToY = false;
        calcH = false;
        calcD = true;
        calcz = false;
        cout << "Propagation of the Dx field.\n";
    }else if(strcmp(component, "Dz")==0) {
        shiftToY = false;
        calcH = false;
        calcD = true;
        calcz = true;
        cout << "Propagation of the Dz field.\n";
    }else if(strcmp(component, "Ex2")==0) {
        shiftToY = false;
        calcH = false;
        calcD = false;
        calcz = false;
        improve_representation = true;
        cout << "Propagation of the Ex (improved) field.\n";
    }else if(strcmp(component, "Ey2")==0) {
        shiftToY = true;
        calcH = false;
        calcD = false;
        calcz = false;
        improve_representation = true;
        cout << "Propagation of the Ey (improved) field.\n";
    }else if(strcmp(component, "Ez2")==0) {
        shiftToY = false;
        calcH = false;
        calcD = false;
        calcz = true;
        improve_representation = true;
        cout << "Propagation of the Ez (improved) field.\n";
    }else {
        postSemaphoreIO();
        throw parsefile_commandError("propagation: unrecognized field"
            " type. Should be {Ex|Ey|Ez|Hx|Hy|Hz|Dx|Dy|Dz|Ex2|Ey2|Ez2}.");
    }

    FILE *f=fopen(filename,"w");
    if(f==NULL) {
        postSemaphoreIO();
        throw parsefile_commandError("propagation: Can not open output"
            " file:-(");
    }

    FILE *f1=NULL;
    // Additional data is recorded only if needed.
    if (p->additional_output_data.should_record_integral &&
        pspd->recordAdditionalData) {
        f1=fopen(p->additional_output_data.window_file_name.c_str(),"w");
        if(f1==NULL) {
            fclose(f);
            postSemaphoreIO();
            throw parsefile_commandError("propagation: Can not open output"
                " file for the integral calculation. The error may be on"
                " one of the previous outdata commands.");
        }
    }

    FILE *f2=NULL;
    if (p->additional_output_data.should_record_generation_rate&&
        pspd->recordAdditionalData) {
        f2=fopen(p->additional_output_data.generation_file_name.c_str(), "w");
        if(f2==NULL) {
            fclose(f);
            if(f1!=NULL) {
                fclose(f1);
            }
            postSemaphoreIO();
            throw parsefile_commandError(
                "propagation: Can not open output file"
                " for the generation rate calculation. The error may be on"
                " one of the previous outdata commands.");
        } else if (calcH) {
            fclose(f);
            fclose(f2);
            if(f1!=NULL) {
                fclose(f1);
            }
            postSemaphoreIO();
            throw parsefile_commandError("Additional output: cannot"
                " compute generation rate while propagating H.");
        } else {
            // in case of generation rate calculation, the absorptance is
            // computed
            // in order to check to corectness.
            p->additional_output_data.
                generation_absorptance_power_interpolated = 0.0;
            p->additional_output_data.generation_nz = 0;
            p->additional_output_data.generation_real_from_z0 = 0.0;
            p->additional_output_data.generation_real_to_z1 = 0.0;

            cout << "Additional output: record the generation rate in the"
                " propagation command.";
            cout << " every "<<
                p->additional_output_data.generation_radial_step_size<<" m"
                    << endl;
        }
    }
    postSemaphoreIO();
    // Launch the calculations for each section of the structure
    for(int sec = 0; sec<p->number_of_sections; ++sec) {
        waitSemaphoreIO();
        cout << "Considering section "<<(sec+1)<<
            " of "<<p->number_of_sections<<".\n";

        cout << "[0%        25%          50%        75%         100%]"
            << endl;
        cout << "[";
        postSemaphoreIO();
        // We create an alias for the current section, for simplicity.
        section &c_section=p->sec_list[sec];

        // We check if the excitation is ok.
        if(c_section.sWp.isEmpty()||c_section.sWm.isEmpty())
            throw parsefile_commandError("propagation: excitation missing.");

        // We compute the Toeplitz epsilon and mu matrices that will be
        // used in the next functions
        db_matrix epsilonxy;
        db_matrix epsz;
        db_matrix muz;

        // We compute the constants epsilon and mu if needed
        getEpsilonAndMu(c_section, calcD, calcz,
            calcH,improve_representation, epsilonxy, epsz, muz);

        // We propagate the field in the current section
        process_section(dx, dy, dz, c_section, alpha_t, z0,
            z,f,rimc, x_t, z_t, nz,c_section.sWp,c_section.sWm,
            shiftToY, calcH, snux, snuy, sizex, sizey, f1,
            calcD, epsilonxy, epsz, muz,calcz,improve_representation,
            f2);

    }

    fclose(f);

    if (f1!=NULL) {
        fclose(f1);
        waitSemaphoreIO();
        cout << "Additional output data: window integral calculation"
            " written on "<<
            p->additional_output_data.window_file_name <<endl;
        postSemaphoreIO();
    }

    if (f2!=NULL) {
        fclose(f2);

        /* The power has to be normalized with the surface in the xy plane.
            since we performed a cylindrical integration up to r1
            double r1 = p->additional_output_data.generation_to_r1;
        */

        waitSemaphoreIO();
        cout << "Additional output data: generation rate calculation "
            "written on "<<
        p->additional_output_data.generation_file_name <<endl;

        cout << "File size: "<<
            p->additional_output_data.generation_to_r1 /
            p->additional_output_data.generation_radial_step_size
            <<" x "<< p->additional_output_data.generation_nz <<" points. "
            <<endl
            << "Absorptance power after interpolation is "
            <<  p->additional_output_data.
                generation_absorptance_power_interpolated
            <<endl;

            p->insertVar("ans",  p->additional_output_data.
                generation_absorptance_power_interpolated);
        postSemaphoreIO();
    }
    waitSemaphoreIO();
    cout << "Field propagation file written: "<< filename <<"\n";
    cout << "File size: "<< sizex<<" x " <<sizey<<" x "<< nz <<" points.\n";
    postSemaphoreIO();

    return parallelize_finish(my_data);
}

/** Here we process all what is needed for the propagation in a section.
    @param f1 the file where the field will be written.
    @param f2 the file where additional output data will be written.
*/
void commands::process_section(double dx, double dy, double dz,
    section &c_section, double &alpha_t, double &z0, double &z,
    FILE *f, enum rimco_e_t rimc, double &x_t, double &z_t, int &nz,
    db_matrix &excitation_p, db_matrix &excitation_m,
    bool applyShift, bool calcH, int snux, int snuy, int sj, int si, 
    FILE *f1,
    bool calcD, db_matrix &epsilonxy, db_matrix &epsz,db_matrix &muz,
    bool calcz,bool improve_representation, 
    FILE *f2)
{

    int shift =0;
    double x_c =0;
    double z_c=0;
    double alpha = 0;
    double tp, tm;
    structure *p;

    p = c_section.father;

    db_matrix pipo; // not used

    int advancement=0;

    // If the section is bent, we need to calculate the coordinates  of the
    // center (x_c, y_c)
    if (c_section.isBent) {
        x_c = x_t + c_section.radius * cos(alpha_t);
        z_c = z_t - c_section.radius * sin(alpha_t);
    }
    // We increment z in order that the sampling is regular.
    for(; z<z0+c_section.tot_z; z+=dz) {
        ++nz;

        // Display the advancement of the section processing:
        while  ((int)((z-z0)/(c_section.tot_z)*50) > advancement) {
            waitSemaphoreIO();
            cout << "*";
            cout.flush();
            postSemaphoreIO();
            advancement++;
        }

        // Distance between the current point and the beginning of the current
        // section
        tp=z-z0;
        // Distance between the current point and the end of the current sect.
        tm=-(c_section.tot_z-tp);

        // Propagate all modes at the given point, by taking into account the
        // propagative as well as the counterpropagative components.
        // Note that to compute Ez, Hx and Hy are needed and
        // in order to compute Hz, Ex and Ey are needed.
        bool needH = (calcH && !calcz) || (!calcH && calcz);
        db_matrix fields = structure::getPropagation(c_section.B, tp, tm,
            excitation_p, excitation_m, needH);

        // compute the field:
        db_matrix out  = getField(c_section, fields,
            applyShift, calcH, snux, snuy, sj, si,
            calcD, epsilonxy,
            epsz, muz, calcz,improve_representation, false, pipo);

        // Angle inside the current section
        if (c_section.isBent)
            alpha = alpha_t + (z-z0)/c_section.radius;

        // Write the calculated field on the output file
        outputfield(out, dx, dy, c_section, alpha, alpha_t, z0,
            z,f,rimc, x_t, z_t, x_c,z_c, f1);

        // compute the generation rate if needed
        if (p->additional_output_data.should_record_generation_rate &&
            f2!=NULL)
        {
            // This trick is to make sort that the first and the last points
            // enter in the calculated range in the integration.
            // However, this solution may be questionable.
            double toll=1e-14;
            
            if (p->additional_output_data.generation_from_z0 < z*(1+toll) &&
                z*(1-toll) < p->additional_output_data.generation_to_z1)
            {
                double coeff = 1.0;
                if(p->additional_output_data.generation_from_z0 >
                    (z-dz)*(1+toll))
                {
                    coeff = 0.5;
                    p->additional_output_data.generation_real_from_z0 = z;
                } else if((z+dz)*(1-toll) > 
                    p->additional_output_data.generation_to_z1)
                {
                    coeff = 0.5;
                    p->additional_output_data.generation_real_to_z1 = z;
                }

                // Calculating the number of points in z is interesting to
                // write how many points have been written in the file at the
                // end of the computation.
                p->additional_output_data.generation_nz ++;
                double dr = p->additional_output_data
                    .generation_radial_step_size;
                double r1 = p->additional_output_data.generation_to_r1;

                // compute 2D generation rate
                 db_matrix generation = getGeneration_rate(c_section,tp, tm,
                    excitation_p, excitation_m,snux,snuy, r1,
                    dr, epsilonxy, epsz, muz);

                double radius;
                double generation_previous = generation(0,0).real();
                // output the field:
                waitSemaphoreIO();
                double integral_calc=0.0;
                fprintf(f2, "%le %le %le\n", 0.0, z, generation(0,0).real());
                for(int j=1; j<generation.getNcol(); ++j) {

                    // calculate the current position
                    radius = j*dr;

                    // write in the file
                    fprintf(f2, "%le %le %le\n", radius, z,
                        generation(0,j).real());

                    // compute the integral of the generation rate to get
                    // the absorptance

                    // Coefficient is 0.5*2.0
                    p->additional_output_data.
                        generation_absorptance_power_interpolated +=
                         coeff*M_PI*
                        (generation_previous*(radius-dr)+generation(0,j).real()
                        *radius)*dr*dz*(h_Planck*CELERITY/p->lambda);

                    generation_previous = generation(0,j).real();
                }
                // gnuplot format:
                fprintf(f2, "\n");
                postSemaphoreIO();
            }
        }
    }
    // Increment the reference position to the beginning of the next section.
    z0+=c_section.tot_z;
    if (c_section.isBent) {
        alpha_t +=  c_section.tot_z/c_section.radius;
        x_t = x_c - c_section.radius * cos(alpha_t);
        z_t = z_c + c_section.radius * sin(alpha_t);
    } else {
        x_t = x_t + c_section.tot_z * sin(alpha_t);
        z_t = z_t + c_section.tot_z * cos(alpha_t);
    }

    // Display the advancement of the section processing:
    while  (50 > (int) advancement) {
        waitSemaphoreIO();
        cout << "*";
        cout.flush();
        postSemaphoreIO();
        advancement++;
    }
    waitSemaphoreIO();
    cout << "]" << endl;
    postSemaphoreIO();
}

/**  Compute the constants epsilon and mu if needed

*/
void commands::getEpsilonAndMu(section &c_section, bool calcD, bool calcz,
    bool calcH, bool improve_representation, db_matrix &epsilonxy,
    db_matrix &epsz, db_matrix &muz)
{
    structure *p;
    p = c_section.father;

    db_matrix epsx;
    db_matrix epsy;
    db_matrix deltax;
    db_matrix deltay;

    // in case of symmetry, determine the factor needed to derive correctly the
    // Toeplitz matrices:
    double factor = 1;
    bool xsymmetry=false, ysymmetry=false;
    if (p->symx == symmetric || p->symx == anti_symmetric)
        xsymmetry = true;

    if (p->symy == symmetric || p->symy == anti_symmetric)
        ysymmetry = true;

    if (p->symx == anti_symmetric || p->symy == symmetric)
        factor = -1;

    if (calcH)
        factor *= -1;

    if (calcD || p->additional_output_data.should_record_generation_rate ||
        improve_representation){
        epsx = c_section.P_fft.toeplitz_sym(xsymmetry, ysymmetry,
                            -factor,-factor, false,false);
        epsy = c_section.Q_fft.toeplitz_sym(xsymmetry, ysymmetry,
                            factor,factor, false,false);

        deltax = c_section.Pm1_fft.toeplitz_sym(xsymmetry, ysymmetry,
                            -factor,-factor, false,false).invert();
        deltax *=-1;
        deltax +=epsx;

        deltay = c_section.Qm1_fft.toeplitz_sym(xsymmetry, ysymmetry,
                            factor,factor, false,false).invert();
        deltay *=-1;
        deltay +=epsy;

        // If the normal field is specified, it is used to get Dx and Dy
        if (c_section.crtr == &c_section.NormalFieldSym) {

            db_matrix Q1 = epsx - deltax*
                c_section.NormalFieldSym.Nxx_fft.toeplitz_sym(xsymmetry,
                ysymmetry, -factor,-factor, false,false);

            db_matrix Q2 = -deltax*
            c_section.NormalFieldSym.Nxy_fft.toeplitz_sym(xsymmetry,
                ysymmetry, factor,factor, true,true);

            db_matrix Q3 = -deltay*
                c_section.NormalFieldSym.Nxy_fft.toeplitz_sym(xsymmetry,
                ysymmetry, -factor,-factor, true,true);

            db_matrix Q4 = epsy - deltay *
                c_section.NormalFieldSym.Nyy_fft.toeplitz_sym(xsymmetry,
                ysymmetry, factor,factor, false,false);

            epsilonxy = db_matrix::mergeMatrixQuad(Q1,Q2,Q3,Q4);
        }
        else if (c_section.crtr == &c_section.NonDevSym || c_section.crtr ==
            &c_section.NonDev) {

            db_matrix O(epsx.getNrow(),epsx.getNcol());

            // In the case of the Lalanne developpement strategy, alpha is used
            //  to compute Dx and Dy
            double alpha = 0;
            if (c_section.crtr == &c_section.NonDevSym)
                alpha = c_section.NonDevSym.getAlpha();
            else
                alpha = c_section.NonDev.getAlpha();

            db_matrix Q1 = epsx - deltax*alpha;

            db_matrix Q4 = epsy - deltay * (1-alpha);
            epsilonxy = db_matrix::mergeMatrixQuad(Q1,O,O,Q4);

        } else {
                db_matrix O(epsx.getNrow(),epsx.getNcol());
                epsilonxy = db_matrix::mergeMatrixQuad(epsx,O,O,epsy);
        }
    }
    if (calcz || p->additional_output_data.should_record_generation_rate) {
        epsz= c_section.R_fft.toeplitz_sym(xsymmetry, ysymmetry,
                factor,-factor, false,false).invert();
    }
    if (calcz && calcH) {
        muz= c_section.O_fft.toeplitz_sym(xsymmetry, ysymmetry,
                -factor,factor, false,false).invert();
    }
}

/** Interpolate the 2D generation rate from r=0 to r=r1
    G(r,z) = integral theta=0..2pi (G(r cos theta, r sin theta, z) d theta
    G(x,y,z)= Im(eps)/ (2 hbar) * |E|^2

    This function returns a 1xN matrix (vector) containing in each column the
    calculated value of the generation rate at the current radius, which is
    discretized.
*/

/*  Those issues introduce minor precision losses in some very particular 
    situations: 
    1. When r1 is an odd number, last point in r-integration can be lost
    depending on the rounding.
    2. When tot_x or tot_y is an odd number, ix might become negative
    depending on the rounding. (ix = shift_x + jx).
*/
db_matrix commands::getGeneration_rate(section &c_section,double tp, double tm,
    db_matrix &excitation_p, db_matrix &excitation_m,
    int snux, int snuy,double r1,double dr,
    db_matrix &epsilonxy,
    db_matrix &epsz,db_matrix &muz)
{
    structure *p;
    p = c_section.father;
    double tol = 1.0e-10;    // Tolerance for comparing floats

    if (r1 < 0){
        waitSemaphoreIO();
        cerr << "generation_rate: Integration radius must be greater than "
                "zero\n";
        postSemaphoreIO();
    }

    int nbOfPoint_x,nbOfPoint_y;

    nbOfPoint_x = p->tot_x / dr;
    nbOfPoint_y = p->tot_y / dr;

    // nbOfPoint must be even so that we have a well defined center:
    if (int(nbOfPoint_x/2) != (double)nbOfPoint_x/2)
        nbOfPoint_x++;
    if (int(nbOfPoint_y/2) != (double)nbOfPoint_y/2)
        nbOfPoint_y++;

    // 1.5 is a way to avoid rounding problems.
    // Equivalent to round(r1/dr+1) where the 1 is needed to include the
    // "zero point" which is the center of the disk.
    int nbOfPoint_radius = (r1/dr + 1.5);

    int shift_x = nbOfPoint_x / 2;   // Index of the origin
    int shift_y = nbOfPoint_y / 2;

    // Propagate all modes at the given point, by taking into account the
    // propagative as well as the counterpropagative components.
    db_matrix fields = structure::getPropagation(c_section.B, tp, tm,
        excitation_p, excitation_m, false);

    // compute Ex and Ey
    db_matrix Ey_cartesian;
    db_matrix Ex_cartesian =  getField(c_section,fields,
        false, false, snux, snuy, nbOfPoint_x, nbOfPoint_y,
        false, epsilonxy,epsz,muz,false,true,
        true,Ey_cartesian);

    // Propagate all modes at the given point, by taking into account the
    // propagative as well as the counterpropagative components.
    // Ez is deduced form Hx and Hy which are calculated here:
    fields = structure::getPropagation(c_section.B, tp, tm,
        excitation_p, excitation_m, true);

    // compute Ez
    db_matrix pipo; // not used yet needed
    db_matrix Ez_cartesian =  getField(c_section,fields,
        false, false, snux, snuy, nbOfPoint_x, nbOfPoint_y,
        false, epsilonxy,epsz,muz,true,true,
        false, pipo);

    // compute the interpolated Generation rate
    db_matrix G2D(1,nbOfPoint_radius);  // Result of the integral in the angle,
                                        // for varying radiuses.
    complex<double> eps;                // Value of epsilon in the current point

    double r;       // radius of the point in real space (m)
    double rad;     // radius of the point in the units of dr
    double x, y;    // Coordinates of the point in real space (m)
    int ix, iy;     // Index of the point

    int ind;    // The index of each point in the radius-array. In other words,
                // the index of the ring, to which this point belongs to.
    int count[nbOfPoint_radius];// Counts the number of points in each ring

    for (int i=0; i<nbOfPoint_radius; ++i)
        count[i]=0;

    int min_int_range=nbOfPoint_radius;
    if(nbOfPoint_x/2<min_int_range)
        min_int_range=nbOfPoint_x/2;
    if(nbOfPoint_y/2<min_int_range)
        min_int_range=nbOfPoint_y/2;

    // We loop through the points and distribute them to different rings
    for (int jx = -min_int_range+1; jx < min_int_range; ++jx){
        for (int jy = -min_int_range+1; jy < min_int_range; ++jy){

            rad = sqrt(pow((double)jx,2) + pow((double)jy,2));
            r = rad*dr;
            // Check if we are inside the integration region
            if (r < r1*(1.0+tol)){
                // Hand-made floor function. Obtain the ring index.
                ind = (int)(rad+0.5);
                // Boundary between the ring #i and #i+1 lies at r=dr*(i+0.5)
                // By construction of rad, ind never has a rounding problem

                ix = shift_x + jx;    // Index and real position of each point
                iy = shift_y + jy;
                x = (double)ix *p->tot_x/nbOfPoint_x - p->tot_x/2.0;
                y = (double)iy *p->tot_y/nbOfPoint_y - p->tot_y/2.0;

                // The refractive index in the point
                eps = EPS_0* c_section.expectedRefractiveIndex(x, y)*
                    c_section.expectedRefractiveIndex(x, y);
                // Contribution of this point to the ring-average
                // (\theta - integral)
                G2D(0,ind) -= 2*M_PI*(pow(abs(Ex_cartesian(iy,ix)),2)
                    + pow(abs(Ey_cartesian(iy,ix)),2)
                    + pow(abs(Ez_cartesian(iy,ix)),2))
                    * eps.imag() /(2.0*h_Planck);
                ++count[ind];
            }
        }
    }
    for (int ir = 0; ir < nbOfPoint_radius; ++ir){
        if (count[ir]==0)
            G2D(0,ir)=0;
        else
            G2D(0,ir) = G2D(0,ir)/(double)count[ir];
    }

    return G2D;
}

/** Compute the field in the traditional space for E, H or D, depending on
    boolean variables.
    If an improved representation is needed, it is possible to compute at the
    same time Ex and Ey in order to save time (using additional Ey and Ey
    variables)
*/
db_matrix commands::getField(section &c_section,db_matrix &fields,
    bool applyShift, bool calcH, int snux, int snuy, int sj, int si,
    bool calcD, db_matrix &epsilonxy,
    db_matrix &epsz,db_matrix &muz,bool calcz,bool improve_representation,
    bool additionnalEy, db_matrix &outEy)
{
    structure *p;

    p = c_section.father;

    db_matrix out;

    if (!improve_representation) {
        // compute the field based on the boolean variable calcH, calcz, calcD
        db_matrix fields_fourier = getFourierField(c_section,fields, calcH,
            calcD, epsilonxy, epsz, muz, calcz);

        // The mode matrix will be constructed according to standard
        // rules followed by FFT.
        db_matrix mode = fields_fourier.vector2fft(p->symx,
             p->symy,snux,snuy,applyShift,calcH,calcz);

        // Zero pad the calculated mode, in order to increase the
        // number of points with which it will be represented.
        out=mode.zero_pad(si,sj).fft2();
    } else {
        // not possible (yet!)
        if (calcD || calcH) {
            throw parsefile_commandError("Calculations on the H and D fields "
                "are not (yet!) implemented for the improved representation.");
            return out;
        }

        if (calcz){
            // In order to improve the representation of field, Dz is computed
            // since it is contiunous. Then Dz_cartesian is divided by epsz:

             // compute Dz
            db_matrix D = getFourierField(c_section, fields, false,
                true, epsilonxy, epsz, muz,true);

            // The mode matrix will be constructed according to standard
            // rules followed by FFT.
            db_matrix mode = D.vector2fft(p->symx, p->symy,snux,snuy,false,
                false, true);
            // Zero pad the calculated mode, in order to increase the
            // number of points with which it will be represented.
            db_matrix Dz_cartesian=mode.zero_pad(si,sj).fft2();

            out = Dz_cartesian;
            double x,y;
            complex<double> epsilon;

            // we compute Ez
            for (int i = 0; i< Dz_cartesian.getNrow(); ++i){
                for (int j = 0; j< Dz_cartesian.getNcol(); ++j){

                    x = ((double)j/(double)Dz_cartesian.getNcol() -0.5) *
                         p->tot_x;
                    y = ((double)i/(double)Dz_cartesian.getNrow() -0.5) *
                        p->tot_y;
                    epsilon = EPS_0* c_section.expectedRefractiveIndex(x, y)*
                        c_section.expectedRefractiveIndex(x, y);

                    out(i,j) = Dz_cartesian(i,j)/ epsilon;
                }
            }
        } else {
            // In order to improve the representation of field, Normal field is
            // used to compute the electric field based on the continuity of
            // the  field: D normal is continuous and E perpendicular is
            // continuous:

             // compute Dx et Dy
            db_matrix D = getFourierField(c_section, fields, false,
                true, epsilonxy, epsz, muz,false);

            // compute Dx
            // The mode matrix will be constructed according to standard
            // rules followed by FFT.
            db_matrix mode = D.vector2fft(p->symx, p->symy,snux,snuy,false,
                false,false);
            // Zero pad the calculated mode, in order to increase the
            // number of points with which it will be represented.
            db_matrix Dx_cartesian=mode.zero_pad(si,sj).fft2();

            // compute Dy
            // The mode matrix will be constructed according to standard
            // rules followed by FFT.
            mode = D.vector2fft(p->symx,p->symy,snux,snuy,true,false, false);
            // Zero pad the calculated mode, in order to increase the
            // number of points with which it will be represented.
            db_matrix Dy_cartesian=mode.zero_pad(si,sj).fft2();

            // free the space
            D.kill();

             // compute Ex et Ey
            db_matrix E = getFourierField(c_section, fields, false,
                false, epsilonxy, epsz, muz,false);

            // compute Ex
            // The mode matrix will be constructed according to standard
            // rules followed by FFT.
            mode = E.vector2fft(p->symx, p->symy,snux,snuy,false,false,false);
            // Zero pad the calculated mode, in order to increase the
            // number of points with which it will be represented.
            db_matrix Ex_cartesian=mode.zero_pad(si,sj).fft2();

            // compute Ey
            // The mode matrix will be constructed according to standard
            // rules followed by FFT.
            mode = E.vector2fft(p->symx,p->symy,snux,snuy,true,false,false);
            // Zero pad the calculated mode, in order to increase the
            // number of points with which it will be represented.
            db_matrix Ey_cartesian=mode.zero_pad(si,sj).fft2();

            // free the space
            E.kill();
            mode.kill();

            double x,y;
            complex<double> epsilon;

            // If the normal field is specified the formula used is:

            // [ Ex_cartesian_improved ] _ [ Nxx_cartesian  Nxy_cartesian ] *
            // [Dx_cartesian/epsx] +
            // [ Ey_cartesian_improved ] - [ Nxy_cartesian  Nyy_cartesian ] *
            // [Dy_cartesian/epsy]

            // [ 1- Nxx_cartesian       - Nxy_cartesian   ] [Ex_cartesian]
            // [ - Nxy_cartesian        1 - Nyy_cartesian ] [Ey_cartesian]
            if (c_section.crtr == &c_section.NormalFieldSym) {
                if (applyShift || additionnalEy) {
                    outEy = Ey_cartesian;
                    // we compute Ey

                    for (int i = 0; i< Ey_cartesian.getNrow(); ++i){
                        for (int j = 0; j< Ey_cartesian.getNcol(); ++j){

                        x = ((double)j /(double)Ey_cartesian.getNcol() -0.5) *
                            p->tot_x;
                        y = ((double)i /(double)Ey_cartesian.getNrow() -0.5) *
                            p->tot_y;

                        epsilon = EPS_0*c_section.expectedRefractiveIndex(x,y)*
                            c_section.expectedRefractiveIndex(x, y);

                        outEy(i,j) = c_section.NormalFieldSym.
                            Nxy_input.expectedNormalField(x,y,p->tot_x,
                            p->tot_y)*
                            Dx_cartesian(i,j) / epsilon +
                            c_section.NormalFieldSym.
                            Nyy_input.expectedNormalField(x,y,p->tot_x,
                            p->tot_y)*
                            Dy_cartesian(i,j) / epsilon -
                            c_section.NormalFieldSym.
                            Nxy_input.expectedNormalField(x,y,p->tot_x,
                            p->tot_y) *
                            Ex_cartesian(i,j)  +
                            (complex<double>(1,0) -c_section.NormalFieldSym.
                            Nyy_input.expectedNormalField(x,y,p->tot_x,
                            p->tot_y))*
                            Ey_cartesian(i,j);
                        }
                    }
                    // We only want Ey
                    if (applyShift)
                        return outEy;
                }
                if (! applyShift)   {
                    out = Ex_cartesian;
                    // we compute Ex
                    for (int i = 0; i< Ex_cartesian.getNrow(); ++i){
                        for (int j = 0; j< Ex_cartesian.getNcol(); ++j){

                        x = ((double)j /(double)Ex_cartesian.getNcol() -0.5) *
                            p->tot_x;
                        y = ((double)i /(double)Ex_cartesian.getNrow() -0.5) *
                            p->tot_y;

                        epsilon = EPS_0*c_section.expectedRefractiveIndex(x,y)*
                            c_section.expectedRefractiveIndex(x, y);

                        out(i,j) = c_section.NormalFieldSym.Nxx_input.
                            expectedNormalField(x,y,p->tot_x,p->tot_y)*
                            Dx_cartesian(i,j) / epsilon +
                            c_section.NormalFieldSym.Nxy_input.
                            expectedNormalField(x,y,p->tot_x,p->tot_y)*
                            Dy_cartesian(i,j) / epsilon +
                            (complex<double>(1,0) - c_section.NormalFieldSym.
                            Nxx_input.expectedNormalField(x,y,p->tot_x,
                            p->tot_y))*
                            Ex_cartesian(i,j)  -
                            c_section.NormalFieldSym.Nxy_input.
                            expectedNormalField(x,y,p->tot_x,p->tot_y)*
                            Ey_cartesian(i,j);
                        }
                    }
                }
            } else if (c_section.crtr == &c_section.NonDevSym
                || c_section.crtr == &c_section.NonDev) {

                // In the case of Lalanne developpement, the same equation
                // can be applied with:
                // Nxx = alpha
                // Nyy = 1 - alpha
                // Nxy = 0
                double alpha = 0;
                if (c_section.crtr == &c_section.NonDevSym)
                    alpha = c_section.NonDevSym.getAlpha();
                else
                    alpha = c_section.NonDev.getAlpha();

                if (applyShift || additionnalEy) {
                    outEy = Ey_cartesian;
                    // we compute Ey

                    for (int i = 0; i< Ey_cartesian.getNrow(); ++i){
                        for (int j = 0; j< Ey_cartesian.getNcol(); ++j){

                        x = ((double)j /(double)Ey_cartesian.getNcol() -0.5) *
                            p->tot_x;
                        y = ((double)i /(double)Ey_cartesian.getNrow() -0.5) *
                            p->tot_y;

                        epsilon = EPS_0*c_section.expectedRefractiveIndex(x,y)*
                            c_section.expectedRefractiveIndex(x, y);

                        outEy(i,j) =    (1-alpha)* Dy_cartesian(i,j) / epsilon-
                                    alpha * Ey_cartesian(i,j);
                        }
                    }
                    // We only want Ey
                    if (applyShift)
                        return outEy;

                }
                if (!applyShift) {
                    out = Ex_cartesian;
                    // we compute Ex
                    for (int i = 0; i< Ex_cartesian.getNrow(); ++i){
                        for (int j = 0; j< Ex_cartesian.getNcol(); ++j){

                            x = ((double)j/(double)Ex_cartesian.getNcol()-0.5)*
                                p->tot_x;
                            y = ((double)i/(double)Ex_cartesian.getNrow()-0.5)*
                                p->tot_y;

                            epsilon = EPS_0* c_section.
                                expectedRefractiveIndex(x, y)*
                                c_section.expectedRefractiveIndex(x, y);

                            out(i,j) = alpha * Dx_cartesian(i,j) / epsilon +
                                (1 - alpha)* Ex_cartesian(i,j);
                        }
                    }
                }
            } else {
                    if (additionnalEy)
                        outEy = Ey_cartesian;
                    if (applyShift)
                        return Ey_cartesian;
                    else
                        return Ex_cartesian;
            }
        }
    }

    return out;
}
/** Compute the Fourier coefficients of fields for E, H or D.
    If we need the magnetic field, the electric field is now transformed
    in the corresponding magnetic field. There is a difference that has been
    taken into account in the getPropagation routine: the sum between the 
    propagative and the counter-propagative component should be done with a
    + sign if we need the E field, and with the - sign if we need the H 
    field.
    Here we are still using an eigenvalue representation, so this task is
    actually quite easy.

    if we want to derive Ez, we first store Hx and Hy since the former is
    obtained from them. The reverse is true for Hz.
*/
db_matrix commands::getFourierField(section &c_section,db_matrix &fields,
    bool calcH,
    bool calcD, db_matrix &epsilonxy,
    db_matrix &epsz,db_matrix &muz,bool calcz)
{
    structure *p;
    p = c_section.father;

    
    db_matrix epsilonz;
    db_matrix mu;
    db_matrix fields_fourier;

    // check wether the symmetry or anti-symmetry is active.
    bool xsymmetry=false, ysymmetry=false;
    if (p->symx == symmetric || p->symx == anti_symmetric)
        xsymmetry = true;

    if (p->symy == symmetric || p->symy == anti_symmetric)
        ysymmetry = true;

    // we obtain now the representation of the propagated field in
    // the Fourier space by multiplying it times the eigenvalues matrix
    if((calcH && !calcz) || (!calcH && calcz)) {
        fields_fourier = c_section.V*fields;
    } else {
        fields_fourier = c_section.W*fields;
        if (calcD && !calcz) {
            fields_fourier = epsilonxy*fields_fourier;
        }
    }

    // In order to compute Ez or Hz, we need to use the following equations
    // [Ez] = 1/w [1/epsz] (Kx  [Hy]  - Ky [Hx])
    // [Hz] = -1/w [1/muz] (Kx  [Ey]  - Ky [Ex])
    // Note: an error was present here and jw was present in the previous
    // equations instead of w. The j is not present since the Kx and Ky
    // matrices do not contain it, it must be added to represent the derivative
    // and thus is simplified with the jw term.

    // See Jérôme Michallon's PhD thesis at page 92.

    if (calcz) {
        double nux=2.0*M_PI/p->tot_x;
        double nuy=2.0*M_PI/p->tot_y;

        // compute the derivative matrices that should be included as input
        db_matrix Kr = c_section.R_fft.toeplitz_deriv_sym(xsymmetry,ysymmetry,
                nux,nuy,1,0);
        db_matrix Kz = c_section.R_fft.toeplitz_deriv_sym(xsymmetry,ysymmetry,
                nux,nuy,0,1);

        db_matrix O2(Kr.getNrow(),Kr.getNcol());

        db_matrix derivative=db_matrix::mergeMatrixQuad(-Kz,Kr,O2,O2);

        fields_fourier = derivative*fields_fourier;

        // compute epsilon z
        if (!calcD) {
            if (calcH){
                mu = db_matrix::mergeMatrixQuad(muz,O2,O2,O2);
                fields_fourier = -mu*fields_fourier;
            } else {
                epsilonz = db_matrix::mergeMatrixQuad(epsz,O2,O2,O2);
                fields_fourier = epsilonz*fields_fourier;
            }
        }
        // multiply by 1/w:
        fields_fourier /= p->getOmega();
        
        /* Version with an error
        // multiply by 1/jw:
        fields_fourier /= complex<double>(0,1)*p->getOmega();*/
    }
    return fields_fourier;
}

/** Write on the given file the field calculated at the given quota.
    @param f1 the file where the field will be written.
*/
void commands::outputfield(db_matrix &out, double dx, double dy,
    section &c_section, double &alpha, double &alpha_t, double &z0, double &z,
    FILE *f, enum rimco_e_t rimc, double &x_t, double &z_t, double &x_c,
    double &z_c, FILE *f1)
{
    double x;
    double y;
    double x_p = 0;
    double z_p = 0;

    double r_a = 0;
    double z_a = 0;
    double norm = 0;

    double ksinthetax=c_section.father->ksinthetax;
    double ksinthetay=c_section.father->ksinthetay;

    for(int i=0; i<out.getNrow(); ++i) {
        for(int j=0; j<out.getNcol(); ++j) {
            // calculate the current position (we are in curvilinear
            // coordinates).
            x = j*dx-(out.getNcol()-1)/2.0*dx;
            y = i*dy-(out.getNrow()-1)/2.0*dy;
            out(i,j)*=exp(complex<double>(0,1)*ksinthetax*x);
            out(i,j)*=exp(complex<double>(0,1)*ksinthetay*y);
            //cout<<"; mv="<<exp(complex<double>(0,1)*ksinthetax*x);

            // If this section is bent or it has inherited an angle from
            // previous bendings, take it into account.

            if (c_section.isBent) {
                r_a = c_section.radius-x;
                x_p = x_c - r_a * cos(alpha);
                z_p = z_c + r_a * sin(alpha);
            } else {
                z_a = z-z0;
                x_p = x_t + z_a * sin(alpha_t) + x * cos(alpha_t);
                z_p = z_t + z_a * cos(alpha_t) - x * sin(alpha_t);
            }

            // If additional data have to be written, do it.
            if (c_section.father->additional_output_data.should_record_integral
                && f1!=NULL &
                abs(x-c_section.father->additional_output_data.window_centerx)
                <
                c_section.father->additional_output_data.window_width/2.0 &&
                abs(y-c_section.father->additional_output_data.window_centerx)
                <
                c_section.father->additional_output_data.window_height/2.0)
            {
                norm += abs(out(i,j)*out(i,j))*dx*dy;
            }

            // Output of the fields.
            waitSemaphoreIO();
            if(rimc==C) {
                fprintf(f, "%le %le %le %le %le\n",
                    x_p, y, z_p, out(i,j).real(),out(i,j).imag());
            } else if(rimc==R) {
                fprintf(f, "%le %le %le %le\n",
                    x_p, y, z_p, out(i,j).real());
            } else if(rimc==I) {
                fprintf(f, "%le %le %le %le\n",
                    x_p, y, z_p, out(i,j).imag());
            } else if(rimc==M) {
                fprintf(f, "%le %le %le %le\n",
                    x_p, y, z_p, abs(out(i,j)));
            }
            postSemaphoreIO();
        }
    }

    // If the field norm should be calculated, write it on the second output
    // file
    if (c_section.father->additional_output_data.should_record_integral &&
        f1!=NULL) {
        norm = sqrt(norm);
        waitSemaphoreIO();
        fprintf(f1, "%le\n", norm);
        postSemaphoreIO();
    }
}
