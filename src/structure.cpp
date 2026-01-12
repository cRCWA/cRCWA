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

    Davide Bucci, 2008-2026
    Jérôme Michallon, 2012-2014
*/

#include <fstream>

#include "block_matrix.h"
#include "structure.h"
#include "section.h"
#include "finterface.h"
#include "phys_constants.h"
#include "parallelize.h"


using namespace std;

structure::structure(numParser *p) : parsefile(p)
{
    sec_list =NULL;
    interf_list =NULL;
    reset();
}

/** Delete the previous structure.

*/
void structure::reset(void)
{
    if(sec_list!=NULL)
        delete[] sec_list;
    if(interf_list!=NULL)
        delete[] interf_list;

    sec_list = new section[INIT_SIZE];
    interf_list = new finterface[INIT_SIZE];
    actual_size = INIT_SIZE;

    for(int i=0; i<actual_size; ++i)
        sec_list[i].setFather(this);

    cur = sec_list;
    number_of_sections = 1;
    excitationCalculated = false;

    HxField = false;
    HyField = false;
    calcPropagation = false;

    ensureConvergence = false;
    number_of_threads_allowed=1;

    dimx = 0;
    dimy = 0;
    lambda = 0;
    tot_x = 0;
    tot_y = 0;
    ksinthetax=0;
    ksinthetay=0;

    // set the default value for symmetry
    symx = previous_version;
    symy = previous_version;

/*    // Set the default matrix creation strategy for PMLs
    crtr = &PMLafterOPT;    */

    additional_output_data.should_record_integral = false;
    additional_output_data.should_record_generation_rate = false;

    bloch_eigvect.kill();
    bloch_eigval.kill();
}


structure::~structure()
{
    if(sec_list!=NULL)
        delete[] sec_list;
    if(interf_list!=NULL)
        delete[] interf_list;
}

/** Sum up all contributions of the propagation modes
*/
db_matrix structure::getPropagation(db_matrix &B, double tp, double tm,
    db_matrix &excitation_p, db_matrix &excitation_m, bool calcH)
{
    db_matrix fields(B.getNrow(),1);

    for (int jj=0; jj<B.getNrow();++jj) {
        // At first, we obtain the propagation constants associated
        // to each eigenmode.

        complex<double>e=B(jj,jj);

        // We then propagate the excitation field
        // For the propagative component, from left to right:
        fields(jj,0)=excitation_p(jj,0)*exp(-complex<double>(0,1)*e*tp);

        // For the counterpropagative component, from right to left:
        if(calcH) {
            fields(jj,0)-=excitation_m(jj,0)*exp(complex<double>(0,1)*e*tm);
        } else {
            fields(jj,0)+=excitation_m(jj,0)*exp(complex<double>(0,1)*e*tm);
        }
    }
    return fields;
}

/** Callback method for running the parallel calculations.
*/
void *structure::calcV(void *threadarg)
{
    struct thread_data *my_data;

    my_data = (struct thread_data *) threadarg;
    int i=my_data->instance_number;
    structure *p=my_data->p;

    section *q= &(p->sec_list[i]);
    waitSemaphoreIO();
    printf("Fabricating V matrix, section: %d, threads in parallel: %d \n", i,
        (*my_data->thread_counter));
    postSemaphoreIO();
    q->V=q->assembleVmatrix();

    return parallelize_finish(my_data);
}

/** Calculate the power flowing in a given area transverse to the propagation
    direction (command "monitor").
    
*/
double structure::do_monitor(const double z, const double wx, const double wy,
    const double px, const double py)
{
    double tp,tm;
    double power;
    // Find the section in which z is included:
    double ztot = 0;    // coordinate at the end of the section considered
    int i = 0;
    db_matrix *e;

    section *q= (sec_list+i);
    while (ztot <= z && i < number_of_sections) {
        q= (sec_list+i);
        ztot += q->tot_z;
        ++i;
    }
    cout << "section considered: " << i << endl;
    
    
    cout << "Considering the surface of size "<<wx<<" m x "<<wy<<
        " m at location "<<px<<" m x "<<py<<" m\n";


    // Distance between the current point and the end of the current sect.
    tm= ztot - z;
    // Distance between the current point and the beginning of the section
    tp= (q->tot_z - tm);

    tm *=-1;

    // Propagate all modes at the given point, by taking into account the
    // propagative as well as the counterpropagative components.
    db_matrix excitationE = structure::getPropagation(q->B, tp, tm,
        q->sWp, q->sWm, false);

    db_matrix excitationH = structure::getPropagation(q->B, tp, tm,
        q->sWp, q->sWm, true);

    // Sum the contribution of all the existing modes.
    db_matrix E=q->W*excitationE;
    db_matrix H=q->V*excitationH;

    power= section::integralPoynting_rectangle(this,E, H,0,wx,wy,px,py);
    cout << "power at z = " << z << " (in W): "<<power<<endl;
    return power;

}

/** Calculate the net power flowing in a given section
    @param z position in the structure.
*/
double structure::do_powerz(const double z)
{
    section *q;
    db_matrix *e;
    double tp,tm;
    double ztot=0;

    // Find the section in which z is included:
    int i = 0;
    q= (sec_list+i);
    while (ztot <= z && i < number_of_sections) {
        q= (sec_list+i);
        ztot += q->tot_z;
        ++i;
    }
    cout << "section considered: " << i << endl;

    // Distance between the current point and the end of the current sect.
    tm= ztot - z;
    // Distance between the current point and the beginning of the section
    tp= (q->tot_z - tm);

    tm *=-1;

    // Propagate all modes at the given point, by taking into account the
    // propagative as well as the counterpropagative components.
    db_matrix excitationE = getPropagation(q->B, tp, tm, q->sWp, q->sWm, false);

    db_matrix excitationH = getPropagation(q->B, tp, tm, q->sWp, q->sWm, true);

    // Sum the contribution of all the existing modes.
    db_matrix E=q->W*excitationE;
    db_matrix H=q->V*excitationH;

    return section::integralPoynting(this,E, H,0)*tot_x*tot_y;

}

/** Implement the "section" command, by setting the length of the current
    section if it is the first one, or by adding a new section to the current
    structure.
    @param tl the length of the section.
*/
void structure::do_section(const double tl)
{
    // Just a check. It may happen the user asks for a propagation
    // in a section having a negative length. This will result in a
    // almost certain convergence problem.
    if(tl<0) {
        throw parsefile_commandError("section: the length should not"
            " be negative.\n");
    }
    // Section is a optional command. We need to determine if we need
    // to add a new section to the array contained in the global
    // structure. We check the total lenght of the current section.
    // If it is different from zero, we need to create a new one.

    if (cur->tot_z==0) {
        // Change the length of the current structure.

        cur->tot_z = tl;
        cout<<"Length of the section 1: "<<tl<<" m. \n";

    } else {

        // We first need to check if we have enough space in the
        // vectors containing all the structure and interface
        // descriptions in order to store the new interface.
        // If not, we should allocate some bigger vectors and
        // then continue.
        if(number_of_sections + 1 >= actual_size) {
            ensureSize(actual_size + INIT_SIZE);
        }

        // Create a new section and operate from it.
        cur=&(sec_list[number_of_sections++]);
        cur->tot_z = tl;
        cout<<"Length of the section "<<number_of_sections<<": "<<tl<<" m.\n";
    }

    // initialise the default matrix development
    if (symx == symmetric || symx == anti_symmetric ||
        symx == no_symmetry ||  symy == symmetric ||
        symy == anti_symmetric || symy == no_symmetry) {

            cur->NonDevSym.setAlpha(-1.0);
            cur->crtr = &cur->NonDevSym;
    } else {
        cur->crtr = &cur->PMLafterOPT;
    }

}


/** Implement the "select" command.
    @param number the number of the section to be selected (starts from 1).
*/
void structure::do_select(unsigned int number)
{

    // At first, we check that everything is ok.

    if (number < 1 || number > number_of_sections)
        throw parsefile_commandError("select: invalid section number.\n");

    // The internal numbering starts from 0, as with C/C++

    cur=&(sec_list[--number]);
    
    cout<<"Selected section "<<(number+1)<<" of "<<number_of_sections
        <<endl;

}

/** Implement the "wants" command, by setting the matrices that should be
    conserved during the calculations.
    @param code the code {Hx|Hy|both|propagation}.
*/
void structure::do_wants(const char *code)
{
    // p contains the current object and should therefore be a structure


    // TODO: DB: I patched the command by setting calcPropagation even in
    // cases where it should not be strictly needed.
    if(strcmp(code,"Hx")==0) {
        HxField = true;
        calcPropagation=true;        // TODO: this is not optimum
        cout<< "Hx field will be available for calculation after solve.\n";
    } else if(strcmp(code,"Hy")==0) {
        HyField = true;
        calcPropagation=true;        // TODO: this is not optimum
        cout<< "Hy field will be available for calculation after solve.\n";
    } else if(strcmp(code,"both")==0) {
        HxField = true;
        HyField = true;
        calcPropagation=true;        // TODO: this is not optimum
        cout << "Hx and Hy fields will be available for calculation after"
            " solve.\n";

    } else if(strcmp(code,"propagation")==0) {
        calcPropagation=true;
        cout << "Matrices will be retained to compute propagation.\n";
    } else
        throw parsefile_commandError("wants: unrecognized field option."
            " Should be {Hx|Hy|both|propagation}.");

}


/** Run the "assemble" operation on the structure, i.e. calculate the S-matrices
    associated to the excitations on each interface as well as the global
    S-matrix.
*/
void structure::do_assemble(void)
{
    // Check if the C matrix is empty
    if (cur->W.isEmpty()) {
        throw parsefile_commandError(
            "assemble: the problem must be set and modes calculated with "
            "solve before assembling matrices for propagation.");
    }
    if (!calcPropagation)
        throw parsefile_commandError("assemble: you should use the "
                            "wants command before assemble to specify you "
                            "want to calculate the field propagation.");
    cout<<"Assembling sections\n";

    // We need to take into account every interface.
    // If we have N sections, there will be N-1 interfaces, since
    // the input and output interface are not counted.


    // We need some temporary matrices. At each iteration, we consider
    // section t.
    db_matrix Wt_m1;   // Inverse of W matrix of the section t
    db_matrix Vt_m1;   // Inverse of V matrix of the section t
    db_matrix Wtp1_m1; // Inverse of W matrix of the section t + 1
    db_matrix Vtp1_m1; // Inverse of V matrix of the section t + 1
    db_matrix Pt;      // Propagation matrix section t
    db_matrix Ptp1;    // Propagation matrix section t + 1

    vector<struct thread_data> td(number_of_sections);
    int used_threads=1;

    // Calculate the V matrices for each section: PARALLELIZED
    for(int i=0; i<number_of_sections; ++i) {
        td.at(i).p=this;
        td.at(i).instance_number=i;
        td.at(i).thread_counter=&used_threads;
        //  p->number_of_threads_allowed
        parallelize(&td.at(i), calcV, i, number_of_threads_allowed,
            number_of_sections);
    }
    wait_for_threads(td, number_of_sections);

    // Calculate the S matrices for each interface

    // In the code, we use variable i instead of t as an index and t=i-1
    for(int i=1; i<number_of_sections; ++i) {
        cout<<"Creating S matrix for interface "<<i<<" of "<<
            (number_of_sections-1) << ".\n";
        cout.flush();

        // We have everything for calculating the interface between
        // the section i-1 and the section i.
        // i.e.    ->>  t                  t+1

        section &t=sec_list[i-1]; // Reference to the current section
        section &tp1=sec_list[i]; // Reference to the next section

        if(i==1) {
            // We must calculate the matrices for the first section.
            Wt_m1=t.W;
            Wt_m1.invert();
            Vt_m1=t.V;
            Vt_m1.invert();
            Pt=t.propagationMatrix();
        } else {
            // We reuse the matrices we have previously calculated.
            Wt_m1=Wtp1_m1;
            Vt_m1=Vtp1_m1;
            Pt=Ptp1;
        }

        db_matrix &Wtp1 = tp1.W;
        db_matrix &Vtp1 = tp1.V;

        Wtp1_m1=tp1.W;
        Wtp1_m1.invert();
        db_matrix &Wt = t.W;

        Vtp1_m1=tp1.V;
        Vtp1_m1.invert();
        db_matrix &Vt = t.V;

        Ptp1=tp1.propagationMatrix();
        interf_list[i-1].createInterface(Wt_m1, Wtp1, Vt_m1, Vtp1,
            Wtp1_m1, Wt, Vtp1_m1, Vt, Pt, Ptp1);
    }

    // We can now compute the S matrix of the complete calculation region

    cout<<"Creating global S matrix.\n";
    globalS = finterface::createSmatrix(interf_list, number_of_sections);

}

/** Run the "solve" operation on the structure, i.e. calculate all "interesting"
    modes supported by all sections and return those of the last section.
*/
vector<double> structure::do_solve(void)
{
    section *q;
    vector<double> results;
    vector<struct thread_data> td(number_of_sections);
    int used_threads=1;

    // Parallelize the calculation for the allowed threads.
    for(int i=0; i<number_of_sections; ++i) {
        td.at(i).p=this;
        td.at(i).instance_number=i;
        td.at(i).thread_counter=&used_threads;
        q=&sec_list[i];
        if((getKsinthetax()!=0 || getKsinthetay()!=0) &&
            (q->crtr->setAngles(getKsinthetax(),getKsinthetay())))
        {
            throw parsefile_commandError("solve: the current matdev"
                " development is"
                " not compatible with nonzero excitation wave angles.\n"
                " If you have used the matdev command, choose a development"
                " that is compatible. If you have not used it, this means"
                " that you are"
                " working with the default settings and you should change"
                " them.");
        }
        parallelize(&td.at(i), calc_eigenvalues, i, number_of_threads_allowed,
            number_of_sections);
    }
    wait_for_threads(td, number_of_sections);

    // Show what it has been obtained for all sections.
    // All modes have been already calculated, this section does not require
    // parallelism since it is very fast.
    
    for(int i=0; i<number_of_sections; ++i) {
        cout <<"Section "<<(i+1)<<" of "<< number_of_sections<< "."<<endl;
        results= sec_list[i].find_interesting();
    }
    return results;
}

/** Execute the bloch command.
*/
vector<double> structure::do_bloch(void)
{
    // We can not calculate Bloch modes for a structure with just one
    // section.
    if (number_of_sections<2) {
        throw parsefile_commandError(
            "bloch: I can not calculate Bloch modes for a structure with"
            " just one section section ");
    }

    // Check if the S matrix is empty
    // This usually happens when the user forgets to call the assemble
    // command before using bloch.
    if (globalS.Rmp.isEmpty()) {
        throw parsefile_commandError(
            "bloch: you should calculate the global S matrix of the whole"
            " structure via the assemble command before using bloch.");
    }

    return calculateBlochModes();
}

/** Used by the multithreaded system. Calculate all the propagation modes
    (eigenvalue/eigenvector expansion) of one section of the structure.
    The idea is that this function should be launched in parallel for
    several sections: thus, no side effects should be present here, apart
    for the section being calculated.
*/
void *structure::calc_eigenvalues(void *threadarg)
{
    struct thread_data *my_data;

    my_data = (struct thread_data *) threadarg;
    int i=my_data->instance_number;
    structure *p=my_data->p;

    section *q= &(p->sec_list[i]);
    // q contains the current section
    waitSemaphoreIO();
    printf("Now considering section %d of %d, threads in parallel: %d \n", i+1,
        p->number_of_sections,
        (*my_data->thread_counter));
    postSemaphoreIO();
    q->modes_calculation();

    return parallelize_finish(my_data);
}

void structure::set_ensureConvergence(bool s)
{
    ensureConvergence=s;
    if(ensureConvergence)
        cout<<"Dust will be swept under the carpet.\n"
            <<"(All eigenvalues having a positive imaginary part will be "
            <<"forced to have it negative).\n";
}

/**
    Set the size of the calculation window
    Parameters:
    sx: total size of the calculation window
    sy: total size of the calculation window
*/
void structure::do_size(double sx, double sy)
{
    if (sx<=0 || sy <=0) {
        throw parsefile_commandError("size: you should specify values"
            " greater than zero.");
    }

    // p contains the current object and should therefore be a structure
    tot_x = sx;
    tot_y = sy;
    cout << "Calculation window size set to: "<<sx<<" m x "<<sy<<" m.\n";
}

/**
    harmonics nx ny

    Parameters:
    nx: number of harmonics on the x axis
    ny: number of harmonics on the y axis
*/
void structure::do_harmonics(int nx, int ny)
{
    if (nx<1 || ny <1) {
        throw parsefile_commandError("harmonics: you should specify values"
            " greater than zero.");
    }

    /* Size of the FFT matrices. This is more useful in calculations
        and take into account the total number of positive and negative
        frequency contributions to be obtained from the FT.
    */
    cout <<"Number of harmonics set to: "<<nx<<" x "<<ny<<"\n";

    dimx=2*nx-1;
    dimy=2*ny-1;
}

/** Set the current wavelength
    @param l the wavelength to be set
*/
void structure::set_wavelength(double l)
{
    lambda=l;
    cout << "Wavelength set to: "<<lambda<<" m.\n";
}

/** Get the current wavelength
    @return the current wavelength
*/
double structure::get_wavelength()
{
    return lambda;
}

/** Increase the size of the section and interface vectors.

*/
void structure::ensureSize(int newsize)
{
    // Save the old pointers and size
    int old_size = actual_size;
    actual_size = newsize;
    section *oldsec = sec_list;
    finterface *oldil = interf_list;

    // Allocate the new pointers
    sec_list = new section[actual_size];
    interf_list = new finterface[actual_size];
    cout<<"ensureSize"<<endl;
    // Copy the old objects in the new position
    for (int i=0; i<old_size; ++i) {
        sec_list[i]=oldsec[i];
        interf_list[i]=oldil[i];
    }

    // Delete the old arrays
    delete[] oldsec;
    delete[] oldil;

    // Set the father of each interface
    for(int i=0; i<actual_size; ++i)
        sec_list[i].setFather(this);
}

double structure::getOrder(complex<double> e, double radius)
{
    complex<double> gg;
    double order;

    double freq=CELERITY/lambda;
    double omega=2.0*M_PI*freq;
    double k_0=omega/CELERITY;

    gg=e/k_0;

    if (radius!=0.0)
        order=2.0*M_PI*abs(radius)/lambda*gg.real();
    else
        order=2.0*M_PI/lambda*gg.real();

    return order;
}

double structure::getLosses_dB_cm(complex<double> e)
{
    return 20.0*e.imag()*log10(exp(1))/100.0;
}

double structure::getQuality(complex<double> e)
{
    complex<double> gg;
    double quality;

    double freq=CELERITY/lambda;
    double omega=2.0*M_PI*freq;
    double k_0=omega/CELERITY;

    gg=e/ (k_0);
    quality = (-gg.real()/gg.imag()/2.0);

    return quality;
}

complex<double> structure::getEffectiveIndex(complex<double> e)
{
    complex<double> g;

    double freq=CELERITY/lambda;
    double omega=2.0*M_PI*freq;
    double k_0=omega/CELERITY;

    g=e/ (k_0);
    return g;
}

/** Set the angles for a non-normal excitation wave.
    Angles in radians.
*/
void structure::setAngles(double n, double thetax, double thetay)
{
    double freq=CELERITY/lambda;
    double omega=2.0*M_PI*freq;
    double k_0=omega/CELERITY;

    ksinthetax=n*k_0*sin(thetax);
    ksinthetay=n*k_0*sin(thetay);
    
}

/** Get the current k*sin(thetax) term. */
double structure::getKsinthetax(void)
{
    return ksinthetax;
}


/** Get the current k*sin(thetay) term. */
double structure::getKsinthetay(void)
{
    return ksinthetay;
}

/** Write the given matrix in the Gnuplot format.

*/
void structure::fileoutputGP(db_matrix out, FILE *f, rimco_e rimc, double dx,
    double dy)
{
    int i,j;
    if(f==NULL)
        throw parsefile_commandError("Can not open output file :-(");

    if(rimc==C) {
        fprintf(f, "# x         y        neff.real       neff.imag\n");
    } else if(rimc==R) {
        fprintf(f, "# x         y        neff.real\n");

    } else if(rimc==I) {
        fprintf(f, "# x         y        neff.imag\n");
    } else if(rimc==M) {
        fprintf(f, "# x         y        abs(neff)\n");
    }
    for(i=0; i<out.getNrow(); ++i) {
        for(j=0; j<out.getNcol(); ++j) {
            if(rimc==C) {
                fprintf(f, "%le %le %le %le\n",
                    j*dx, i*dy,out(i,j).real(),out(i,j).imag());
            } else if(rimc==R) {
                fprintf(f, "%le %le %le\n",
                    j*dx, i*dy,out(i,j).real());
            } else if(rimc==I) {
                fprintf(f, "%le %le %le\n",
                    j*dx, i*dy,out(i,j).imag());
            } else if(rimc==M) {
                    fprintf(f, "%le %le %le\n",
                        j*dx, i*dy,abs(out(i,j)));
            }
        }
        fprintf(f,"\n");
    }
    fclose(f);
}

/** Write the given matrix in the Optiwave format.
*/
void structure::fileoutputOW(db_matrix out, FILE *f, double dx, double dy)
{
    int i,j;
    if(f==NULL)
        throw parsefile_commandError("Can not open output file :-(");
    fprintf(f, "BCF3DCX 3.0\n");
    fprintf(f, "%d %d\n", out.getNcol(), out.getNrow());
    fprintf(f, "0 %le 0 %le\n", out.getNcol()*dx, out.getNrow()*dy);
    for(i=0; i<out.getNrow(); ++i) {
        for(j=0; j<out.getNcol(); ++j) {
            fprintf(f, "%le, %le\n", out(i,j).real(),out(i,j).imag());
        }
    }
    fclose(f);
}

/** Calculate the radian frequency associated to the current wavelength.
*/
double structure::getOmega(void)
{
    double freq=CELERITY/lambda;
    double omega=2.0*M_PI*freq;

    return omega;
}

/** Get the total structure length
*/
double structure::getTotalLength(void)
{
    double totlen=0;
    for(int i=0; i<number_of_sections;++i)
        totlen+=sec_list[i].tot_z;

    return totlen;
}

/** Calculate the Bloch modes once the current structure is supposed to be
    the fundamental period of a periodic structure.
*/
vector<double> structure::calculateBlochModes()
{

    if(globalS.Rmp.getNrow()==0)
        throw parsefile_commandError("The S matrix should be created before"
        " attempting to calculate Bloch modes.");

    cout << "Creating matrices for generalized eigenvector "
        "decomposition."<<endl;

    db_matrix zero(globalS.Rmp.getNrow(), globalS.Rmp.getNcol());

    db_matrix one=db_matrix::createUnitMatrix(globalS.Rmp.getNrow(),
        globalS.Rmp.getNcol());

    int first_section=0;
    int last_section=number_of_sections-1;

    complex<double>e;
    double len;

    /*
    Developpements D. Bucci
avec "carpet"
Bloch-mode with the real part of r.i. closest to 1.582 is (1.58196,-0.00230479)
Bloch-mode with the real part of r.i. closest to 1.609 is (1.61436,5.6659)

sans "carpet"
Bloch-mode with the real part of r.i. closest to 1.582 is (1.58196,-0.00225252)
Bloch-mode with the real part of r.i. closest to 1.609 is (1.61436,5.6659)


*/
    db_matrix A = db_matrix::mergeMatrixQuad(
        sec_list[last_section].propagationMatrix()*globalS.Tpp,
        zero,
        sec_list[first_section].propagationMatrix()*globalS.Rmp,
        -one);


    /*db_matrix B = db_matrix::mergeMatrixQuad(
        one,
        -sec_list[last_section].propagationMatrix()*globalS.Rpm,
        zero,
        -sec_list[first_section].propagationMatrix()*
        globalS.Tmm);   */


    db_matrix toto1=sec_list[last_section].propagationMatrix();
    toto1.scalemult(complex<double>(-1.0,0.0), globalS.Rpm);
    db_matrix toto2=sec_list[first_section].propagationMatrix();
    toto2.scalemult(complex<double>(-1.0,0.0), globalS.Tmm);

    db_matrix B = db_matrix::mergeMatrixQuad(one, toto1,    zero, toto2);
    toto1.kill();
    toto2.kill();

    cout << "Launched the eigenvector calculation."<<endl;

    db_matrix eigval;

    A.eigGen(B,&eigval);

    bloch_eigvect=A;

    db_matrix d(A.getNcol(),1);

    bloch_eigval=d;

    vector<double> results;
    complex<double> g;

    for(int i=0; i<A.getNcol();++i) {
        g=eigval(i,i);
        bloch_eigval(i,0)=g;
        results.insert(results.end(),g.real());
        results.insert(results.end(),g.imag());
    }

    return results;
}

/** outBlochFile outputs all the calculated Bloch eigenvalues on a file
    @param pp the integer number of 2*pi to be considered and added for phase
        unwrapping
    @param name the file name.
*/

void structure::outBlochFile(double pp, char *name)
{
    double totlen = getTotalLength();
    double k0=(2.0*M_PI)/lambda;
    double realmin = 1e38;
    complex<double> calc;
    complex<double> closest=0;
    bool toPrint=false;
    FILE *fout=NULL;

    fout=fopen(name,"w");
    if(fout==NULL){
        throw parsefile_commandError("Can not create file.");
    }
    fprintf(fout,
        "#eig_re       eig_im         \t  n_eff_re       n_eff_im\n");

    for(int i=0; i<bloch_eigvect.getNcol();++i) {

        // Calculation of the phase
        calc=complex<double>(0,1)*(log(bloch_eigval(i,0))-
            2.0*M_PI*pp*complex<double>(0,1))/(k0*totlen);

        toPrint=false;
        fprintf(fout, "%le %le  \t  %le %le\n",
            bloch_eigval(i,0).real(),
            bloch_eigval(i,0).imag(),
            calc.real(),
            calc.imag());
    }
    fclose(fout);
    cout<<"File: "<<name<<" written."<<endl;
}

/** outBlochClosest returns the Bloch eigenvalue having an effective index
    closest to pa.
    @param pp the integer number of 2*pi to be considered and added for phase
        unwrapping.
    @param pa
    @return the effective index of the Bloch eigenvector whose real part is
        the closest to the given pa.
*/
complex<double> structure::outBlochClosest(double pp, complex<double> pa)
{
    double totlen = getTotalLength();
    double k0=(2.0*M_PI)/lambda;
    double realmin = 1e38;
    complex<double> calc;
    complex<double> closest=0;

    for(int i=0; i<bloch_eigvect.getNcol();++i) {
        // Calculation of the phase
        calc=complex<double>(0,1)*(log(bloch_eigval(i,0))-
            2.0*M_PI*pp*complex<double>(0,1))/(k0*totlen);

        if(realmin>abs(calc-pa)) {
            closest=calc;
            realmin = abs(closest-pa);
        }
    }
    return closest;
}

int structure::indexBlochClosest(double pp, double pa)
{
    double totlen = getTotalLength();
    double k0=(2.0*M_PI)/lambda;
    double realmin = 1e38;
    complex<double> calc;
    complex<double> closest=0;
    int imin=-1;

    for(int i=0; i<bloch_eigvect.getNcol();++i) {
        // Calculation of the phase
        calc=complex<double>(0,1)*(log(bloch_eigval(i,0))-
            2.0*M_PI*pp*complex<double>(0,1))/(k0*totlen);

        if(realmin>abs(calc.real()-pa)) {
            closest=calc;
            realmin = abs(closest.real()-pa);
            imin=i;
        }

    }
    cout<<"Selected effective index: "<<closest<<" index: "<<imin<<endl;

    return imin;

}

/**
    @param pp the integer number of 2*pi to be considered and added for phase
        unwrapping
*/
void structure::outBlochMaxImag(double pp, double pa)
{
    double totlen = getTotalLength();
    double k0=(2.0*M_PI)/lambda;
    double realmin = 1e38;
    complex<double> calc;
    complex<double> closest=0;
    bool toPrint=false;


    for(int i=0; i<bloch_eigvect.getNcol();++i) {

        // Calculation of the phase
        calc=complex<double>(0,1)*(log(bloch_eigval(i,0))-
            2.0*M_PI*pp*complex<double>(0,1))/(k0*totlen);

        if(abs(calc.imag())<pa) {
            cout << bloch_eigval(i,0)<< " ------> " << calc;
            cout <<endl;
        }
    }

}
