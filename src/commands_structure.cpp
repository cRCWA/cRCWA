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
#include <vector>

#include <unistd.h>

#include "block_matrix.h"
#include "parsefile.h"
#include "structure.h"
#include "section.h"
#include "finterface.h"
#include "commands.h"
#include "parallelize.h"


/***************************************************************************
*     Structure-level commands                                             *
****************************************************************************/


/*  SIZE: Define calculation window size.
    Usage:
    size sx sy

    Parameters:
    sx: total size of the calculation window
    sy: total size of the calculation window
*/
int commands::c_size(parsefile *obj, int argc,char *argv[])
{
    double sx;
    double sy;
    structure *p;
    p = wantedParameters(argv[0], argc, 3, obj);

    if(sscanf(argv[1], "%20lf", &sx)!=1 || sscanf(argv[2], "%20lf", &sy)!=1){
        cerr<<"size: error while reading parameters.\n";
        return ERROR;
    }else{
        p->do_size(sx, sy);
    }

    return 0;
}


/*  HARMONICS: Define the number of space harmonics to be used. We count the
    number of harmonics taken into account in the DFT of the input structure.
    The result will be obtained on half harmonics (not counting the zero term).

    Usage:
    harmonics nx ny

    Parameters:
    nx: number of harmonics on the x axis
    ny: number of harmonics on the y axis
*/
int commands::c_harmonics(parsefile *obj, int argc,char *argv[])
{
    int nx;
    int ny;
    structure *p;
    p = wantedParameters(argv[0], argc, 3, obj);

    if(sscanf(argv[1], "%20d", &nx)!=1 || sscanf(argv[2], "%20d", &ny)!=1){
        cerr<<"harmonics: error while reading parameters.\n";
        return ERROR;
    }else{
        p->do_harmonics(nx,ny);
    }

    return 0;
}


/*  SYMMETRY: Define  if the structure is symmetrical or not.

    Usage:
    symmetry sanx sany

    Parameters:
    sanx: Set up the mirror symmetry perpendicular of the x axis, can be
        {s,a,n}
    sany: Set up the mirror symmetry perpendicular of the y axis, can be
        {s,a,n}

    s stand for symmetrycal
    a stand for anti-symmetrical
    n stand for none
*/
int commands::c_symmetry(parsefile *obj, int argc,char *argv[])
{
    char sanx;
    char sany;
    symmetry_e symx;
    symmetry_e symy;
    structure *p;
    p = wantedParameters(argv[0], argc, 3, obj);

    if(sscanf(argv[1], "%c", &sanx)!=1 || sscanf(argv[2], "%c", &sany)!=1){
        cerr<<"symmetry: error while reading parameters.\n";
        return ERROR;
    }else{
        if ((sanx != 's'  && sanx != 'a' && sanx != 'n') ||
         (sany != 's'  && sany != 'a' && sany != 'n')){
            throw parsefile_commandError("symmetry: you should specify"
            " {s|a|n}");
        }

        // p contains the current object and should therefore be a structure

        switch (sanx){
            case 's':
                cout << "symmetry along x has been set, ";
                symx = symmetric;
                break;
            case 'n':
                cout << "no symmetry along x has been set, ";
                symx = no_symmetry;
                break;

            case 'a':
                cout << "anti-symmetry along x has been set, ";
                symx = anti_symmetric;
                break;
        }

        switch (sany) {
            case 's':
                cout << "symmetry along y has been set." << endl;
                symy = symmetric;
                break;
            case 'n':
                cout << "no symmetry along y has been set." << endl;
                symy = no_symmetry;
                break;

            case 'a':
                cout << "anti-symmetry along y has been set." << endl;
                symy = anti_symmetric;
                break;
        }

        // In principle, this should be solved in the future. There is to
        // adjust the construction of the W matrix (section::calcWmatrix
        // method).
        if (symx==symy && symx!=no_symmetry) {
            throw parsefile_commandError("symmetry: you can not (yet!) have"
                " the same symmetry in x and y.");
        }

        p->symx = symx;
        p->symy= symy;

/*      // Set the default matrix development strategy
        p->crtr = &p->NonDevSym;
        p->NonDevSym.setAlpha(-1.0);*/
    }

    return 0;
}

/*  WAVELENGTH: Define the wavelength to be used.

    Usage:
    wavelength lambda

    Parameters:
    lambda: the wavelength to be used
*/
int commands::c_wavelength(parsefile *obj, int argc,char *argv[])
{
    double l;
    structure *p;
    p = wantedParameters(argv[0], argc, 2, obj);

    if(sscanf(argv[1], "%20lf", &l)!=1){
        cerr<<"wavelength: error while reading parameter.\n";
         return ERROR;
    }else{
        // p contains the current object and should therefore be a structure

        p->set_wavelength(l);
        p->insertVar("ans",l);
    }

    return 0;
}

/*  ANGLES: Define excitation angles different from 0.

    Usage:
    angles n0 thetax thetay

    Parameters:
    n0: the refractive index of the excitation section to be used for the angle
        calculation
    thetax: the angle of the excitation wave in the xz plane, with respect to
        the propagation direction (in radians).
    thetay: the angle of the excitation wave in the yz plane, with respect to
        the propagation direction (in radians).
*/
int commands::c_angles(parsefile *obj, int argc,char *argv[])
{
    double n0,thetax,thetay;
    structure *p;
    p = wantedParameters(argv[0], argc, 4, obj);

    if(sscanf(argv[1], "%20lf", &n0)!=1){
        cerr<<"angles: error while reading n0.\n";
         return ERROR;
    }if(sscanf(argv[2], "%20lf", &thetax)!=1){
        cerr<<"angles: error while reading thetax.\n";
         return ERROR;
    } if(sscanf(argv[3], "%20lf", &thetay)!=1){
        cerr<<"angles: error while reading thetay.\n";
         return ERROR;
    } else{
        // p contains the current object and should therefore be a structure

        p->setAngles(n0,thetax, thetay);
        cout << "Angles set to: "<<thetax<<" rad and "<<thetay<<
            " rad in a section with refractive index "<<n0<<"\n";
    }

    return 0;
}


/** SOLVE: Find "interesting" modes. A mode is considered interesting if the
    real part of its effective index is comprised between the substrate index
    and the maximum value of the refractive index used in the structure.
    The substrate refractive index can be varied using the lowindex command.
    An alternative way of selecting interesting modes can be selected using the
    order command in bent waveguides. The software will calulate the azimuthal
    order and ensure it is in the interval specified. This is useful when
    working with bent resonating structures (micro disks, micro rings,
    micro spheres...).
    The calculations are iterated for all defined sections.


    Usage:
    solve

    Parameters:

    Notes:
    Everything should be defined...
*/
int commands::c_solve(parsefile *obj, int argc,char *argv[])
{

    structure *p;
    vector<double> results;
    cout<<endl;

    p = wantedParameters(argv[0], argc, 1, obj);
    if(!p->cur->isSubstrateSet) {
        throw parsefile_commandError("solve: The substrate should be set "
            "before trying to find modes.");
    }
    if (p->lambda==0) {
        cerr<<"solve: wavelength not defined\n";
        return ERROR;
    }
    if (p->tot_x==0 || p->tot_y==0) {
        cerr<<"solve: calculation window size not defined\n";
        return ERROR;
    }
    results=p->do_solve();
    double *s =new double[results.size()];
    for(int j=0; j<results.size();++j) {
        s[j]=results[j];
    }
    p->nP->insertArray("ans", results.size(), s);
    delete[] s;
    return 0;
}

/** WANTS: save info for later field calculations. This command should be
    used before using solve. This prevents solve eliminate matrices which will
    then be useful for the H calculation or for the propagation.

    Usage:
    wants h

    Parameters:
    h: can be either "Hx" or "Hy" or "both" or "propagation". Specifies which
    field component will then be requested by a call to outgmodes or if a
    propagation will be calculated. This command is important as afmm by
    default discards some matrices that become useful only if certain
    field components are calculated. Thus, by using wants, afmm is instructed
    to not to discard them.

*/
int commands::c_wants(parsefile *obj, int argc,char *argv[])
{
    structure *p;
    if(argc==2){
        if((p=dynamic_cast<structure *>(obj))==NULL)
            throw parsefile_commandError("Incorrect dynamic cast: "
                "programming error:-(");
        p->do_wants(argv[1]);
    } else {
        cerr<<"wants: invalid number of parameters.\n";
        return ERROR;
    }
    return 0;
}

/*  ASSEMBLE: prepares all the S matrices needed for the propagation.

    Usage:
    assemble

    Parameters:
    none

*/

int commands::c_assemble(parsefile *obj, int argc,char *argv[])
{
    structure *p;

    cout<<endl;

    if(argc==1){
        if((p=dynamic_cast<structure *>(obj))==NULL)
            throw parsefile_commandError("Incorrect dynamic cast: programming"
                " error:-(");
        p->do_assemble();
    } else {
        cerr<<"assemble: invalid number of parameters.\n";
        return ERROR;
    }
    return 0;

}

/** BLOCH: calculates the Bloch modes of the current structure, if it is
    periodized by an infinite repetition on the propagation axis.

    Usage:
    bloch
*/
int commands::c_bloch(parsefile *obj, int argc,char *argv[])
{
    structure *p;

    if(argc==1){
        if((p=dynamic_cast<structure *>(obj))==NULL)
            throw parsefile_commandError("Incorrect dynamic cast: programming"
                " error:-(");

        vector<double> results= p->do_bloch();

        double *s=new double[results.size()];
        for(int j=0; j<results.size();++j) {
            s[j]=results[j];
        }

        p->nP->insertArray("ans", results.size(), s);
        delete[] s;

        cout<<"Done."<<endl;
    } else {
        cerr<<"bloch: invalid number of parameters.\n";
        return ERROR;
    }
    return 0;
}

/** OUTBLOCH: output the results of the Bloch mode calculation (see the BLOCH
        command.

    Usage:
        outbloch phase selection d

    Parameters:
    phase       number of complete 2*pi contributions to be added for phase
                unwrapping
    selection   selection strategy:
        il      Imaginary lower than a given threshold.  The maximum absolute
                value of imaginary part for printing modes must be less than
                parameter.
        cr      Closest real. Just one mode closest to the specified parameter.
        fi      File. All effective indices as well as coefficients will be
                written in a file specified by the parameter d.
    d           depends on selection
*/
int commands::c_outbloch(parsefile *obj, int argc,char *argv[])
{
    structure *p;

    if(argc==4){
        if((p=dynamic_cast<structure *>(obj))==NULL)
            throw parsefile_commandError("Incorrect dynamic cast: programming"
                " error:-(");

        if (p->bloch_eigvect.isEmpty()) {
            throw parsefile_commandError(
                "outbloch: you should calculate the Bloch modes of the"
                " structure by using the bloch command, before calling"
                " outbloch.");
        }

        double mi=0;
        char *name=NULL;

        double pp=0;
        if(sscanf(argv[1], "%20lf", &pp)!=1) {
            throw parsefile_commandError(
                "outbloch: I can not read the parameter for phase"
                " unwrapping.");
        }

        if(strcmp(argv[2],"il")==0) {
            //type = structure::MAX_IMAGINARY;
            if(sscanf(argv[3], "%20lf", &mi)!=1) {
                throw parsefile_commandError(
                    "outbloch: I can not read the parameter for imaginary"
                    " part selection.");
            }
            p->outBlochMaxImag(pp, mi);
        } else if(strcmp(argv[2],"cr")==0) {
            //type = structure::CLOSEST_REAL;
            if(sscanf(argv[3], "%20lf", &mi)!=1) {
                throw parsefile_commandError(
                    "outbloch: I can not read the parameter for closest real"
                    " part selection.");
            }
            cout << "Bloch-mode with the real part of r.i. closest to "
              << mi << " is "<< p->outBlochClosest(pp, complex<double>(mi,0))
              <<endl;

        } else if(strcmp(argv[2],"ci")==0) {
            //type = structure::CLOSEST_IMAGINARY;
            if(sscanf(argv[3], "%20lf", &mi)!=1) {
                throw parsefile_commandError(
                    "outbloch: I can not read the parameter for closest real"
                    " part selection.");
            }
            cout << "Bloch-mode with the real part of r.i. closest to "
              << mi << " is "<< p->outBlochClosest(pp, complex<double>(0, mi))
              <<endl;

        } else if(strcmp(argv[2],"fi")==0) {
            //type = structure::FILE_SEL;
            name=argv[3];
            p->outBlochFile(pp, name);
        } else {
            throw parsefile_commandError(
                    "outbloch: Unrecognized output strategy.");
        }

    } else {
        cerr<<"outbloch: invalid number of parameters.\n";
        return ERROR;
    }
    return 0;
}

/**  EXCITATION: defines the excitation of the calculation region.

    Usage:
    excitation fb t r i a b

    Parameters:
    fb      can be {f|b} where f indicates the forward excitation where b
            stands for backward excitation. The forward excitation is applied
            to the beginning of the calculation region in section 0 and
            propagates towards the last section defined.
            The backward excitation is applied to the end of the calculation
            region, i.e. at the right of the last section which has been
            defined and propagates towards the left.
    t       defines the excitation types and can be {m|s|fx|fy|cx|cy|bl}.
            The letter m
            selects the fundamental mode exciation, which is the mode
            recognized as "interesting" which has the highest real part for
            the effective index. In this case, the parameter a is ignored.
            The letter s stands for the selection of the "interesting"
            eigenmode having the real part of the effective index the closest
            possible to the value indicated by the parameter a.
            If t=fx, the excitation field is read from a file and considered
            as the Ex electric field component. If i=fy, the program consider
            the file as the Ey component. When the field is read from a file,
            a represent the multiplication coefficient (see indfile command)
            whereas b is the file name. The size specified in the excitation
            file must match with the size of the calculation window.
            cx, cy means that a plane wave field excitation is provided.
            In this case, if a and b are provided, it will be considered as
            the indexes to be used for the angle (see the manual for details).
            bl indicates a Bloch-mode excitation (see the bloch command).
            The Bloch eigenmodes should have been already calculated and the
            wanted one is selected via the closest real part to the specified
            b parameter, while a is the phase unwrapping parameter (see
            outbloch for more information about this detail).
            t=id means that the mode index should be given (see the modeindex
            command).
    r,i     those parameters give the real and imaginary part of the excitation
            of the selected mode. This is useful to control the phase
            relation between two modes or the forward and backward excitation.
    a,b     additional parameters whose signification depends to the
            parameter t which has been provided to the command.

*/
int commands::c_excitation(parsefile *obj, int argc,char *argv[])
{
    enum exc_d {forward, backward} direction;
    enum exc_t {modal, select, file, constant, bloch, index} type;

    bool do_the_shift = false;

    structure *p;

    // This is the number of the section to which the excitation will be
    // applied
    int sec=0;
    db_matrix *exc;

    if(argc>4){
        if((p=dynamic_cast<structure *>(obj))==NULL)
            throw parsefile_commandError("Incorrect dynamic cast: programming"
                " error:-(");

        // If an excitation command is executed, the excitation vectors at
        // each interface MUST be calculated during the upcoming propagation
        // command.
        p->excitationCalculated=false;

        // We first begin by analyzing the command parameters given by the user
        if(strcmp(argv[1],"f")==0) {
            direction = forward;
        } else if(strcmp(argv[1],"b")==0) {
            direction = backward;
        } else {
            throw parsefile_commandError("excitation: unrecognized direction"
                " of the excitation. It must be {b|f}.\n");
        }

        if(strcmp(argv[2],"m")==0) {
            type = modal;
        } else if(strcmp(argv[2],"s")==0) {
            type = select;
        } else if(strcmp(argv[2],"fx")==0) {
            type = file;
            do_the_shift = false;
        } else if(strcmp(argv[2],"fy")==0) {
            type = file;
            do_the_shift = true;
        } else if(strcmp(argv[2],"cx")==0) {
            type = constant;
            do_the_shift = false;
        } else if(strcmp(argv[2],"cy")==0) {
            type = constant;
            do_the_shift = true;
        } else if(strcmp(argv[2],"bl")==0) {
            type = bloch;
            do_the_shift = false;
        } else if(strcmp(argv[2],"id")==0) {
            type = index;
            do_the_shift = false;
        } else {
            throw parsefile_commandError("excitation: output type. It must be"
                " {m|s|fx|fy|cx|cy|bl|id}.\n");
        }

        // Basically, we just needs to know in which section and on which
        // vector we are working to determine how the excitation will be made.

        if(direction == forward) {
            sec = 0;
            exc = &p->sec_list[sec].sWp;
            cout << "Forward excitation:\n";
        } else if(direction == backward) {
            sec = p->number_of_sections-1;
            exc = &p->sec_list[sec].sWm;
            cout << "Backward excitation:\n";

        } else {
            throw parsefile_commandError("excitation: unrecognized direction"
                " of the excitation. Programming erroror incomplete???\n");
        }

        complex<double> phase;
        double re, im;

        // Here we read the real and imaginary part to be used as the
        // coefficient for the phase of the excitation.

        if(sscanf(argv[3], "%20lf", &re)!=1) {
            throw parsefile_commandError("excitation: can not read the real"
                " part of the excitation coefficient.\n");
        }
        if(sscanf(argv[4], "%20lf", &im)!=1) {
            throw parsefile_commandError("excitation: can not read the"
                " imaginary part of the excitation coefficient.\n");
        }

        phase = complex<double>(re, im);

        // Here we fabricate the excitation.
        if(type == modal) {
            *exc = (p->sec_list[sec]).create_excitation(
                section::SELECT_FUNDAMENTAL,0.0, phase,0.0,0.0);
        } else if (type == constant) {
            int index_xz=0;
            int index_yz=0;
            if(argc>5 && (sscanf(argv[5], "%20d", &index_xz)!=1
                || sscanf(argv[6], "%20d", &index_yz)!=1)) {
               throw parsefile_commandError("excitation: can not read the"
                " specified angle indexes.\n");
            }
            if (do_the_shift)
                *exc = (p->sec_list[sec]).create_excitation_const(true,
                    0.0, phase,index_yz,index_xz); // exc. for Ey
            else
                *exc = (p->sec_list[sec]).create_excitation_const(false,
                    0.0, phase,index_yz,index_xz); // exc. for Ex
        } else if (type == select) {
            // If we are in the select mode, we should first read the
            // parameter a as a double precision constant and then use it for
            // searching of the modes.
            double value;
            if (argc<6) {
                throw parsefile_commandError("excitation: invalid number of"
                    " parameters.\n");
            }
            if(sscanf(argv[5], "%20lf", &value)!=1) {
                throw parsefile_commandError("excitation: error while reading"
                    " the parameter 'a'.\n");
            }

            db_matrix newexc = (p->sec_list[sec]).create_excitation(
                section::SELECT_REAL_INDEX,value, phase,0.0,0.0);

            // When the excitation already exists, the new requested excitation
            // Is added to the previous one, allowing us to excite multiple
            // modes

            if ((*exc).getNcol() > 0) {
                *exc += newexc;
                cout << "The mode excitation is added to the previous"
                        " existing excitation\n";
            } else {
                *exc = newexc;
            }
        } else if(type==file) {

            // This version of the command requires 4 parameters.
            if (argc<5) {
                throw parsefile_commandError("excitation: invalid number of"
                    " parameters.\n");
            }

            double mult;
            // If we are in the select mode, we should first read the
            // parameter a as a double precision constant and then use it for
            // the multiplicator.

            if(sscanf(argv[5], "%20lf", &mult)!=1) {
                throw parsefile_commandError("excitation: error while reading"
                    " the parameter 'a'.\n");
            }

            // In this case, the excitation will be read from a file specified
            // by the last parameter

            double tot_x;
            double tot_y;

            // The file is in the Optiwave f3d format.
            db_matrix rd = read_file(argv[6], tot_x, tot_y, mult);

            // The sizes specified in the file must match with the size of the
            // calculation window. Eps is the tolerance to be used to compare
            // the two double precision values.

            double eps = 1e-15;

            if(abs(tot_x -p->tot_x)>eps*tot_x ||
                abs(tot_y- p->tot_y)>eps*tot_y) {
                cerr.precision(17);
                cerr<<tot_x<<"!="<<p->tot_x<<endl;
                cerr<<tot_y<<"!="<<p->tot_y<<endl;
                throw parsefile_commandError("excitation: the size specified"
                    " in the excitation file must match with the size of the"
                    " calculation window.\n");
            }

           *exc = (p->sec_list[sec]).create_excitation_from_file(
            do_the_shift, rd);

           cout << "File " << argv[6] <<" read correctly.\n";
           cout.flush();

        } else if (type==bloch) {
            // If we are in the bloch mode, we should first read the
            // parameter a as a double precision constant and then use it for
            // searching of the modes.
            double phase_unwr;
            double re;
            if (argc<7) {
                throw parsefile_commandError("excitation: invalid number of"
                    " parameters.\n");
            }
            if(sscanf(argv[5], "%20lf", &phase_unwr)!=1) {
                throw parsefile_commandError("excitation: error while reading"
                    " the parameter 'a' (phase unwrapping).\n");
            }
            if(sscanf(argv[6], "%20lf", &re)!=1) {
                throw parsefile_commandError("excitation: error while reading"
                    " the parameter 'b' (real part of neff).\n");
            }
            if(p->bloch_eigvect.isEmpty()) {
                throw parsefile_commandError("excitation: you should first use"
                    " the bloch command to calculate Bloch-mode of the"
                    " structure.\n");
            }
            db_matrix column_front(p->bloch_eigvect.getNrow()/2,1);
            db_matrix column_back(p->bloch_eigvect.getNrow()/2,1);
            int shift_f=p->bloch_eigvect.getNrow()/2;

            int i=p->indexBlochClosest(phase_unwr, re);

            for(int j=0; j<p->bloch_eigvect.getNrow()/2; ++j) {
                column_front(j,0)=p->bloch_eigvect(j,i);
                column_back(j,0)=p->bloch_eigvect(j+shift_f,i)*
                    p->bloch_eigval(i,0);
            }

            p->sec_list[0].sWp = column_front;
            p->sec_list[p->number_of_sections-1].sWm = column_back;
        } else if(type==index) {
            // If we are in the select mode, we should first read the
            // parameter a as an integer constant and then use it for
            // the mode index
            int idx;
            if (argc<6) {
                throw parsefile_commandError("excitation: invalid number of"
                    " parameters.\n");
            }
            if(sscanf(argv[5], "%20d", &idx)!=1) {
                throw parsefile_commandError("excitation: error while reading"
                    " the parameter 'a' (mode index).\n");
            }
            // test jerome:

            db_matrix newexc = (p->sec_list[sec]).create_excitation(
                section::SELECT_MODE_POSITION,(double)idx, phase,0.0,0.0);

            // When the excitation already exists, the new requested excitation
            // Is added to the previous one, allowing us to excite multiple
            // modes

            if ((*exc).getNcol() > 0) {
                *exc += newexc;
                cout << "The mode excitation is added to the previous"
                        " existing excitation\n";
            } else {
                *exc = newexc;
            }
        } else {
            throw parsefile_commandError("excitation: unrecognized excitation"
                " type. Programming error or incomplete???\n");
        }
    } else {
        cerr<<"excitation: invalid number of parameters.\n";
        return ERROR;
    }
    return 0;
}

/** OUTDATA: set up for the calculation of additional output data. This command
    should be used before the propagation command.

    Usage:
    outdata t par1 par2 ...

    Parameters:
    t   type of the output data which should be recorded
        if t==i integral of the field over a rectangular window. The field
        component which will be integrated depends on which field component
        will constitute the output for the propagation command.

        outdata i cx cy w h filename

        where:
            cx  is the x center of the window to be used for the calculation
            cy  is the y center of the window to be used for the calculation
            w is the width of the window
            h is the height of the window
            filename is the name of the file on which the results should be
                recorded

        if t == g compute the 2D generation rate cylindrincally interpolated
        from the 3D generation rate.

        outdata g z0 z1 r1 dr filename

        where:
            dr is the step for the radius
            filename is the name of the file on which the results should be
                recorded
            it will be computed from z0 to z1 and from r=0 to r=r1.

        NOTE that the step in z in given by the command propagation.
*/
int commands::c_outdata(parsefile *obj, int argc,char *argv[])
{
    structure *p;

    if(argc>2){
        if((p=dynamic_cast<structure *>(obj))==NULL)
            throw parsefile_commandError("Incorrect dynamic cast: programming"
                " error:-(");

        // We first begin by analyzing the command parameters given by the user
        // which defines the type of additional output data needed

        if(strcmp(argv[1],"i")==0) {
            if(argc<7) {
                throw parsefile_commandError("outdata i: some parameters are"
                    " missing.\n");
            }

            p->additional_output_data.should_record_integral=true;

            if(sscanf(argv[2], "%20lf", &p->additional_output_data.
                window_centerx)!=1){
                throw parsefile_commandError("outdata: could not read the x"
                    " center of the window to be used for the calculation.\n");
            }
            if(sscanf(argv[3], "%20lf", &p->additional_output_data.
                window_centery)!=1){
                throw parsefile_commandError("outdata: could not read the y"
                    " center of the window to be used for the calculation.\n");
            }
            if(sscanf(argv[4], "%20lf", &p->additional_output_data.
                window_width)!=1){
                throw parsefile_commandError("outdata: could not read the"
                    " width of the window to be used for the calculation.\n");
            }
            if(sscanf(argv[5], "%20lf", &p->additional_output_data.
                window_height)!=1){
                throw parsefile_commandError("outdata: could not read the"
                    " height of the window to be used for the calculation.\n");
            }
            p->additional_output_data.window_file_name=argv[6];

            cout << "Additional output: record the field integral in the next"
                " propagation command.\n";
            cout << "Window centered in "<<p->additional_output_data.
                window_centerx<<" m, "<<
                p->additional_output_data.window_centery<<" m, size "<<
                p->additional_output_data.window_width<<
                " m x "<<p->additional_output_data.window_height<<" m.\n";


        } else if(strcmp(argv[1],"g")==0) {
            if(argc<7) {
                throw parsefile_commandError("outdata g: some parameters are"
                    " missing.\n");
            }

            p->additional_output_data.should_record_generation_rate=true;

            if(sscanf(argv[2], "%20lf", &p->additional_output_data.
                generation_from_z0)!=1){
                throw parsefile_commandError("outdata: could not read the z0 "
                    " coordinate for the generation rate calculation.\n");
            }

            if(sscanf(argv[3], "%20lf", &p->additional_output_data.
                generation_to_z1)!=1){
                throw parsefile_commandError("outdata: could not read the z1"
                    " coordinate for the generation rate calculation.\n");
            }

            if(sscanf(argv[4], "%20lf", &p->additional_output_data.
                generation_to_r1)!=1){
                throw parsefile_commandError("outdata: could not read the"
                    " radius for the generation rate calculation.\n");
            }

            if(sscanf(argv[5], "%20lf", &p->additional_output_data.
                generation_radial_step_size)!=1){
                throw parsefile_commandError("outdata: could not read the"
                   " radial step size for the generation rate calculation.\n");
            }

            p->additional_output_data.generation_file_name=argv[6];


            if(p->additional_output_data.generation_to_r1 > p->tot_x ||
                    p->additional_output_data.generation_to_r1 > p->tot_y){
                throw parsefile_commandError("outdata: it is not possible to"
                " use 'outdata g' if the the computational window is smaller"
                " than the radius.\n");
            }

        }else {
            throw parsefile_commandError("outdata: unrecognized parameter"
                " type.\n");
        }

        return 0;
    } else {
        throw parsefile_commandError("outdata: some parameters are needed.\n");
    }
}

/** CARPET sweep some convergence problems under the carpet: Force all modes to
    have a negative imaginary part for the effective index. This command should
    be used before SOLVE.

    Usage:
    CARPET

    Parameters:

*/
int commands::c_carpet(parsefile *obj, int argc,char *argv[])
{
    structure *p;
    p = wantedParameters(argv[0], argc, 1, obj);

    p->set_ensureConvergence(true);
    return 0;

}
