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

#include <iostream>
#include <iomanip>
#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "block_matrix.h"
#include "parsefile.h"
#include "structure.h"
#include "section.h"
#include "finterface.h"
#include "commands.h"


/* Version history


 [1] J.P. Hugonin et al. "Fourier modal methods for modelling optical
    dielectric waveguides," Optical and quantum electronics (2005) 37:107-119

 [2] J.P. Hugonin, P. Lalanne "Perfectly matched layers as nonlinear coordinate
     transforms: a generalized formalization," J. Opt. Soc. Am. A, vol. 22, n.
     9, september 2005

 [3] B. Martin "Etude et réalisation d'un spectromètre compact en optique
     intégrée sur verre," PhD report, INP Grenoble, 23 janvier 2009

 [4] D. Bucci, B. Martin, A. Morand "Application of the three-dimensional
    aperiodic Fourier modal method using arc elements in curvilinear
    coordinates," J. Opt. Soc. Am. A 29, 367-373 (2012)

*/
commands::commands(structure &s, numParser *p): parsefile(p)
{
    init(s);
}


// Register all commands in the hash table
void commands::init(structure &s)
{
    // Added in v. 1.0
    insert_command("#",&s,&commands::c_diesis);
    insert_command("substrate",&s,&commands::c_substrate);
    insert_command("size",&s,&commands::c_size);
    insert_command("harmonics",&s,&commands::c_harmonics);
    insert_command("wavelength",&s,&commands::c_wavelength);
    insert_command("rectangle",&s,&commands::c_rectangle);
    insert_command("pml",&s,&commands::c_pml);
    insert_command("solve",&s,&commands::c_solve);

    // Added in v. 1.0.2
    insert_command("inpstruct",&s,&commands::c_inpstruct);
    insert_command("outgmodes",&s,&commands::c_outgmodes);

    // Added in v. 1.0.3
    insert_command("lowindex",&s,&commands::c_lowindex);
    insert_command("highindex",&s,&commands::c_highindex);

    // Added in v. 1.1.0
    insert_command("pml_transf",&s,&commands::c_pml_transf);
    insert_command("indfile",&s,&commands::c_indfile);

    // Added in v. 1.1.0 (bent waveguides)
    insert_command("bend",&s,&commands::c_bend);

    // Added in v. 1.2.2 (H field calculation)
    insert_command("wants",&s,&commands::c_wants);

    // Added in v. 1.2.3
    insert_command("spectrum",&s,&commands::c_spectrum);

    // Added in v. 1.3 (propagation)
    insert_command("propagation",&s,&commands::c_propagation);
    insert_command("assemble",&s,&commands::c_assemble);
    insert_command("section",&s,&commands::c_section);
    insert_command("order",&s,&commands::c_order);
    insert_command("excitation", &s, &commands::c_excitation);
    insert_command("quit", &s, &commands::c_quit);
    insert_command("help", &s, &commands::c_help);
    insert_command("coefficient", &s, &commands::c_coefficient);
    insert_command("select", &s, &commands::c_select);
    insert_command("print", &s, &commands::c_print);

    // Added in v. 1.3.1
    insert_command("carpet",&s,&commands::c_carpet);
    insert_command("matdev",&s,&commands::c_matdev);

    // Added in v. 1.3.2
    insert_command("clear",&s,&commands::c_clear);

    // Added in v. 1.3.3
    insert_command("norm",&s,&commands::c_norm);
    insert_command("outdata",&s,&commands::c_outdata);
    // added also the system command in the parsefile class

    // Added in v. 1.3.4
    // label and goto in the parsefile class

    // Added in v. 1.3.5
    insert_command("eigenvc",&s,&commands::c_eigenvc);
    insert_command("eigenam",&s,&commands::c_eigenam);
    insert_command("memocc",&s,&commands::c_memocc);

    // Added in v. 1.4.1
    insert_command("power",&s,&commands::c_power);
    insert_command("powerZ",&s,&commands::c_powerZ);
    insert_command("monitor",&s,&commands::c_monitor);
    insert_command("symmetry",&s,&commands::c_symmetry);
    insert_command("let",&s,&commands::c_diesis);

    // Added in v. 1.4.2
    insert_command("bloch",&s,&commands::c_bloch);
    insert_command("outbloch",&s,&commands::c_outbloch);

    // Added in v. 1.4.3
    insert_command("parallel",&s,&commands::c_parallel);

     // Added in v. 1.4.4
    insert_command("selmodes",&s,&commands::c_selmodes);

    // Added in v. 1.4.5
    insert_command("modepos",&s,&commands::c_modepos);

    // Added in v. 1.4.6
    insert_command("angles",&s,&commands::c_angles);
    insert_command("powerz",&s,&commands::c_powerZ); // accept version with "z"


    // When adding new commands, please remember to update c_help

}

/*  HELP: show an help

    Usage:
    help

    Parameters:
    none
*/

int commands::c_help(parsefile *obj, int argc,char *argv[])
{
    cout << "I am AFMM, Aperiodic Fourier Modal Method full vectorial 3D propagation. \n"
        << "I can calculate propagation modes of a waveguide as well as the\n"
        << "electromagnetic field propagation through a 3D structure represented\n"
        << "using Fourier series.\n\n"
        << "Here is a list of the commands I can understand. Refer to the manual\n"
        << "for a detailed description of each of them.\n"
        << "\n"
        << " General or structure-level definition commands:\n"
        << "\n"
        << "  size         Define the size of the calculation window.\n"
        << "  harmonics    Define the number of harmonics to be taken.\n"
        << "  wavelength   Define the working wavelength.\n"
        << "  solve        Solve for eigenmodes for all sections.\n"
        << "  wants        Declare that the user wants to calculate H field or propagation.\n"
        << "  propagation  Propagate the field in the structure.\n"
        << "  assemble     Assemble the different sertions to create a structure.\n"
        << "  excitation   Specify the field excitation of the structure.\n"
        << "  carpet       Sweep some numerical dust under the carpet.\n"
        << "  outdata      Set up for the calculation of additional output data.\n"
        << "  symmetry     Set up if the structure is symmetrical.\n"

        << "\n"
        << "Section-level commands:\n"
        << "  substrate    Define a substrate of a given refractive index.\n"
        << "  rectangle    Define a rectangular waveguide.\n"
        << "  pml          Enter an anisotropic PML.\n"
        << "  pml_transf   Enter a coordinate transform PML region.\n"
        << "  indfile      Read refractive index distribution from an input file.\n"
        << "  inpstruct    Write on a file the current index distribution.\n"
        << "  outgmodes    Write on files the interesting modes found.\n"
        << "  lowindex     Define the lower bounds for refractive index.\n"
        << "  highindex    Define the upper bounds for refractive index.\n"
        << "  bend         Define the bending radius of the center of the section.\n"
        << "  spectrum     Write on a file the effective indices of all modes found.\n"
        << "  section      Create a section of the given length.\n"
        << "  order        Specify the modal search based on azimuthal order.\n"
        << "  coefficient  Get the mode coefficient on the current section.\n"
        << "  select       Select the current section.\n"
        << "  matdev       Specify the wanted strategy for the matrix developments.\n"
        << "  norm         Calculates the L2 norm of the excitation vector.\n"
        << "  eigenvc      Writes on a file all the modal excitation coefficients.\n"
        << "  eigenam      Write the eigenvector matrix of the current section.\n"
        << "  power        Calculate the power of forward or backward wave of the current section.\n"
        << "  powerZ       Calculate the total power at a given depth.\n"
        << "  monitor      Calculate the total power at a given depth and for a given surface.\n"
        << "  bloch        Calculate Bloch modes of the periodicized structure.\n"
        << "  outbloch     Output the  Bloch modes calculated with bloch.\n"
        << "  selmodes     Write modes using a certain selection rule.\n"
        << "  modepos      Load 'ans' with the positions of interesting modes.\n"

        << "\n"
        << "Other general usage commands:\n"
        << "  #            Comment.\n"
        << "  if           |\n"
        << "  else         | Conditional test.\n"
        << "  endif        |\n"
        << "  help         This help.\n"
        << "  load         Open and execute all commands contained in a file.\n"
        << "  quit         Quit immediately from AFMM.\n"
        << "  clear        Clear all matrices in memory.\n"
        << "  system       Execute a system command (must be activated with -e option).\n"
        << "  label        Create a label.\n"
        << "  goto         Jump to a label.\n"
        << "  for... next  Cycle.\n"
        << "  memocc       Print some data about the current memory occupation.\n"
        << "  print        Print a message.\n"
        << "  let          Evaluate expressions.\n"
        << "  parallel     Launch multiple threads in parallel (when possible).\n";

    return 0;
}

/* # comment */
int commands::c_diesis(parsefile *obj, int argc,char *argv[])
{ /* do nothing */  return 0; }



/*  QUIT: exit immediately. This is useful in the interactive mode as well
    when debugging a script.

    Usage:
    quit

    Parameters:
    none
*/

int commands::c_quit(parsefile *obj, int argc,char *argv[])
{
    structure* p = wantedParameters(argv[0], argc, 1, obj);
    throw parsefile_stop();

    return 0;
}

/** CLEAR Reset the structure, the wavelength and all definitions, except
    variables.

    Usage:
    clear
*/
int commands::c_clear(parsefile *obj, int argc,char *argv[])
{
    structure *p;
    p = wantedParameters(argv[0], argc, 1, obj);
    p->reset();
    cout<<"The structure definition has been cleared.\n";
    db_matrix::leaks();
    return 0;
}

/** MEMOCC Print some info about the current memory occupation

    Usage:
    memocc
*/
int commands::c_memocc(parsefile *obj, int argc,char *argv[])
{
    structure *p;
    section *q;

    p = wantedParameters(argv[0], argc, 1, obj);
    long wsize;
    long vsize, bsize;
    long totalsize=0;
    int nr, nc;
    int xx, yy;
    int fftsize;
    int excitation;
    int otherdata;
    int bloch;

    cout <<"------------------------- Memory report --------------------------"
         << endl;

    for(int i=0; i<p->number_of_sections; ++i) {
        q= (p->sec_list+i);
        // q contains the current section
        cout<<"Now considering section "<<(i+1)<<" of "<<
            p->number_of_sections << ".\n";

        wsize = q->W.getSize();
        vsize = q->V.getSize();
        bsize = q->B.getSize();
        totalsize+=wsize+vsize+bsize;

        nr=q->W.getNrow();
        nc=q->W.getNcol();

        cout<<"Size of the W matrix ("
            <<nr<<"x"<<nc<<"): "<<wsize<<" bytes."<<endl;
        nr=q->V.getNrow();
        nc=q->V.getNcol();

        cout<<"Size of the V matrix ("
            <<nr<<"x"<<nc<<"): "<<vsize<<" bytes."<<endl;

        nr=q->B.getNrow();
        nc=q->B.getNcol();
        cout<<"Size of the B matrix ("
            <<nr<<"x"<<nc<<"): "<<bsize<<" bytes."<<endl;



        xx = q->X1.getSize();
        xx+= q->X2.getSize();
        xx+= q->X3.getSize();
        xx+= q->X4.getSize();
        totalsize+=xx;

        yy = q->Y1.getSize();
        yy+= q->Y2.getSize();
        yy+= q->Y3.getSize();
        yy+= q->Y4.getSize();
        totalsize+=yy;

        cout<<"Size of the stored X matrices: "<<xx<<" bytes."<<endl;
        cout<<"Size of the stored Y matrices: "<<yy<<" bytes."<<endl;

        fftsize = q->M_fft.getSize();
        fftsize+= q->N_fft.getSize();
        fftsize+= q->O_fft.getSize();
        fftsize+= q->P_fft.getSize();
        fftsize+= q->Q_fft.getSize();
        fftsize+= q->R_fft.getSize();

        fftsize+= q->Mm1_fft.getSize();
        fftsize+= q->Nm1_fft.getSize();
        fftsize+= q->Om1_fft.getSize();
        fftsize+= q->Pm1_fft.getSize();
        fftsize+= q->Qm1_fft.getSize();
        fftsize+= q->Rm1_fft.getSize();

        totalsize+=fftsize;

        cout<<"Coeff. of permittivity/permeability: "<<fftsize
            <<" bytes."<<endl;

        excitation = q->sWp.getSize();
        excitation+= q->sWm.getSize();

        totalsize+=excitation;
        cout<<"Modal excitation coefficients: "<<excitation<<" bytes."<<endl;

        otherdata = sizeof(section);
        totalsize+=otherdata;

        cout<<"Other stored data: "<<otherdata<< " bytes."<<endl;

    }
    cout <<"------------------------------------------------------------------"
         << endl;

    bloch=q->father->bloch_eigvect.getSize();
    bloch+=q->father->bloch_eigval.getSize();

    totalsize+=bloch;
    cout<<"Bloch eigenvalues and eigenvectors: "<<bloch<<" bytes."<<endl;

    cout <<"------------------------------------------------------------------"
         << endl;

    cout<<"Total memory (approximate): "<<totalsize<<" bytes."<<endl;
    cout<< setprecision(2);
    cout<<"                         : "<<(float)totalsize/(1024*1024)<<
        " MiB."<<endl;
    cout<< setprecision(6);
    cout <<"------------------------------------------------------------------"
        <<endl;
    cout<<"Note: AFMM might have required more memory for storing some\n"
          "intermediate results when matrices are being created. \nA more"
          " realistic estimation might give a peak"
          " memory occupation\nthat can exceed by 200% what shown"
          " here."<<endl;
    cout <<"----------------------- End of Memory report ---------------------"
         << endl;
    return 0;
}

/* LET evaluate expressions
    All the expressions will be evaluated thanks to the calculator included
    in the parser.

*/
int commands::c_let(parsefile *obj, int argc,char *argv[])
{ /* do nothing */  return 0; }


/** Check that the number of parameters is correct and convert the generic
    pointer to a structure pointer

    @param cname the name of the command
    @param argc the number of parameters received
    @param w the wanted number of parameters
    @param obj the generic pointer to be converted
*/
structure * commands::wantedParameters(char *cname,
    int argc, int w, parsefile *obj)
{
    structure *p;
    if(argc!=w) {
        cerr<<cname;
        throw parsefile_commandError(": invalid number of parameters.");
    }
    if((p=dynamic_cast<structure *>(obj))==NULL)
        throw parsefile_commandError("Incorrect dynamic cast: programming"
            " error:-(");
    return p;
}

/** Read the contents of a file in the Optiwave format and store it in a
    db_matrix which will have the number of points defined in the file.
    The order of the point in the file is row-based:
    (x0,y0)
    (x1,y0)
    (x2,y0)
    (x3,y0)
    (x0,y1)
    (x1,y1)
    (x2,y1)
    (x3,y1)
    (x0,y2)
    ...

    This routine can be used for both .rid (refractive index) and .f3d file
    formats, since they are almost identical, except for the header.

    @param filename the path and the name of the file to be opened.
    @param tot_x a reference to the width of the calculation window. It will be
        changed by the call to read_file and the new value will be the size
        indicated in the file.
    @param tot_y a reference to the height of the calculation window. It will
        be changed by the call to read_file and the new value will be the size
        indicated in the file.
    @param mult the multiplier to be used for all sizes specified in the file
        (this is useful since Optiwave often specifies all dimensions in
        micrometers whereas AFMM uses meters).


    @return the db_matrix containing the points specified in the file.
*/
db_matrix commands::read_file(char* filename, double &tot_x, double &tot_y,
    double mult)
{
    int nx, ny;
    double xmin,xmax;
    double ymin,ymax;
    int i,j;

    FILE *f=fopen(filename,"r");
    if(f==NULL)
        throw parsefile_commandError("Can not open input file:-(");

    // skip a line

    fscanf(f, "%*[^\n]");
    if(fscanf(f,"%20d%20d",&nx,&ny)!=2) {
        fclose(f);
        throw parsefile_commandError("Unable to read the number of x and y"
            " points in the given file.");
    }
    if(fscanf(f,"%lf%lf%lf%lf",&xmin,&xmax, &ymin, &ymax)!=4) {
        fclose(f);
        throw parsefile_commandError("Unable to read the x and y ranges in the"
            " given file.");
    }
    tot_x=(xmax-xmin)*mult;
    tot_y=(ymax-ymin)*mult;

    cout << "The file contains "<<nx<< " x  " << ny << " points.\n";
    cout << "The x range is "<<tot_x << " m.\n";
    cout << "The y range is "<<tot_y << " m.\n";

    db_matrix lect(ny, nx);

    double re, im;
    double norm = 0;

    for(i=0; i<ny; ++i) {
        for(j=0; j<nx; ++j) {
            if(fscanf(f,"%lf,%lf",&re,&im)!=2) {
                fclose(f);
                throw parsefile_commandError("Unable to read a value in the"
                    " input file.");
            }
            lect(i,j)=complex<double>(re,im);
            norm += abs(lect(i,j)*lect(i,j));
        }
    }
    norm *= (tot_x/nx);
    norm *= (tot_y/ny);
    norm = sqrt(norm);

    cout << "L2 norm: " << norm << endl;
    fclose(f);

    return lect;
}


/** PARALLEL: when possible, allows AFMM to launch multiple processes in
    parallel.

    usage:
    parallel n

    parameters:
        n the number of maximum processes AFMM tries to create (unrelated to
            parallel implementations of LAPACK or BLAS, so check it out!).

*/
int commands::c_parallel(parsefile *obj, int argc,char *argv[])
{
    structure *p;
    int np;

    p = wantedParameters(argv[0], argc, 2, obj);

    if(sscanf(argv[1], "%20d", &np)!=1) {
         throw parsefile_commandError("parallel: error while reading "
             "the number of threads.\n");
     }

     if(np<1) {
         throw parsefile_commandError("parallel: how am I supposed to work without any process?\n");
     }

    p->number_of_threads_allowed=np;
    cout << "Multithreading allowed with "<<np<<" threads."<<endl;

    return 0;
}