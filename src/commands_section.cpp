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
#include "phys_constants.h"


/***************************************************************************
*     Section-level commands                                               *
****************************************************************************/


/*  SUBSTRATE: Define substrate relative refractive index index.
    Usage:
    substrate nr ni

    Parameters:
    nr: substrate refractive index, real part
    ni: substrate refractive index, imaginary part

    Notes:
    This command should be used once the size of the calculation window and
    the number of space harmonics to be taken has been defined.
    Executing this command blanks the structure!
*/
int commands::c_substrate(parsefile *obj, int argc,char *argv[])
{
    double nr;
    double ni;
    structure *p;

    p = wantedParameters(argv[0], argc, 3, obj);

    if(sscanf(argv[1], "%20lf", &nr)!=1 || sscanf(argv[2], "%20lf", &ni)!=1) {
        cerr<<"substrate: error while reading parameters.\n";
        return ERROR;
    }else{
        if((p=dynamic_cast<structure *>(obj))==NULL)
            throw parsefile_commandError("Incorrect dynamic cast: programming"
                " error:-(");
        p->cur->set_substrate(complex<double>(nr,ni));

    }

    return 0;
}

/*  RECTANGLE: Add a rectangular waveguide to the structure

    Usage:
    rectangle ng ngi wx wy px py

    Parameters:
    ngr: core refractive index (real part)
    ngi: core refractive index (imaginary part)
    wx: waveguide width
    wy: waveguide height
    px: x position of the center of the waveguide (referred to the center of
        the calculation window)
    py: y posizion of the center of the waveguide (referred to the center of
        the calculation window)

    Notes:
    This command should be used once the size of the calculation window and
    the number of space harmonics to be taken has been defined.
    The substrate index should have already been defined.

*/
int commands::c_rectangle(parsefile *obj, int argc,char *argv[])
{
    complex<double> n_g;
    double n_gi=0;
    double n_gr=1;
    double wx=0;
    double wy=0;
    double px=0;
    double py=0;
    structure *p = wantedParameters(argv[0], argc, 7, obj);

    if(sscanf(argv[1], "%20lf", &n_gr)!=1 ||
        sscanf(argv[2], "%20lf",&n_gi)!=1 ||
        sscanf(argv[3], "%20lf", &wx)!=1 ||
        sscanf(argv[4], "%20lf", &wy)!=1 || sscanf(argv[5], "%20lf", &px)!=1 ||
        sscanf(argv[6], "%20lf", &py)!=1) {
        throw parsefile_commandError("rectangle: error while reading"
            " parameter.\n");
    }else{
        n_g = complex<double>(n_gr, n_gi);
        p->cur->add_rectangle(n_g, wx, wy, px, py);
    }
    return 0;

}

/*  PML: Add perfectly matched layers to the waveguide structure.
         See "Use of grating theories in integrated optics"
         Silberstein et al., J. Opt. Soc. Am. A, Vol. 18, No. 11, November 2001

    Usage:
    pml alphar alphai wx wy

    Parameters:
    alphar: real part of the alpha parameter
    alphai: imaginary part of the alpha parameter


    Notes:
    This command should be used once the size of the calculation window and
    the number of space harmonics to be taken has been defined.
    The substrate index should have already been defined.
*/
int commands::c_pml(parsefile *obj, int argc,char *argv[])
{
    double alphar;
    double alphai;
    double wx;
    double wy;
    structure *p;

    p = wantedParameters(argv[0], argc, 5, obj);

    if(sscanf(argv[1], "%20lf", &alphar)!=1 ||
        sscanf(argv[2], "%20lf", &alphai)!=1  ||
        sscanf(argv[3], "%20lf", &wx)!=1 || sscanf(argv[4], "%20lf", &wy)!=1) {
            cerr<<"pml: error while reading parameter.\n";
            return ERROR;
    }else{

        p->cur->set_pml(wx, wy, complex<double>(alphar, alphai));
    }
    return 0;

}


/*  INPSTRUCT: write on a file the permettivity or permeability chart of the
        input structure (stands for INPut STRUCTure)

    Usage:
    inpstruct type size_x size_y file


    Parameters:
    type    must be {i|im|imo|ex|ey|ez|mux|muy|muz}, i indicates the refractive
            index (it will take er_x squared), while ex, ey or ez will indicate
            the x, y or z component of the permittivity tensor. On the same
            way, mux, muy or muz will indicate the x, y or z component of the
            permeability tensor. The im case is a little bit special, since
            it will give the refractive index calculated *without* representing
            the structure via the Fourier series. It is thus the "ideal case"
            in some way. imo is the same, but in the Optiwave format.

    size_x  number of points to be plot in the x direction

    size_y  number of points to be plot in the y direction

    file    filename which should be written

    File format:
        The file format used is extremely simple and is organized in order to
        be directly compatible with Gnuplot way of life.

*/
int commands::c_inpstruct(parsefile *obj, int argc,char *argv[])
{

    /*  Note May, 14 2013. The output (at least in Y) is not oriented correctly
        and the origin of coordinates is not correct?
        June 18, 2015: INDEED!!!
    */
    structure *p;

    p = wantedParameters(argv[0], argc, 5, obj);

    if(!p->cur->isSubstrateSet) {
        throw parsefile_commandError("inpstruct: The substrate should be set"
            " before trying to define structures.");
    }

    // p contains the current object and should therefore be a structure

    int sj =0;
    int si =0;

    sscanf(argv[2],"%20d", &sj);
    sscanf(argv[3],"%20d", &si);
    db_matrix out = p->cur->do_inpstruct(sj, si,argv[1]);
    
    double dx=p->tot_x / sj;
    double dy=p->tot_y / si;

    FILE *f=fopen(argv[4],"w");
    if(f==NULL)
        throw parsefile_commandError("inpstruct: I can not open the output"
            " file:-(");

    if(strcmp(argv[1],"imo")==0) {
        fprintf(f, "UPI3DRI 3.0\n");
        fprintf(f, "%d %d\n", out.getNcol(), out.getNrow());
        fprintf(f, "%lf %lf %lf %lf\n", 0.0, p->tot_x*1e6, 0.0, p->tot_y*1e6);
        for(int i=0; i<out.getNrow(); ++i) {
            for(int j=0; j<out.getNcol(); ++j) {
                fprintf(f, "%le, %le\n",
                    out(i,j).real(),
                    out(i,j).imag());
            }
            fprintf(f,"\n");
        }
    } else {
        fprintf(f, "# x         y        %s.real       %s.imag\n", argv[1],
            argv[1]);
        for(int i=0; i<out.getNrow(); ++i) {
            for(int j=0; j<out.getNcol(); ++j) {
                fprintf(f, "%le %le %le %le\n",
                    j*dx,
                    i*dy,
                    out(i,j).real(),
                    out(i,j).imag());
            }
            fprintf(f,"\n");
        }
    }
    fclose(f);
    cout << "Structure file written: "<< argv[4]<<"\n";

    return 0;
}

/** SELMODES: selects the modes using a defined rule and write on a file
        their shape. See OUTGMODES for a description of the file formats
        as well as the rules for file names.
        
    Usage:
    selmodes type rimco rule value1 value2 size_x size_y fileb

    Parameters:
    type    see OUTGMODES
    rimco   see OUTGMODES
    rule    must be {exco} (EXCitation Order).
            In this case the modes written by SELMODES are those carrying
            most of the power. The parameter value1 in this case can be
            {f|b} to indicate wether the forward power will have to be
            considered, or the backward power. The power will be calculated
            by considering the excitation coefficient, times the normalization
            coefficient (the same calculated by EIGENVC).
            value2 in this case specifies how many modes should be used.
            For example:
            
            selmodes Ex o exco f 3 101 101 modes
            
            will write on files modes_Ex_o_0.f3d, modes_Ex_o_1.f3d,
            modes_Ex_o_3.f3d the eigenfunction of the three modes which
            in the current section are carrying most of the power. Each file
            will be in the Optiwave f3d file format and will be composed by
            101x101 points. Of course, since the excitation coefficients need
            to be known, this has to be used after a propagation.
            NOTE: This rule must be used with a great care when there are modes
            which are not orthogonal in the problem being considered.
            
     value1 see rule
     value2 see rule
*/
int commands::c_selmodes(parsefile *obj, int argc,char *argv[])
{
    structure *p;

    field_e field;
    enum rimco_e_t rimc;


    if((p=dynamic_cast<structure *>(obj))==NULL)
        throw parsefile_commandError("Incorrect dynamic cast: programming"
            " error:-(");


    double f=CELERITY/p->lambda;
    double omega=2.0*M_PI*f;
    double k_0=omega/CELERITY;

    bool shiftToY = false;
    bool calcH = false;
    bool calcD = false;
    bool calcz = false;
    bool improve_representation=false;

    enum rule_t {exco} rule;

    if(argc==9){
        if(!p->cur->isSubstrateSet) {
            throw parsefile_commandError("selmodes: The substrate should be"
                " set before trying to write modal eigenfields.");
        }

        // p contains the current object and should therefore be a structure

        int sj =0;
        int si =0;

        int snux=p->dimx/2+1;
        int snuy=p->dimy/2+1;

        sscanf(argv[6],"%20d", &sj);
        sscanf(argv[7],"%20d", &si);
        if (si<snuy) si=snuy;
        if (sj<snux) sj=snux;

        cout << "Fields represented with: "<<sj<<" x "<<si<<" points.\n";

        // TODO: more or less the same code in propagation routines. Should
        // be kept together, if possible.
        if(strcmp(argv[1],"Ex")==0) {
            shiftToY = false;
            calcH = false;
            calcz = false;
        } else if(strcmp(argv[1], "Ey")==0) {
            shiftToY = true;
            calcH = false;
            calcz = false;
        } else if(strcmp(argv[1], "Ez")==0) {
            shiftToY = false;
            calcH = false;
            calcz = true;
        } else if(strcmp(argv[1],"Hx")==0) {
            shiftToY = false;
            calcH = true;
            calcz = false;
        } else if(strcmp(argv[1], "Hy")==0) {
            shiftToY = true;
            calcH = true;
            calcz = false;
        } else if(strcmp(argv[1], "Hz")==0) {
            shiftToY = false;
            calcH = true;
            calcz = true;
        } else if(strcmp(argv[1], "Dy")==0) {
            shiftToY = true;
            calcH = false;
            calcD = true;
            calcz = false;
        } else if(strcmp(argv[1], "Dx")==0) {
            shiftToY = false;
            calcH = false;
            calcD = true;
            calcz = false;
        }else if(strcmp(argv[1], "Dz")==0) {
            shiftToY = false;
            calcH = false;
            calcD = true;
            calcz = true;
        }else if(strcmp(argv[1], "Ex2")==0) {
            shiftToY = false;
            calcH = false;
            calcD = false;
            calcz = false;
            improve_representation = true;
        }else if(strcmp(argv[1], "Ey2")==0) {
            shiftToY = true;
            calcH = false;
            calcD = false;
            calcz = false;
            improve_representation = true;
        }else if(strcmp(argv[1], "Ez2")==0) {
            shiftToY = false;
            calcH = false;
            calcD = false;
            calcz = true;
            improve_representation = true;
        }else {
            throw parsefile_commandError("selmodes: unrecognized field type."
                " It should be {Ex|Ey|Ez|Ex2|Ey2|Ez2|Hx|Hy|Hz|Dx|Dy|Dz}.");
        }

        if(strcmp(argv[2],"r")==0) {
            rimc=R;
        } else if(strcmp(argv[2],"i")==0) {
            rimc=I;
        } else if(strcmp(argv[2],"m")==0) {
            rimc=M;
        } else if(strcmp(argv[2],"c")==0) {
            rimc=C;
        } else if(strcmp(argv[2],"o")==0) {
            rimc=O;
        } else {
            throw parsefile_commandError("selmodes: unrecognized rimco. It "
                "should be {r|i|m|c|o}.");
        }

        if(strcmp(argv[3], "exco")==0) {
            rule=exco;
        } else {
            throw parsefile_commandError("selmodes: unrecognized rule. It "
                " should be {exco}.");
        }

        // Check if the C matrix is empty
        if (p->cur->W.isEmpty()) {
            throw parsefile_commandError(
                "selmodes: the problem must be set and modes calculated with"
                " solve before trying to describe them on a file.");
        } else if(p->cur->B.isEmpty()) {
            throw parsefile_commandError(
                "selmodes: rule=exco needs to have the excitation coefficients"
                " calculated with a propagation.");
        }
        double dx=p->tot_x / sj;
        double dy=p->tot_y / si;

        db_matrix out(si, sj);

        int mi=0;
        int mm, ll;



        // compute epsilon and mu constants if needed
        db_matrix epsilonxy;
        db_matrix epsz;
        db_matrix muz;

        getEpsilonAndMu(*p->cur, calcD, calcz,
            calcH,improve_representation, epsilonxy, epsz, muz);


        // The V matrix will be used if we want to compute Ez,Ez2,Dz, Hx and Hy
        // if the V matrix does not exist, it is created
        if ((calcH && !calcz) || (calcz && !calcH)) {
            if (p->HxField && p->HyField) {
                if (p->cur->V.getNcol() == 0 && p->cur->V.getNrow() == 0) {
                    p->cur->V = p->cur->assembleVmatrix();
                }
            } else {
                cerr<<"selmodes: You need to use 'wants both' "
                        " to output the required mode" << endl;
                return ERROR;
            }
        }

        // Matrix containing the modal power of each mode.
        double powermodes[p->cur->B.getNrow()];

        double integral;
        int nmodes=0;

        if(sscanf(argv[5], "%d", &nmodes)!=1) {
            throw parsefile_commandError(
                "selmodes: could not read the number of modes to be written.");
        } else {
            if (nmodes<0 || nmodes > p->cur->B.getNrow()) {
                throw parsefile_commandError(
                "selmodes: Invalid value of number of modes.");
            }
        }

        // Make sort that e points to the right excitation vector.
        
        db_matrix *ee;
        cout<<"Sorting modes on the basis of the power which "
            " is coupled on them in the ";
        if(strcmp(argv[4],"f")==0) {
            cout<<"forward ";
            ee = &(p->cur->sWp);
        } else if(strcmp(argv[4],"b")==0) {
            cout<<"backward ";
            ee = &(p->cur->sWm);
        }
        cout << "direction."<<endl;

        bool writemodes[p->cur->B.getNrow()];
        complex<double>e;
        complex<double>g;
        double order;
        double quality;
        double r=0;

        // Calculate the power carried by each mode.
        for (int jj=0; jj<p->cur->B.getNrow();++jj) {
            writemodes[jj]=false;
            e=ee->operator()(jj,0);
            integral = abs(section::integralPoynting(p, p->cur->W,p->cur->V, jj)
                *p->tot_x*p->tot_y);

            powermodes[jj]=abs(e)*abs(e)*integral;
        }

        // Now that the power is known, we can search for the maximum and
        // flag it as a mode to be written. We continue in such a way
        // nmodes times, each time searching for the maximum avoid what
        // it has already been flagged.

        int modesorder[nmodes];

        double max=0;
        int maxindex=-1;

        // This sorting method would be pretty inefficient for big values
        // of nmodes. However, since the calculations are not very
        // heavy, and nmodes will be in most cases very low (1-5), this should
        // be ok.
        for (int k=0; k<nmodes; ++k) {
            max=0;
            maxindex=-1;
            for (int jj=0; jj<p->cur->B.getNrow();++jj) {
                if(powermodes[jj]>max && !writemodes[jj]) {
                    maxindex=jj;
                    max=powermodes[jj];
                }
            }
            if(maxindex>=0) {
                modesorder[k]=maxindex;
                writemodes[maxindex]=true;
            }
            cout << "mode "<<k<<", power "<<max<<endl;
        }


        // Write the selected modes in the right order.
        for (int k=0; k<nmodes; ++k) {
            // Pick up the modes in the right order.
            int jj=modesorder[k];
            db_matrix select(p->cur->W.getNrow(), 1);
            e=p->cur->B(jj,jj);
            g = p->getEffectiveIndex(e);
            order = p->getOrder(e, p->cur->radius);
            quality = p->getQuality(e);



            // Select the jj-order mode
            // so that W times jj selects the mode
            select(jj, 0)=1;

            // We compute the field:
            db_matrix unused;

            out = getField(*p->cur, select,
                    shiftToY, calcH, snux, snuy, sj, si,
                    calcD, epsilonxy, epsz, muz,calcz,improve_representation,
                    false, unused);
            size_t buf_size=strlen(argv[8])+strlen(argv[1])+strlen(argv[2])+40;
            char *s=new char[buf_size];

            // We can then write the results on the output file
            if (rimc == O) {
                // Optiwave format
                snprintf(s, buf_size, "%s_%s_%s_%d.f3d", argv[8], argv[1],
                        argv[2], mi++);
                FILE *f=fopen(s,"w");
                p->fileoutputOW(out, f, dx*1e6, dy*1e6);
                cout << "Mode file written (Optiwave format): "<< s<<"\n";
            } else {
                // Gnuplot format. Put a small header.
                snprintf(s, buf_size, "%s_%s_%s_%d.mode", argv[8], argv[1],
                     argv[2], mi++);
                FILE *f=fopen(s,"w");

                if(p->cur->useAzimuthalOrder) {
                    fprintf(f, "# order=%lf, Q= %lf, power=%e\n",order,
                        quality, powermodes[jj]);
                } else {
                    fprintf(f, "# neff=(%lf, %e), power=%e\n",g.real(),
                        g.imag(), powermodes[jj]);
                }
                p->fileoutputGP(out, f, rimc, dx, dy);
                cout << "Mode file written (Gnuplot format): "<< s<<"\n";
            }
            if(s) delete[] s;

        }
    } else {
        cerr<<"selmodes: invalid number of parameters.";
        return ERROR;
    }
    return 0;
}



/** OUTGMODES: write on a file the calculated field associated to all
        interesting modes. The user must specify the output field
        characteristics and the program will automatically append modal
        informations to the file name provided.

    Usage:
    outgmodes type rimco size_x size_y fileb

    Parameters:
    type    must be {Ex|Ey|Ez|Ex2|Ey2|Ez2|Hx|Hy|Hz|Dx|Dy|Dz} and indicates
            which transverse field component
            must be drawn. The number 2 refers as the improved representation
            of the fields

    rimco   must be {r|i|m|c|o} and indicates if we want to output the real
            part (r), the imaginary part (i), the magnitude (m), the complex
            magnitude (c) of the fields or (o) in the Optiwave format.

    size_x  number of points to be plot in the x direction

    size_y  number of points to be plot in the y direction

    fileb   filename base which should be written. The file name will be
            written in the form fileb_type_rimc_nn.mode
            where nn is the progressive mode number and, fileb, type and rimc
            are the calling parameters of outgmodes

    File format:
        If the output format is r|i|m|c, the file format used is extremely
        simple and is organized in order to be directly compatible with
        Gnuplot.
        The extension used for the file will be mode.
        Have a look at the header:
        % neff = (3.233212, 2.1233242e-8)
        % x         y        Ex.real       Ex.imag
        0.000000e+00 0.000000e+00 -2.881303e-04 6.743971e-05
        1.500000e-08 0.000000e+00 -2.261833e-04 5.277011e-05
        3.000000e-08 0.000000e+00 -7.094755e-05 1.610860e-05
        ... and so on

        If the "o" format is used, the Optiwave file format is instead used,
        which can be useful since some commands of AFMM understand this format.
        The extension used for the file will be f3d.
        Here is a typical file header:
        BCF3DCX 3.0
        128 128
        0 30.000000 0 30.000000
        1.440000,0
        1.440000,0
        1.440000,0
        ... and so on

        The first line contains a fixed header. This is for example plainly
        ignored by AFMM when performs a read operation on this file format.
        The second line gives the number of points in x and y.
        The third line gives first the interval in x direction (in microns)
        and in the y direction (in microns). In this case, we will have a grid
        of 128 x 128 points which will cover the intervals between 0 and 30
        microns in x and in y.
        Then follows the field or refractive index data, line by line.
*/

int commands::c_outgmodes(parsefile *obj, int argc,char *argv[])
{
    structure *p;

    if((p=dynamic_cast<structure *>(obj))==NULL)
        throw parsefile_commandError("Incorrect dynamic cast: programming"
            " error:-(");

    if(argc==6) {
        int sj=0;
        int si=0;
        enum rimco_e_t rimc;
        sscanf(argv[3],"%20d", &sj);
        sscanf(argv[4],"%20d", &si);
        if(strcmp(argv[2],"r")==0) {
            rimc=R;
        } else if(strcmp(argv[2],"i")==0) {
            rimc=I;
        } else if(strcmp(argv[2],"m")==0) {
            rimc=M;
        } else if(strcmp(argv[2],"c")==0) {
            rimc=C;
        } else if(strcmp(argv[2],"o")==0) {
            rimc=O;
        } else {
            throw parsefile_commandError("outgmodes: unrecognized rimc. Should"
                " be {r|i|m|c|o}.");
        }
        p->cur->do_outgmodes(sj, si, argv[1], rimc, argv[5]);
    } else {
        cerr<<"outgmodes: invalid number of parameters.";
        return ERROR;
    }
    return 0;
}


/*  LOWINDEX: Define the lowest value allowed for the effective index to be
    considered a guided mode. If this command is not used, this value is assumed
    to be the substrate's refractive index. This command can be useful when
    dealing with asymmetric structures, in order to avoid considering modes
    which are not guided.


    Usage:
    lowindex lr li

    Parameters:
    lr lowest value of the real part of the refractive index
    li lowest value of the imaginary part of the refractive index

*/
int commands::c_lowindex(parsefile *obj, int argc,char *argv[])
{
    double lr;
    double li;
    structure *p;
    if(argc==3){
        if(sscanf(argv[1], "%20lf", &lr)!=1 ||sscanf(argv[2], "%20lf", &li)!=1){
            cerr<<"lowindex: error while reading parameter.\n";
            return ERROR;
        }else{
            if((p=dynamic_cast<structure *>(obj))==NULL)
                throw parsefile_commandError("Incorrect dynamic cast:"
                    " programming error:-(");

            // p contains the current object: it should be a structure

            p->cur->set_lowindex(complex<double>(lr,li));
        }
    } else {
        cerr<<"lowindex: invalid number of parameters.\n";
        return ERROR;
    }
    return 0;

}


/*  HIGHINDEX: Define the highest value allowed for the effective index to be
    considered an interesting mode. If this command is not used, this value is
    assumed to be the highest refractive index used for waveguide definition.
    This command can be useful when the automatic calculations are not working,
    for example because of several waveguide sections overlapping.

    Usage:
    highindex hr hi

    Parameters:
    hr highest value of the real part of the refractive index
    hi highest value of the imaginary part of the refractive index

*/
int commands::c_highindex(parsefile *obj, int argc,char *argv[])
{
    double lr;
    double li;
    structure *p;
    if(argc==3){
        if(sscanf(argv[1], "%20lf", &lr)!=1 || sscanf(argv[2],
            "%20lf", &li)!=1){
            cerr<<"highindex: error while reading parameter.\n";
            return ERROR;
        }else{
            if((p=dynamic_cast<structure *>(obj))==NULL)
                throw parsefile_commandError("Incorrect dynamic cast:"
                    " programming error:-(");


            // p contains the current object and it should be a structure
            p->cur->set_highindex(complex<double>(lr,li));
        }
    } else {
        cerr<<"highindex: invalid number of parameters.\n";
        return ERROR;
    }
    return 0;

}

/*  PML_TRANSF Use PMLs with coordinate transform as described in [1]

    Usage:
    pml_transf qx qy gamma

    Parameters:
    qx: size of the transform region in x (should be <= to the x size of the
        calculation window.
    qy: size of the transform region in y (should be <= to the y size of the
        calculation window.
    gamma: complex factor as described in [2]
*/
int commands::c_pml_transf(parsefile *obj, int argc,char *argv[])
{
    double qdx;
    double qdy;
    double qgammar;
    double qgammai;
    structure *p;
    if(argc==5){
        if(sscanf(argv[1], "%20lf", &qdx)!=1 || sscanf(argv[2],
            "%20lf", &qdy)!=1
           || sscanf(argv[3], "%20lf", &qgammar)!=1 ||
           sscanf(argv[4], "%20lf", &qgammai)!=1){
            cerr<<"pml_transf: error while reading parameters.\n";
            return ERROR;
        }else{
            if((p=dynamic_cast<structure *>(obj))==NULL)
                throw parsefile_commandError("Incorrect dynamic cast:"
                    " programming error:-(");

            p->cur->set_pml_transf(qdx,qdy,complex<double>(qgammar, qgammai));
        }
    } else {
        cerr<<"PML_TRANSF: invalid number of parameters.\n";
        return ERROR;
    }
    return 0;

}

/*  INDFILE: read refractive index distribution from an input file. This
    command overwrites the previous settings regarding the size of the
    calculation window (command SIZE), as well as the substrate settings
    (command SUBSTRATE). For this reason, this command can be used just after
    the HARMONICS command, since the definition of the calculation window is
    read from the input file. This command does NOT define the minimum and
    maximum value of the indices to be used for the guided mode search. The
    user should define them by using the LOWINDEX and HIGHINDEX commands.

    Usage:
    indfile mult filename

    Parameters:
    mult: a multiplier to be considered for length calculations. If the sizes
        in the input file are specified in micrometers, mult should be 1e-6.
        If they are in meters, it should be 1 and so on.
    filename: the name of the file to be used (it should not contain spaces).

    Notes:
    The file format to be used will be as follows. The first line is ignored.
    This can be an header or a description of the file format.
    The second line contains the number of x and y subdivisions, separated by
    a space.
    The third line specifies the minimum and maximum ranges of the calculation
    window, for the x and y axis. Remember that the "mult" parameter should be
    set accordingly to the convenction used in your file.
    The input file specifies a sampled index distribution, point by point, line
    by line. The refractive index is specified with its real part and imaginary
    part separated by a comma.

*/
int commands::c_indfile(parsefile *obj, int argc,char *argv[])
{
    complex<double> n_g;
    structure *p;
    double mult;

    int nx;
    int ny;

    double tot_x;
    double tot_y;

    int i,j;

    if(argc==3){
        if(sscanf(argv[1], "%20lf", &mult)!=1){
            cerr<<"indfile: error while reading parameter.\n";
            return ERROR;
        }else{
            if((p=dynamic_cast<structure *>(obj))==NULL)
                throw parsefile_commandError("Incorrect dynamic cast:"
                    " programming error:-(");

            cout << "Opening index file: " <<argv[2]<<"\n";
            cout.flush();
            try {
            db_matrix rd = read_file(argv[2], tot_x, tot_y, mult);
            // Set the size of the calculation window as set by the file
            // settings.

            p->tot_x=tot_x;
            p->tot_y=tot_y;

            p->cur->store_refractive_index(rd);
            cout << "Just read a file index distribution: "<< argv[2] << "\n";
            } catch (db_matrix_error E) {
                throw parsefile_commandError("indfile: Could not allocate"
                    " memory for index structure. Is there an error in the"
                    " file? Normally modern computer's memories should be"
                    " enough for most practical cases...");
            }
        }
    } else {
        cerr<<"indfile: invalid number of parameters.\n";
        return ERROR;
    }
    return 0;
}


/*  BEND: Define the bending ratio.

    Usage:
    bend r

    Parameters:
    r: curvature radius (around y axis) of the CENTER of the calculation region
    This value is used also when calculating the effective index of the modes.
    For this reason, the user should choose where to place the waveguide
    in the calculation window. Some use the internal radius of the waveguide as
    the curvature radius, others the central radius, other the external radius.
    Be aware of the choice made!
*/
int commands::c_bend(parsefile *obj, int argc,char *argv[])
{
    double l;
    structure *p;
    if(argc==2){
        if(sscanf(argv[1], "%20lf", &l)!=1){
            cerr<<"bend: error while reading parameter.\n";
            return ERROR;
        }else{
            if((p=dynamic_cast<structure *>(obj))==NULL)
                throw parsefile_commandError("Incorrect dynamic cast:"
                    " programming error:-(");

            // p contains the current object and should be a structure

            p->cur->set_bend(l);
        }
    } else {
        cerr<<"bend: invalid number of parameters.\n";
        return ERROR;
    }
    return 0;
}

/*  SPECTRUM: writes on a file the real and imaginary part of ALL effective
        indices of ALL modes calculated.

    Usage:
    spectrum filename

    Parameters:
    filename the name of the file on which to write
*/
int commands::c_spectrum(parsefile *obj, int argc,char *argv[])
{
    structure *p;
    if(argc==2){
        if((p=dynamic_cast<structure *>(obj))==NULL)
            throw parsefile_commandError("Incorrect dynamic cast: "
                "programming error:-(");

        // Check if the C matrix is empty
        if (p->cur->W.isEmpty()) {
            throw parsefile_commandError(
                "spectrum: the problem must be set and modes calculated with"
                "solve before trying to write the spectrum on a file.");
        }

        FILE *f=fopen(argv[1],"w");
        if(f==NULL)
            throw parsefile_commandError("spectrum: Can not open output"
                " file:-(");

        fprintf(f, "# real       imag\n");
        for(int i=0; i<p->cur->B.getNrow(); ++i) {
            complex<double>e=p->cur->B(i,i);
            e = p->getEffectiveIndex(e);

            fprintf(f, "%le %le\n", e.real(), e.imag());
        }
        fclose(f);
        cout << "Spectrum file written: "<< argv[1]<<" ("<<
            p->cur->B.getNrow()<<" entries)\n";

    } else {
        cerr<<"spectrum: invalid number of parameters.\n";
        return ERROR;
    }
    return 0;

}

/*  SECTION: Define a section of the given length

    Usage:
    section tl

    Parameters:
    tl total length in meters
*/
int commands::c_section(parsefile *obj, int argc,char *argv[])
{
    double tl;
    structure *p;
    cout<<endl;
    if(argc==2){
        if(sscanf(argv[1], "%20lf", &tl)!=1){
            cerr<<"section: error while reading parameters.\n";
            return ERROR;
        }else{
            if((p=dynamic_cast<structure *>(obj))==NULL)
                throw parsefile_commandError("Incorrect dynamic cast:"
                    " programming error:-(");

            p->do_section(tl);
        }
    } else {
        cerr<<"section: invalid number of parameters.\n";
        return ERROR;
    }
    return 0;
}

/*  ORDER: specify that the propagation constants will be selected in terms of
    azimuthal order, instead of of effective indices. This can be particularly
    useful when the curvature radius should not be used for any reason (for
    example, it is set at zero).

    Usage:
    order min max

    Parameters:
    min the minimum azimuthal order retained
    max the maximum azimuthal order retained
*/
int commands::c_order(parsefile *obj, int argc,char *argv[])
{
    double min;
    double max;

    structure *p;

    if(argc==3){
        if(sscanf(argv[1], "%20lf", &min)!=1 || sscanf(argv[2],
            "%20lf", &max)!=1)
        {
            cerr<<"order: error while reading parameter.\n";
            return ERROR;
        } else {
            if((p=dynamic_cast<structure *>(obj))==NULL)
                throw parsefile_commandError("Incorrect dynamic cast:"
                    " programming error:-(");

            p->cur->do_order(min, max);
        }
    } else {
        cerr<<"order: invalid number of parameters.\n";
        return ERROR;
    }
    return 0;
}

/** COEFFICIENT determine the transmission coefficient of a mode in the
    structure, specified by the real part of the effective index.

    Usage:
    coefficient fb a

    Parameters:
    fb  can be {f|b} where f stands for forward and b for backward
    a   the real part of the refractive index of the mode to analyze

    or:
    coefficient fb "id" pos

    Parameters:
    fb  can be {f|b} where f stands for forward and b for backward
    id  indicates that the selection must be done on the mode position index
        (see the modepos command).
    pos the position index of the wanted mode.

    Return value:
    ans will be an array of three elements:
    - ans[0] is the real part of the coefficient
    - ans[1] is the imaginary part of the coefficient
    - ans[2] is the squared module of the coefficient

*/
int commands::c_coefficient(parsefile *obj, int argc,char *argv[])
{
    double value;
    int iMode;      // modal position index resulting from the selection
    structure *p;
    complex<double> e;


    if(strcmp(argv[2], "id")==0) {
        p = wantedParameters(argv[0], argc, 4, obj);
        sscanf(argv[3], "%20d", &iMode);
    } else {
        p = wantedParameters(argv[0], argc, 3, obj);
        if(sscanf(argv[2], "%20lf", &value)!=1) {
            throw parsefile_commandError("coefficient: error while reading the"
                " refractive index.\n");
        }
        iMode = p->cur->select_real(value);
    }

    if (iMode<0) {
        cout<<"I could not find the given mode.\n";
        return 0;

    }
    cout << "Selected mode: ";
    cout << p->getEffectiveIndex(p->cur->B(iMode,iMode)) << "\n";

    if(strcmp(argv[1],"f")==0) {
        if(p->cur->sWp.isEmpty())
            throw parsefile_commandError("coefficient: excitation coefficients"
                " have not been calculated yet.\n");
        else {
            cout<<"Forward ";
            e = p->cur->sWp(iMode, 0);
        }
    } else if(strcmp(argv[1],"b")==0) {
        if(p->cur->sWm.isEmpty())
            throw parsefile_commandError("coefficient: excitation coefficients"
                " have not been calculated yet.\n");
        else {
            cout<<"Backward ";
            e = p->cur->sWm(iMode, 0);
        }
    } else {
        throw parsefile_commandError("coefficient: unrecognized direction "
            "of the excitation. It must be {b|f}.\n");
    }
    cout << "mode coefficient: " << e << "\n";
    cout << "Squared module: " << abs(e)*abs(e) << "\n";

    // Save the calculated data as an array in the output variable "ans".
    double *s =new double[3];
    s[0]=e.real();
    s[1]=e.imag();
    s[2]=abs(e)*abs(e);

    p->nP->insertArray("ans", 3, s);
    delete[] s;
    return 0;
}

/** EIGENVC writes on a file all the excitation coefficients of the modes
    propagating in the current section. This command should be used after
    propagation.

    Usage:
    eigenvc fba name

    Parameters:
    fba can be {f|b|a} where f stands for forward and b for backward
        If this parameter is a, a more complex set of
        information is given for each line:
        - real and imaginary part of the effective index,
        - real and imaginary part of the forward coefficients
        - real and imaginary part of the backward coefficients
        - the integral of the Poynting vector in the calculation window
    name the filename.

*/
int commands::c_eigenvc(parsefile *obj, int argc,char *argv[])
{
    double value;
    structure *p;
    db_matrix *e;
    bool all=false;
    double integral;

    p = wantedParameters(argv[0], argc, 3, obj);

    if(strcmp(argv[1],"f")==0) {
        cout<<"Forward ";
        e = &(p->cur->sWp);
    } else if(strcmp(argv[1],"b")==0) {
        cout<<"Backward ";
        e = &(p->cur->sWm);
    } else if(strcmp(argv[1],"a")==0) {
        cout<<"All ";
        all = true;
    } else {
        throw parsefile_commandError("eigenvc: unrecognized direction "
            "of the excitation. It must be {b|f|a}.\n");
    }
    // At first, we check if the needed matrices are existing or not
    if (p->cur->sWp.getNcol()==0 && p->cur->sWm.getNcol()==0) {
        // We close the output file, since that will not be done processing
        // the exception.
        throw parsefile_commandError("eigenvc: this command should be"
            " used after a field propagation.\n");
    }

    FILE *f=fopen(argv[2],"w");
    if(f==NULL)
        throw parsefile_commandError("eigenvc: Can not open output file:-(");

    if(!all) {
        fprintf(f, "# real       imag\n");
    } else {

        fprintf(f,"# n_real\tn_imag\tf_real\tf_imag\tb_real\tb_imag\tint_S\n");
    }
    for(int i=0; i<p->cur->B.getNrow(); ++i) {
        if (!all) {
            complex<double>s=(*e)(i,0);
            fprintf(f, "%le %le\n", s.real(), s.imag());
        } else {
            complex<double>e=p->cur->B(i,i);
            e = p->getEffectiveIndex(e);

            integral = section::integralPoynting(p, p->cur->W,p->cur->V, i)
                *p->tot_x*p->tot_y;

            fprintf(f, "%le %le %le %le %le %le %le\n",
                e.real(), e.imag(),
                p->cur->sWp(i,0).real(),p->cur->sWp(i,0).imag(),
                p->cur->sWm(i,0).real(),p->cur->sWm(i,0).imag(),
                integral);
        }
    }
    fclose(f);
    cout << "Eigenvc file written: "<< argv[2]<<" ("<<p->cur->B.getNrow()
        <<" entries)\n";

    return 0;
}

/** EIGENAM writes on a file all the Fourier coefficients of the eigenvalues
        (matrices W or V).

    Usage:
    eigenam eh name

    Parameters:
    eh can be {e|h} specifies the electric (e) or the magnetic (h) field
    name the filename.

*/
int commands::c_eigenam(parsefile *obj, int argc,char *argv[])
{
    double value;
    structure *p;
    db_matrix *e;

    p = wantedParameters(argv[0], argc, 3, obj);
    cout << "Writing all the eigenvalues of the ";
    if(strcmp(argv[1],"e")==0) {
        cout<<"electric field\n";
        e = &(p->cur->W);
    } else if(strcmp(argv[1],"h")==0) {
        cout<<"magnetic field\n";
        e = &(p->cur->V);
    } else {
        throw parsefile_commandError("eigenam: unrecognized field. It must be"
            " {e|h}.\n");
    }

    FILE *f=fopen(argv[2],"w");
    if(f==NULL)
        throw parsefile_commandError("eigenam: Can not open output file:-(");

    for(int i=0; i<p->cur->W.getNrow(); ++i) {
        for(int j=0; j<p->cur->W.getNcol(); ++j) {
            fprintf(f, "(%le, %le) ", (*e)(i,j).real(), (*e)(i,j).imag());
        }
        fprintf(f, "\n");
    }
    fclose(f);
    cout << "Eigenam file written: "<< argv[2]<<" ("<<p->cur->W.getNrow()<<
        " x " <<p->cur->W.getNcol()<<" entries)\n";

    return 0;

}


/** NORM: calculate the L2 norm of the vector of excitation coefficients for
    the current section

    Usage:
    norm fb

    Parameters:
    fb  can be {f|b} where f stands for forward and b for backward

*/
int commands::c_norm(parsefile *obj, int argc,char *argv[])
{
    structure *p;
    double norm=0;

    p = wantedParameters(argv[0], argc, 2, obj);
    db_matrix A;
    if(strcmp(argv[1],"f")==0) {
        cout<<"Norm calculation forward ";
        A=p->cur->W*p->cur->sWp;

    } else if(strcmp(argv[1],"b")==0) {
        cout<<"Norm calculation backward ";
        A=p->cur->W*p->cur->sWm;
    } else {
        throw parsefile_commandError("norm: unrecognized direction "
            "of the excitation. It must be {b|f}.\n");
    }

    for(int i=0; i<p->cur->sWp.getNrow(); ++i) {
        norm += abs(A(i,0))*abs(A(i,0));
    }
    norm = sqrt(norm/p->cur->sWp.getNrow());

    cout << "L2 norm: " << norm << "\n";

    p->insertVar("ans",norm);

    return 0;
}


/** SELECT allow the user specify the current section.

    Usage:
    select a

    Parameters:
    a   the number of the wanted section, numbered starting from 1.

*/
int commands::c_select(parsefile *obj, int argc,char *argv[])
{
    unsigned int number;
    structure *p;
    p = wantedParameters(argv[0], argc, 2, obj);

    if(sscanf(argv[1], "%20d", &number)!=1) {
        throw parsefile_commandError("select: error while reading "
            "the refractive index.\n");
    }

    p->do_select(number);

    return 0;

}


/** MATDEV Determine the matrix development strategy

    Usage:
    matdev type

    Parameters:
    type    can be {af|an|be|nd|la|nf|nfs|las} where

    af means multiplication after construction strategy,

    an indicates still a multiplication after construction as in v. 1.2.3

    be indicates multiplication before construction as in v. 1.3.

    nd is as litterature

    la is Lalanne type developments. The alpha parameter is thus required:
        matdev la alpha

    nf indicates the normal vector field to be loaded:
        matdev nf mult file_x.f3d file_y.f3d

developpement for symmetry
    las is Lalanne type developments. The alpha parameter is thus required:
        matdev la alpha

    nfs indicates the normal vector field to be loaded:
        matdev nfs mult file_x.f3d file_y.f3d


    The default developments type which is used is af. This command should be
    used before SOLVE.

*/
int commands::c_matdev(parsefile *obj, int argc,char *argv[])
{
    structure *p;
    if(argc<=1) {
        cerr<<"pml_type";
        throw parsefile_commandError(": invalid number of parameters.");
    }
    if((p=dynamic_cast<structure *>(obj))==NULL)
        throw parsefile_commandError("Incorrect dynamic cast: programming"
            " error:-(");

    cout << "Matrix developments: ";
    if(strcmp(argv[1],"af")==0) {
        cout<<"af; multiplication after construction.\n";
        p->cur->crtr = &p->cur->PMLafterOPT;
    } else if(strcmp(argv[1],"an")==0) {
        cout<<"an; multiplication after construction (non optimized, like "
        "in v. 1.2.3).\n";
        p->cur->crtr = &p->cur->PMLafter;
    } else if(strcmp(argv[1],"be")==0) {
        cout<<"be; multiplication before construction.\n";
        p->cur->crtr = &p->cur->PMLbefore;
    } else if(strcmp(argv[1],"nd")==0) {
        cout<<"nd; compact matrix development.\n";
        p->cur->NonDev.setAlpha(-1.0);
        p->cur->crtr = &p->cur->NonDev;
    } else if(strcmp(argv[1],"la")==0) {
        if(argc<3) {
            throw parsefile_commandError("matdev: specify the alpha parameter"
                " for the Lalanne matrix developments.\n");
        }
        double alpha;

        if(sscanf(argv[2], "%20lf", &alpha)!=1){
            cerr<<"matdev: error while reading parameter.\n";
            return ERROR;
        }
        cout<<"la; Lalanne-type with alpha = "<<alpha<<".\n";

        p->cur->NonDev.setAlpha(alpha);
        p->cur->crtr = &p->cur->NonDev;

    } else if(strcmp(argv[1],"nf")==0) {

        cout<<"nf: normal field read from "<<argv[3]<<" and "<<argv[4]<<endl;
        double mult;
        if(argc<5) {
            throw parsefile_commandError("matdev: specify mult, file_x.f3d and"
                " file_y.f3d\n");
        }
        if(sscanf(argv[2], "%20lf", &mult)!=1){
            cerr<<"matdev: error while reading parameter mult.\n";
            return ERROR;
        }

        p->cur->crtr = &p->cur->NormalField;
        db_matrix file_x = read_file(argv[3], p->tot_x, p->tot_y, mult);

        file_x *= 1.0/file_x.getNcol()/file_x.getNrow();

        db_matrix file_x_fft = file_x.ifft2();

        p->cur->NormalField.Nr_fft = fft2cent(p->dimx, p->dimy, file_x_fft);

        db_matrix file_y = read_file(argv[4], p->tot_x, p->tot_y, mult);
        file_y *= 1.0/file_y.getNcol()/file_y.getNrow();
        db_matrix file_y_fft = file_y.ifft2();

        p->cur->NormalField.Nz_fft = fft2cent(p->dimx, p->dimy, file_y_fft);

    } else if(strcmp(argv[1],"las")==0) {
        if(argc<3) {
            throw parsefile_commandError("matdev: specify the alpha parameter"
                " for the Lalanne (symmetric) matrix developments.\n");
        }
        double alpha;

        if(sscanf(argv[2], "%20lf", &alpha)!=1){
            cerr<<"matdev: error while reading parameter.\n";
            return ERROR;
        }
        cout<<"la: (symmetric) Lalanne-type with alpha = "<<alpha<<".\n";

        p->cur->NonDevSym.setAlpha(alpha);
        p->cur->crtr = &p->cur->NonDevSym;

    } else if(strcmp(argv[1],"nfs")==0) {

        cout<<"nfs: normal field (symmetric) read from "<<argv[3]<<
            " and "<<argv[4]<<endl;
        double mult;
        if(argc<5) {
            throw parsefile_commandError("matdev: specify mult, file_x.f3d and"
                " file_y.f3d\n");
        }
        if(sscanf(argv[2], "%20lf", &mult)!=1){
            cerr<<"matdev: error while reading parameter mult.\n";
            return ERROR;
        }

        p->cur->crtr = &p->cur->NormalFieldSym;
        db_matrix file_x = read_file(argv[3], p->tot_x, p->tot_y, mult);

        db_matrix file_y = read_file(argv[4], p->tot_x, p->tot_y, mult);

        db_matrix Nxx  = file_x;
        db_matrix Nyy = file_y;
        db_matrix Nxy = file_x;

        Nxx.hadamard(file_x);
        Nyy.hadamard(file_y);
        Nxy.hadamard(file_y);

        // Save the matrices:
        p->cur->NormalFieldSym.Nxx_input = Nxx;
        p->cur->NormalFieldSym.Nxy_input = Nxy;
        p->cur->NormalFieldSym.Nyy_input = Nyy;

        // normalise the matrices for the fft:

        Nxx *= 1.0/file_x.getNcol()/file_x.getNrow();
        Nyy *= 1.0/file_x.getNcol()/file_x.getNrow();
        Nxy *= 1.0/file_x.getNcol()/file_x.getNrow();

        db_matrix Nxx_fft = Nxx.fft2();
        db_matrix Nxy_fft = Nxy.fft2();
        db_matrix Nyy_fft = Nyy.fft2();

        p->cur->NormalFieldSym.Nxx_fft = fft2cent(p->dimx, p->dimy, Nxx_fft);
        p->cur->NormalFieldSym.Nyy_fft = fft2cent(p->dimx, p->dimy, Nyy_fft);
        p->cur->NormalFieldSym.Nxy_fft = fft2cent(p->dimx, p->dimy, Nxy_fft);

    } else {
        throw parsefile_commandError("matdev: unrecognized development. "
            "It must be {af|an|be|la|nf|nfs|las}.\n");
    }

    if(!p->symmetryNotTakenIntoAccount() &&
        p->cur->crtr != &p->cur->NormalFieldSym &&
        p->cur->crtr != &p->cur->NonDevSym )
    {
        throw parsefile_commandError("matdev: when symmetry is taken into"
            " account, development type must be {nfs|las}.\n");
    }
    return 0;
}

/** POWER: calculate the power flowing in a given section
    This commands computes the power carried by all the modes propagating in
    a given section in a given direction (forward or backward).
    A mode is considered as propagating when the real part of its effective
    refractive index exceeds a certain threshold.

    Usage:
    power fb threshold

    Parameters:
        fb can be {f|b} if the direction is forward or backward.
        threshold is the threshold for the real part of the refractive index
            for a propagation mode to be considered as propagating.
*/
int commands::c_power(parsefile *obj, int argc,char *argv[])
{
    structure *p;
    db_matrix *e;
    double power=0;
    double threshold=0;

    p = wantedParameters(argv[0], argc, 3, obj);

    if(strcmp(argv[1],"f")==0) {
        cout<<"Forward ";
        e = &(p->cur->sWp);
    } else if(strcmp(argv[1],"b")==0) {
        cout<<"Backward ";
        e = &(p->cur->sWm);
    } else {
        throw parsefile_commandError("power: unrecognized direction "
            "of the excitation. It must be {b|f}.\n");
    }

    if(sscanf(argv[2], "%20lf", &threshold)!=1) {
        throw parsefile_commandError("power: error while reading "
            "the threshold.\n");
    }

    cout<<"Filtering modes with a real part of neff > "<<threshold<<endl;

    db_matrix excitation(p->cur->B.getNrow(),1);

    // Filter the modes that are not propagating.
    for(int i=0; i<p->cur->B.getNrow(); ++i) {

        // Obtain the effective index of the given mode
        complex<double>s=p->cur->B(i,i);
        s = p->getEffectiveIndex(s);

        // Check if the current mode is propagating or not.
        if (abs(s.real())>threshold) {
            excitation(i,0)=e->operator()(i,0);
        } else {
            excitation(i,0)=0.0;
        }
    }

    // Sum the contribution of all the existing modes. Those which are
    // not propagating are not being taken into account since their
    // excitation has been put equal to zero.
    db_matrix E=p->cur->W*excitation;
    db_matrix H=p->cur->V*excitation;

    power = section::integralPoynting(p,E, H,0)*p->tot_x*p->tot_y;;

    p->insertVar("ans",power);

    cout << "power in the current section (in W): "<<power<<endl;

    return 0;
}


/** POWERZ: calculates the power flowing in a given section
    This commands computes the power carried by all the modes propagating at a
    given z. All modes are taken into account. Measured power is considered
    positive when it is flowing in the forward direction, negative otherwise.
    If there is a flux of power in the two directions, the net result will
    be calculated.

    Usage:
    powerz z

    Parameters:
        z coordinate at which we compute the total power
*/
int commands::c_powerZ(parsefile *obj, int argc,char *argv[])
{
    structure *p;
    double power=0;
    double z=0;

    p = wantedParameters(argv[0], argc, 2, obj);

    if(sscanf(argv[1], "%20lf", &z)!=1) {
         throw parsefile_commandError("power: error while reading "
             "the z value.\n");
     }

    power = p->do_powerz(z);
    p->insertVar("ans",power);

    cout << "power at z = " << z << " (in W): "<<power<<endl;

    return 0;
}

/** MONITOR
    This function calculates the power flowing in a rectangular surface at
    a given depth z. The surface is transverse to the propagation axis z. The
    calculation is done as follows. The longitudinal component of the Poyinting
    vector is evaluated after the transverse electric and magnetic fields.
    Then, the power is calculated as a surface integral. This function can be
    useful in case of PMLs, since one can specify a rectangle which excludes
    them.

    Usage:
        monitor z wx wy px py

    Parameters:
        z   coordinate at which compute the total power
        wx  rectangle width
        wy  rectangle height
        px  x position of the center of the surface (referred to the center of
            the calculation window)
        py  y position of the center of the surface (referred to the center of
            the calculation window)
*/
int commands::c_monitor(parsefile *obj, int argc,char *argv[])
{
    structure *p;
    double z=0, power;
    double wx,wy,px,py;


    p = wantedParameters(argv[0], argc, 6, obj);

    // get the paramters:
    if(sscanf(argv[1], "%20lf", &z)!=1) {
         throw parsefile_commandError("monitor: error while reading "
             "the z value.\n");
     }
    if(sscanf(argv[2], "%20lf", &wx)!=1) {
         throw parsefile_commandError("monitor: error while reading "
             "the wx value.\n");
     }
    if(sscanf(argv[3], "%20lf", &wy)!=1) {
         throw parsefile_commandError("monitor: error while reading "
             "the wy value.\n");
     }
    if(sscanf(argv[4], "%20lf", &px)!=1) {
         throw parsefile_commandError("monitor: error while reading "
             "the px value.\n");
     }
    if(sscanf(argv[5], "%20lf", &py)!=1) {
         throw parsefile_commandError("monitor: error while reading "
             "the py value.\n");
     }

    power=p->do_monitor(z,wx,wy,px,py);
    p->insertVar("ans",power);

    return 0;
}


// NOTE: this should be in db_matrix instead?
db_matrix fft2cent(int dimx, int dimy, db_matrix dd)
{
    int nx = dd.getNcol();
    int ny = dd.getNrow();

    db_matrix res(dimy, dimx);
    int shx=(dimx-1)/2+1;
    int shy=(dimy-1)/2+1;

    int i,j;

    // positive x and y frequencies
    for(i=0; i<(dimy+1)/2; ++i) {
        for(j=0; j<(dimx+1)/2; ++j) {
            res(i,j)=dd(i,j);
        }
    }

    // positive x, negative y
    for(i=0; i<(dimy-1)/2; ++i) {
        for(j=0; j<(dimx+1)/2; ++j) {
            res(i+shy,j)=dd(i+ny-shy+1,j);
        }
    }

    // negative x, positive y
    for(i=0; i<(dimy+1)/2; ++i) {
        for(j=0; j<(dimx-1)/2; ++j) {
            res(i,j+shx)=dd(i,j+nx-shx+1);
        }
    }

    // negative x, negative y
    for(i=0; i<(dimy-1)/2; ++i) {
        for(j=0; j<(dimx-1)/2; ++j) {
            res(i+shy,j+shx)=dd(i+ny-shy+1,j+nx-shx+1);
        }
    }

    return res;
}

/** MODEPOS: to be called after SOLVE, this command loads in the "ans"
    variable the list of all indices of the interesting modes of the current
    section
    
    Usage:
    modeindex
    
    Return value:
    ans array loaded with the position indices
*/
int commands::c_modepos(parsefile *obj, int argc,char *argv[])
{
    structure *p;
    section *q;

    cout<<endl;
    
    p = wantedParameters(argv[0], argc, 1, obj);
    
    // Get the current section
    q=p->cur;
    
    // Check if the modes have been calculated for the current section.
     if (q->W.isEmpty()) {
        throw parsefile_commandError(
            "modeindex: the problem must be set and modes calculated with "
            "solve before assembling matrices for propagation.");
    }


    vector<double> results;
    vector<int> idx;
    q->find_interesting_g(results, idx);

    double *s =new double[idx.size()];
    
    int j=0;
    cout << "Indices are: ";    
    for (vector<int>::iterator it = idx.begin(); it != idx.end(); ++it) {
        cout<<*it<<" ";
        s[j++]=*it;
    }
    cout <<endl;
    p->nP->insertArray("ans", idx.size(), s);
    delete[] s;
    return 0;
}
