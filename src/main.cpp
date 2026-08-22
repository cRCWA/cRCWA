/***************************************************************************
*     FMM 3D full vectorial mode solver and propagator                     *
*     Davide Bucci, CROMA                                                  *
*     Jérôme Michallon, CROMA                                              *
*     MINATEC-INPG, 3, parvis Louis Neel                                   *
*     38016, Grenoble CEDEX, France                                        *
*                                                                          *
*     bucci@minatec.grenoble-inp.fr                                        *
*                                                                          *
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

using namespace std;

#include "structure.h"
#include "block_matrix.h"
#include "commands.h"
#include "compileinfo.h"



void process(int argc, char * const argv[])
{
    numParser nP;

    structure waveguide(&nP);
    commands co(waveguide, &nP);
    bool externalCommands = false;
    int fileIndex = 1;

    // Let's check if we need to allow the execution of external commands.
    // This is given by the "-e" option

    if(argc>1 && strcmp(argv[1],"-e")==0) {
        externalCommands=true;
        co.allow_system_command(externalCommands);
        cout<<"Execution of external commands allowed."<<endl;
        ++fileIndex;
    }

    // We determine if the user wants to enter in the interactive mode by
    // watching if he has specified the name of the file to be processed or
    // not.

    try {
        if (argc>fileIndex) {
            string fileName = argv[fileIndex];
            co.read(fileName, false);
        } else
            co.read("stdin (interactive mode)", true);

    } catch(parsefile_commandError P) {
        cerr<<endl <<"Exc. in main."<<endl;
        cerr<<P.getMess()<<"\n";
        cerr.flush();
    }


}
static double sinc(double alpha) {return (alpha==0)?1.0:sin(alpha)/alpha;}
	double tot_x = 1e-6 ;
	double tot_y = 1e-6 ;
db_matrix transfPML(bool isX, double q, double qgammar, double qgammai) {
    int ii, jj;
    int i, j;

    double d;               // period (size) of the window
    //double q;               // size of the mapped region
    int Nx;                 // Size of the resulting matrix
    int Ny;

    int k;
    bool tozero=false;

    Nx=256 ; //father->dimx;
    Ny=256 ; //father->dimy;

	complex<double> gamma=complex<double>(qgammar, qgammai);
	
    if(isX) {
        d=tot_x;
        //q=qx;
    } else {
        d=tot_y;
        //q=qy;
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

            if (isX){  // x case
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
                    
                //A(ii,jj)= 1-exp(-pow(k-0.5*d,2)/(q*q))-exp(-pow(k+0.5*d,2)/(q*q)) ;
                                  
        
            }
       }
    }

    return A;
}
void output(const char * filename, db_matrix &out) {
	double x,y ;
    double dx=tot_x / out.getNcol();
    double dy=tot_y / out.getNrow();
    FILE *f=fopen(filename,"w");
    for(int i=0; i<out.getNrow(); ++i) {
        for(int j=0; j<out.getNcol(); ++j) {
            x = j*dx-(out.getNcol()-1)/2.0*dx;
            y = i*dy-(out.getNrow()-1)/2.0*dy;

			fprintf(f, "%le %le %le %le\n",
				x, y, out(i,j).real(),out(i,j).imag());
		}
		fprintf(f, "\n");
	}
fclose(f);  
}

/** Specify that for the current section, a coordinate-transform PML strategy
    has to be applied. See J. Hugonin, P. Lalanne, I. Villar, and I. Matias,
    "Fourier modal methods for modeling optical dielectric waveguides,"
    Optical and quantum electronics, vol. 37, no. 1, pp. 107-119, 2005.*/
void test(void)
{
	section *sec = new section ;
	bool isX = true ;
	db_matrix A=transfPML(isX, 50e-9, 0.5, -0.5) ;
	db_matrix AA ;
	db_matrix x_coord ;
	db_matrix y_coord ;

	// the coordinate change are stored:
	int Mx=257;
	int My=257;
cout << "computing x_coord:\n";
	// the coordinate change are computed with:
	// int 1/A(ii,jj)
	if(isX) {
// JM_iFFT
		AA = A.zero_pad(1,Mx).ifft2();
		x_coord = AA;		
        double dx=tot_x / Mx;
		int shiftx=Mx/2 ;
		int i=0 ;
		//for(int i=0; i<My; ++i) {

			x_coord(i,shiftx) = 0 ;
		    for(int j=1; j<=shiftx; ++j) {
				x_coord(i,j+shiftx) = x_coord(i,j+shiftx-1) + 1.0/AA(i,j+shiftx)*dx;
				x_coord(i,shiftx-j) = x_coord(i,shiftx-j+1) - 1.0/AA(i,shiftx-j)*dx;
			}
		//}
	} else {
// JM_iFFT
		AA = A.zero_pad(My,Mx).ifft2();
		y_coord = AA;	
        double dy=tot_y / My;
		int shifty=My/2 ;
		for(int j=0; j<Mx; ++j) {
			y_coord(shifty,j) = 0 ;
		    for(int i=1; i<=shifty; ++i) {
				y_coord(shifty+i,j) = y_coord(shifty+i-1,j) + 1.0/AA(shifty+i,j)*dy;
				y_coord(shifty-i,j) = y_coord(shifty-i+1,j) - 1.0/AA(shifty-i,j)*dy;
			}
		}
	}

	// save to file:
	output("FFT_der.txt", A) ;
	output("derm1.txt", AA) ;
	output("x_coord.txt", x_coord) ;
	
    /*unsigned int row=5;
    unsigned int col=5;
    
    db_matrix A(row,col);

    for(unsigned int i=0;i<row;++i)
        for (unsigned int j=0; j<col;++j)
            A(i,j)=complex<double>(i,j);

    A.printMatrix();
    cout<<endl;
    A.fft2().printMatrix();
    cout<<endl;
    unsigned int nr=7;
    unsigned int nc=7;
    db_matrix BB=A.zero_pad(nr,nc);
    BB.printMatrix();
    cout<<endl;
    db_matrix DD=A.fft2().zero_pad(nr,nc);
    DD.printMatrix();
    cout<<endl;
    db_matrix CC=DD.ifft2()*(1.0/row/col);
    CC.printMatrix();*/
    
    /*FILE* pFin=fopen("test.fmm", "rb");
    // Can not open the given file.
    if (pFin==NULL){
		cerr << "File not found\n";
    }
    fseek(pFin,0,0);
    int line=0;
    char c='a' ;
    char string[256] ;
    int i=0 ;
    long int position[10] ;
    char buffer[BUFDIM+1];
    bool cont ;
    while (line<= 5) {
		position[line]=ftell(pFin);
		cont = (fgets(buffer, BUFDIM, pFin)!=NULL);
		cout << "at " << position[line] << " " <<line << "> " << buffer << endl ;
		line++;
    }
    cout << "reading back\n";
    i=1 ;
    fseek(pFin, position[i],0);
    cont = (fgets(buffer, BUFDIM, pFin)!=NULL);
    cout << "at " << position[i] << " " <<i << "> " << buffer << endl ;
    
    i=4 ;
    fseek(pFin, position[i],0);
    cont = (fgets(buffer, BUFDIM, pFin)!=NULL);
    cout << "at " << position[i] << " " <<i << "> " << buffer << endl ;
    
    
    i=3 ;
    fseek(pFin, position[i],0);
    cont = (fgets(buffer, BUFDIM, pFin)!=NULL);
    cout << "at " << position[i] << " " <<i << "> " << buffer << endl ; */
}

int main (int argc, char * const argv[])
{
    cout << " ***************************************************************************\n"
         << " *      Aperiodic Fourier Modal Method full vectorial 3D propagation       *\n"
         << " *                            version 1.5.4                                *\n"
         << " *                                                                         *\n"
         << " *     Build date: " << __DATE__<<  "                                             *\n"
         << " *     Source revision: "
         << setw(30)
         << setfill(' ')
         << left
         << _SVNGLOBALVERSION
         << std::setw(0)
         << "                     *\n"
         << " *                                                                         *\n"
         << " *     Davide Bucci, CROMA March 2008 - current                            *\n"
         << " *     Jérôme Michallon, CROMA May 2012 - February 2014                    *\n"
         << " *     MINATEC-Grenoble INP, 3, parvis Louis Neel                          *\n"
         << " *     38016, Grenoble CEDEX, France                                       *\n"
         << " *                                                                         *\n"
         << " *     davide.bucci@grenoble-inp.fr                                        *\n"
         << " *                                                                         *\n"
         << " ***************************************************************************\n";
    init_semaphore_FFTW();
    process(argc, argv);
    //test();
    //db_matrix::leaks();

    delete_semaphore_FFTW();
    return 0;
}


