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
#include <complex>
#include <stdio.h>

typedef enum rimco_e_t{R,I,M,C,O} rimco_e;

using namespace std;

class commandError {
    string errMess;
public:

    commandError(){errMess="Error";}
    commandError(string m){errMess=m;}
    ~commandError(){}

    string getMess(){return errMess;}
};

/** Read the contents of a file in the Optiwave format and store it in a
    matrix which will have the number of points defined in the file.
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
    @param filename the path and the name of the file to be opened.
    @param mult the multiplier to be used for all sizes specified in the file
        (this is useful since Optiwave often specifies all dimensions in
        micrometers whereas AFMM uses meters).
    @return the matrix containing the points specified in the file.
*/
complex<double> *read_optiwave(char* filename, int &nx, int &ny,
    double &xmin, double &xmax, double &ymin, double &ymax, double mult)
{

    double tot_x;
    double tot_y;
    int i,j;


    FILE *f=fopen(filename,"r");
    if(f==NULL)
        throw commandError("Can not open input file :-(");

    // skip a line
    fscanf(f, "%*[^\n]");

    if(fscanf(f,"%20d%20d",&nx,&ny)!=2) {
        fclose(f);
        throw commandError("Unable to read the number of x and y points in the"
            " given file.");
    }

    if(fscanf(f,"%20lf%20lf%20lf%20lf",&xmin,&xmax, &ymin, &ymax)!=4) {
        fclose(f);
        throw commandError("Unable to read the x and y ranges in the given"
        " file.");
    }

    xmin*=mult;
    xmax*=mult;
    ymin*=mult;
    ymax*=mult;
    tot_x=(xmax-xmin);
    tot_y=(ymax-ymin);

    cout << "The file contains "<<nx<< " x  " << ny <<
            " points.\n";
    cout << "The x range is "<<tot_x << " m.\n";
    cout << "The y range is "<<tot_y << " m.\n";

    complex<double> *lect= new complex<double>[ny*nx];

    double re, im;
    double norm = 0;
    int ch;

    for(i=0; i<ny; ++i) {
        for(j=0; j<nx; ++j) {
            while ((ch=fgetc(f))=='('||ch==')'||ch=='\n'||ch=='\r')
                /*does nothing*/;
            if(ch==EOF) {
                fclose(f);
                throw commandError("File terminated abruptly.");
            }

            ungetc(ch, f);
            if(fscanf(f,"%20lf,%20lf",&re,&im)!=2) {
                fclose(f);
                throw commandError("Unable to read a value in the input file.");
            }
            lect[i*nx+j]=complex<double>(re,im);
        }
    }

    fclose(f);

    return lect;
}

void write_gnuplot(char *filename, complex<double> *out, int nx, int ny,
    double dx, double dy, double x0, double y0, enum rimco_e_t rimc)
{
    double x;
    double y;

    FILE *f=fopen(filename,"w");
    if(f==NULL)
        throw commandError("Can not open output file :-(");
    cout<<"nx="<<nx<<"  ny="<<ny<<endl;
    for(int i=0; i<ny; ++i) {
        for(int j=0; j<nx; ++j) {

            x = j*dx-(nx-1)/2.0*dx+x0;
            y = i*dy-(ny-1)/2.0*dy+y0;

            // Output
            if(rimc==C) {
                fprintf(f, "%le %le %le %le\n",
                    x, y, out[i*nx+j].real(),out[i*nx+j].imag());
            } else if(rimc==R) {
                fprintf(f, "%le %le %le\n",
                    x, y, out[i*nx+j].real());
            } else if(rimc==I) {
                fprintf(f, "%le %le %le\n",
                    x, y, out[i*nx+j].imag());
            } else if(rimc==M) {
                fprintf(f, "%le %le %le\n",
                    x, y, abs(out[i*nx+j]));
            }
        }
        fprintf(f,"\n");
    }
    fclose(f);
}

int main(int argc, char **argv)
{
    if (argc<3) {
        cerr << "Optiwave to Gnuplot converter. \nUsage "<<
            argv[0]<< " file_in.optiwave file_out.gnuplot"<<endl;
        return 1;
    }
    double xmin, xmax, ymin, ymax;
    int nx, ny;
    try{
        // First read the matrix in the Optiwave format
        complex<double> *A = read_optiwave(argv[1], nx, ny, xmin, xmax, ymin,
            ymax, 1.0);
        // Then write the contents in the Gnuplot format
        write_gnuplot(argv[2], A, nx, ny, (xmax-xmin)/nx, (ymax-ymin)/ny,
            (xmax+xmin)/2.0, (ymax+ymin)/2.0, C);
    } catch (commandError E) {
        cerr << "Error: "<<E.getMess()<<endl;
        return 1;
    }
    return 0;
}

