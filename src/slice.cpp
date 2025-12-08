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

#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

typedef enum{XY, XZ, YZ, S} SliceType;

/** Extracts the xy slice from the field matrix

*/
void slice(SliceType s, char *fin_name, double quota, double tol,
    char *fout_name, bool lx, bool ly, bool lz)
{
    FILE *fin=fopen(fin_name,"r");
    FILE *fout=fopen(fout_name,"w");
    char line[128];

    if(fin==NULL) {
        cerr << "Unable to open the input file "<< fin_name<<endl;
        exit(1);
    }
    if(fout==NULL) {
        cerr << "Unable to open the output file "<< fout_name<<endl;
        exit(1);
    }

    // skip a line (IS IT ALWAYS NEEDED THIS???)

    //fscanf(fin, "%*[^\n]");

    double x, y, z;
    double oldx, oldy, oldz;

    // set a very very large default value
    oldx=oldy=oldz = 1e307;

    bool ch = false;
    bool firstLine=true;
    int numberCounter =0;
    int numberValue =-1;

    fprintf(fout, "# x          y            z             v\n");
    while(true) {
        if(fgets(line, sizeof line, fin) == NULL) {
            cout<<"End of file"<<endl;
            break;
        }
        if(sscanf(line,"%20lf%20lf%20lf",&x,&y,&z)!=3) {
            break;
        }
        switch (s) {
            case XY:
                ch=(abs(z-quota)<tol);
                break;

            case XZ:
                ch=(abs(y-quota)<tol);
                break;

            case YZ:
                ch=(abs(x-quota)<tol);
                break;

            case S:
                cerr<<"Uh???"<<endl;
                return;
        }

        if (ch) {
            if (numberValue<0 && ((lx && oldx!=x) || (ly && oldy!=y) ||
                (lz && oldz!=z))) {
                if(!firstLine) {
                    fprintf(fout, "\n");

                    numberValue=numberCounter;
                    cout <<"Scanline found: "<< numberValue<<" points."<<endl;
                    numberCounter=0;
                } else {
                    firstLine=false;
                    oldx=x;
                    oldy=y;
                    oldz=z;
                }
            } else if(numberCounter==numberValue) {
                numberCounter =0;
                fprintf(fout, "\n");
            }

            ++ numberCounter;
            fprintf(fout, "%s",line);

            oldx=x;
            oldy=y;
            oldz=z;
        }

    }

    fclose (fin);
    fclose (fout);
    cout << fout_name << " has been written.\n";

}

/** Extracts the wanted propagation step from the field matrix

*/
void propstep(char *fin_name, char *fout_name, int step)
{
    FILE *fin=fopen(fin_name,"r");
    FILE *fout=fopen(fout_name,"w");

    if(fin==NULL) {
        cerr << "Unable to open the input file "<< fin_name<<endl;
        exit(1);
    }
    if(fout==NULL) {
        cerr << "Unable to open the output file "<< fout_name<<endl;
        exit(1);
    }

    // skip a line

    //fscanf(fin, "%*[^\n]");

    double x, y, z;
    char line[128];
    double oldx, oldy, oldz;

    // set a very very large default value
    oldx=oldy=oldz = 1e307;

    bool firstLine=true;
    int number_of_points=0;
    int counter=0;
    int step_counter=0;
    bool ch=false;
    bool startpoint=true;


    fprintf(fout, "# x          y            z             v\n");
    while(true) {
        if(fgets(line, sizeof line, fin) == NULL) {
            break;
        }

        if(sscanf(line,"%20lf%20lf%20lf",&x,&y,&z)!=4) {
            break;
        }
        if(startpoint) {
            oldz=z;
            startpoint=false;
        }
        ch=false;
        if(firstLine) {
            if(oldz!=z) {
                cout <<"Scanline found: "<<number_of_points<<" points."<<endl;
                firstLine = false;
            } else {
                ++number_of_points;

                if (step==0)
                    ch=true;
            }
            counter=0;
        } else {
            if(counter>number_of_points) {
                ++step_counter;
                cout << ".";
                cout.flush();
                counter=0;
            }
        }
        if(step_counter==step)
            ch=true;
        if(step_counter>step){
            cout<<endl;
            break;
        }

        if (ch)
            fprintf(fout, "%s",line);
        ++counter;
        oldz=z;
    }

    fclose (fin);
    fclose (fout);
    cout << fout_name << " has been written.\n";
}

int main (int argc, char * const argv[])
{
    cout << " ***************************************************************************\n"
         << " *      Aperiodic Fourier Modal Method full vectorial 3D Mode Solver       *\n"
         << " *                             version 1.3.5                               *\n"
         << " *                            (slicer utility)                             *\n"
         << " *                                                                         *\n"
         << " *     Davide Bucci, CROMA     march 2008 - march 2012                     *\n"
         << " *     MINATEC-INPG, 3, parvis Luis Neel                                   *\n"
         << " *     38016, Grenoble CEDEX, France                                       *\n"
         << " *                                                                         *\n"
         << " *     bucci@minatec.grenoble-inp.fr                                       *\n"
         << " *                                                                         *\n"
         << " ***************************************************************************\n\n";

    if (argc < 6) {
        cerr << "Not enough parameters on the command line\n";
        cerr << "Usage: "<<argv[0]<<" [-lx] [-ly] [-lz] {xy|xz|yz} quota tolerance input_file output_file\n";
        return 1;
    }

    double quota;

    int eat_options = 1;

    bool lx = false;
    bool ly = false;
    bool lz = false;


    while(argv[eat_options][0]=='-') {
        if(strcmp(argv[eat_options],"-lx")==0) {
            lx = true;
            cout << "Line skip each x step.\n";
        } else if(strcmp(argv[eat_options],"-ly")==0) {
            ly = true;
            cout << "Line skip each y step.\n";
        } else if(strcmp(argv[eat_options],"-lz")==0) {
            lz = true;
            cout << "Line skip each z step.\n";
        } else {
            cerr << "Unrecognized option.\n";
            return 1;
        }
        ++eat_options;
    }

    int type_arg = eat_options;
    int quota_arg = type_arg+1;
    int tol_arg = quota_arg+1;
    int inp_file_arg = tol_arg+1;
    int out_file_arg = inp_file_arg+1;

    double tolerance;
    int propagation_step=-1;


    SliceType s=XY;
    if(strcmp(argv[type_arg],"xy")==0) {
        s = XY;
        cout << "z = ";
    } else if(strcmp(argv[type_arg],"xz")==0) {
        s = XZ;
        cout << "y = ";
    } else if(strcmp(argv[type_arg],"yz")==0) {
        s = YZ;
        cout << "x = ";
    } else if(strcmp(argv[type_arg], "s")==0) {
        s = S;
        inp_file_arg = type_arg+2;
        out_file_arg = inp_file_arg+1;
        if(sscanf(argv[type_arg+1], "%20d", &propagation_step)<1) {
            cerr << "I could not read the wanted propagation step number.\n";
            return 1;
        }
        cout << "Search for propagation step "<<propagation_step<<endl;
    } else {
        cerr << "Slice "<<argv[type_arg]<< " is not recognized.\n";
        cerr << "Allowed slices are {xy|xz|yz}\n";
        return 1;
    }

    if(propagation_step<0) {
        cout << "Slicing "<< argv[inp_file_arg] << " at quota ";

        if(sscanf(argv[quota_arg], "%20lf", &quota)<1) {
            cerr << "I could not read the wanted quota.\n";
            return 1;
        }

        if(sscanf(argv[tol_arg], "%20lf", &tolerance)<1) {
            cerr << "I could not read the wanted tolerance.\n";
            return 1;
        }
        cout << quota << " m with tolerance "<<tolerance << " m\n";
        slice(s, argv[inp_file_arg], quota, tolerance, argv[out_file_arg],
            lx, ly, lz);
    } else {
        propstep(argv[inp_file_arg], argv[out_file_arg],propagation_step);
    }

    return 0;
}
