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
        cerr<<"Exc. in main."<<endl;
        cerr<<P.getMess()<<"\n";
        cerr.flush();
    }


}

void test(void)
{
    unsigned int row=5;
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
    CC.printMatrix();
}

int main (int argc, char * const argv[])
{
    cout << " ***************************************************************************\n"
         << " *      Aperiodic Fourier Modal Method full vectorial 3D propagation       *\n"
         << " *                            version 1.5                                  *\n"
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
         << " *     Davide Bucci, CROMA     March 2008 - current                        *\n"
         << " *     Jérôme Michallon, CROMA     May 2012 - February 2014                *\n"
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


