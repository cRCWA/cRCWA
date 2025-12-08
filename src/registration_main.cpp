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

#include <stdio.h>
#include <iostream>
#include <string>
#include <stdlib.h>


#include "lock_processor.h"

using namespace std;

void upgradeLicense(int key)
{
    FILE *file;
    char * home = getenv ("HOME");
    string fn= string(home)+"/.afmm_registry";
    file = fopen(fn.c_str(), "w");
    if (file==NULL) {
        cerr<<"Can not write on the license file."<<endl;
        return;
    }
    fprintf(file, "%d", key);
    fclose(file);
    cout<<"License file correctly updated."<<endl;
}

int main(void)
{
    char machineFile[]="/proc/cpuinfo";
    int hash = hashFile(machineFile);
    char passphrase[] = {118, 97, 110, 105, 108, 108, 97,0};
    int key=getKey(passphrase);

    cout << " ***************************************************************************\n"
         << " *      Aperiodic Fourier Modal Method full vectorial 3D propagation       *\n"
         << " *                                                                         *\n"
         << " *     Davide Bucci, CROMA     March 2008 - January 2013                   *\n"
         << " *     Jérôme Michallon, CROMA     May 2011 - January 2013                 *\n"
         << " *     MINATEC-Grenoble INP, 3, parvis Luis Neel                           *\n"
         << " *     38016, Grenoble CEDEX, France                                       *\n"
         << " *                                                                         *\n"
         << " *     bucci@minatec.grenoble-inp.fr                                       *\n"
         << " *     jerome.michallon@minatec.grenoble-inp.fr                            *\n"
         << " *                                                                         *\n"
         << " ***************************************************************************\n";
    cout <<endl<<"License registration utility. This utility should not be "
        "distributed along\n with the executable (and locked) version of afmm."
        <<endl<<endl;

    FILE *file= fopen(machineFile, "r");
    if (file==NULL) {
        cout<<"WARNING: the machine description file "<<machineFile<<
            " is not available."<<endl;
        cout<<"         Afmm can still be licensed, but the license file will "
            "not depend"<<endl
            <<"         upon a particular machine."<<endl<<endl;
    } else {
        fclose(file);
    }

    cout <<"Key on this computer: "<<key<<endl;

    if(checkKey(passphrase)) {
        cout<<"License check passed."<<endl;
    } else {
        cout<<"The license file is absent, or it does not contain the correct "
            "license code"<<endl;

        cout<<"Do you want to update or create the license file? [y/N] ";
        int c=getchar();
        if(c=='y' || c=='Y') {
            upgradeLicense(key);
        }
    }
}