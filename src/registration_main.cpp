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