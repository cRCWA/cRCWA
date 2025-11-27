#include <stdio.h>
#include <iostream>
#include <iomanip>

#include <string>
#include <stdlib.h>
#include "compileinfo.h"


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
    //char machineFile[]="/proc/cpuinfo";
    int hash = hashFile(NULL);
    int key;

    cout << " ***************************************************************************\n"
         << " *      Aperiodic Fourier Modal Method full vectorial 3D propagation       *\n"
         << " *                                                                         *\n"
         << " *     Davide Bucci, CROMA     March 2008 - October 2013                   *\n"
         << " *     Jérôme Michallon, CROMA     May 2011 - January 2013                 *\n"
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
         << " *     MINATEC-Grenoble INP, 3, parvis Luis Neel                           *\n"
         << " *     38016, Grenoble CEDEX, France                                       *\n"
         << " *                                                                         *\n"
         << " *     bucci@minatec.grenoble-inp.fr                                       *\n"
         << " *     jerome.michallon@minatec.grenoble-inp.fr                            *\n"
         << " *                                                                         *\n"
         << " ***************************************************************************\n";
    cout <<endl<<"License registration utility. "
        <<endl<<endl;

    cout <<"Contact bucci@minatec.grenoble-inp.fr by specifying the\n";
    cout <<"ID of this computer: "<<hash<<endl;

    cout<<"Do you want to update or create the license repository? [y/N] ";
    int c=getchar();
    if(c=='y' || c=='Y') {
        cout << "Insert the license key: ";
        cin >> key;
        upgradeLicense(key);
    }

}