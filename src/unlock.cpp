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


// #include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string>

#include <unistd.h>


#include <stdlib.h>
#include "compileinfo.h"

#include "lock_processor.h"

using namespace std;

int main(void)
{
    int key;
    int hash;

    string passphrase;


    cout << " ***************************************************************************\n"
         << " *      Aperiodic Fourier Modal Method full vectorial 3D propagation       *\n"
         << " *                                                                         *\n"
         << " *     Davide Bucci, CROMA     March 2008 - present                        *\n"
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
         << " *                                                                         *\n"
         << " ***************************************************************************\n";
    cout <<endl<<"License registration utility. "
        <<endl<<endl;

    cout << "Insert the computer ID: ";
    cin >> hash;
    char *pass=getpass("Insert the passphrase: ");
    if(pass==NULL) {
        cerr << "Could not read the passphrase"<<endl;
        return 1;
    }

    passphrase += pass;

    cout << "Access code is: "<<getKey(hash, passphrase.c_str())<<endl;

    return 0;
}