/***************************************************************************
*     Draw utility                 	                                        *
*	  Jérôme Michallon		2012-20							   				*
*                                  	                                        *
****************************************************************************

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

    Davide Bucci, 2008-2026
    Jérôme Michallon, 2012-2014 2025-2026
*/
#include "draw.h" 
#include "block_matrix.h"


int main (int argc,char *argv[])
{

    if (argc<2) {
        cout << "Draw Utility part of cRCWA sofware suite. \n" <<
        		"Use ./draw --help for usage options or answer the following questions:\n "<<endl;
   }
    

	draw dw;
    try {
		if (argc >1) {
			if (strncmp(argv[1], "--legacy",8)==0){
				dw.parse_legacy(argc-1,&argv[1]) ;
			} else if (strncmp(argv[1], "--help",6)==0){
				dw.help();
			} else 
				dw.parse_standalone(argc,argv) ;
		} else
			dw.parse_standalone(argc,argv) ;
    } catch(parsefile_commandError P) {
    cerr<<endl <<"Exc. in draw."<<endl;
    cerr<<P.getMess()<<"\n";
    cerr.flush();
    }
	

}
