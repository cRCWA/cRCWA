/***************************************************************************
*     CLASS MATERIAL                                                        *
****************************************************************************
*    C++ class containing methods for refractive index variation 			*
*			with wavelength 												*
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
    Jérôme Michallon, 2012-2014 & 2026
*/

// Prevent multiple includes
#ifndef MATERIAL_H
#define MATERIAL_H

#include "block_matrix.h"

// This is needed since material contains a pointer to the structure class
// We do not need to define the class here, but just provide a placeholder
class structure;

class material_data {
public:
	
	// arrays for wavelength, n and k
	double * wavelength_r;
	double *nr ;
	double *ni ;
	int nb_wavelength_r ;
	// wavelengthi is used for ni when nr and ni are not on the same wavelength
	double *wavelength_i ;	
	int nb_wavelength_i ;
	
	// material name
	char* name ;
	
};

class material {
private:
	// array of material:
	material_data *mat ;
	int number_of_material ;
	int actual_size ;
	
    // read refractive index from file:
    void read_from_file(char* filename, double mult, 
        	bool lineCountOnly, material_data *m);
    
   
    // display some information on already registered materials:
    void print_material(void) ;
    
    // Increase the size of material_data array.
    void ensureSize(int newsize) ;
    
    // return refractive index for a given material at a given wavelength
    void get_refractive_index(char *name, double wavelength, 
		double *nr, double *ni) ;
		
    // get data_material index for a given material name
    int get_material_id(char* name) ;  
public:

    // Standard constructor
    material();
    // Copy constructor
    material(material &m);
    // Copy operator
    //material& operator= (const material &m);

    // Destructor
    ~material();

    // Resets the current section
    void reset(void);
    
    // public functions
    
    // create a new material from file 
    void create_from_file(char* filename, double mult, char* name) ;

	// update material variables:
	void update_variables(structure *p, double wavelength) ;
		    

};

#endif
