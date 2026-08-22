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

#include <cstring>
#include <fstream>
#include <iostream>

#include "material.h"
#include "structure.h"

using namespace std;


material::material()
{
	mat=NULL ;
    reset();
}

material::material(material &m)
{
    material::operator= (m);
}

/*material& material::operator= (const material &m)
{
    // Some ideas from there:
    // http://stackoverflow.com/questions/12423058

    // Prevent self-assignment
    if(this==&m)
        return *this;

    //material_data::operator= (m);

    return *this;
}*/

material::~material()
{
    if(mat!=NULL)
        delete[] mat ;
}

/** Increase the size of the material_data.

*/
void material::ensureSize(int newsize) {

    // Save the old pointers and size
    int old_size = actual_size;
    actual_size = newsize;
    material_data *oldmat = mat;

    // Allocate the new pointers
    mat = new material_data[actual_size];

    // Copy the old objects in the new position
    for (int i=0; i<old_size; ++i) {
        mat[i]=oldmat[i];
    }

    // Delete the old arrays
    delete[] oldmat;

}


void material::reset(void)
{
	if(mat!=NULL)
        delete[] mat;
	number_of_material = 0 ;
	actual_size = INIT_SIZE*0+1 ;
    mat = new material_data[actual_size];
}

/* read refractive index from file
parameters:
	filename: path to refractive index file to read
	mult is a multiplicator factor to apply to the wavelength in the file to 
		convert it to meter
	name is the material name to be stored as a reference for further use.
	lineCountOnly is set to true to only count the number of lines and not store 
		into arrays.
	m is a pointer to the material data where to store read data

Two file format are supported (see https://refractiveindex.info/):
Both files start with a headerline with names of column
File format with 2 columns : wavelength nr possibly followed by wavelength ni
File format with 3 columns: wavelength nr ni

After reading, adds the new material in the database.
*/
void material::read_from_file(char* filename, double mult, 
	bool lineCountOnly, material_data *m)
{
    double wavelength, nr, ni ;
	int ncols=0, nread ;
	char c ;

    FILE *f=fopen(filename,"r");
    if(f==NULL)
        throw parsefile_commandError("Can not open input file:-(");
	
	m->nb_wavelength_r = 0;
	m->nb_wavelength_i= 0 ;
	
    // skip a line corresponding to header   
    if (fscanf(f, "%*[^\n]")!=0)
        throw parsefile_commandError("error in reading an empty line");
        
    // parse the first line to determine the number of columns
    bool data_found = false ;
    bool EOL = false ;
    c = getc(f);
    while (not EOL) {
    	data_found = false ;
		c = getc(f);

    	// check if there is data (not separators)
    	if (c != '\t' and c !=' ' and c != ',' and c != ';') {
    		data_found = true ;
    		// remove extra data (not separators)
			while(c != '\t' and c !=' ' and c != ',' and c != ';' 
					and c != '\n' and c != '\r')
				c = getc(f);			
		}

		// look for separator
		if (c == '\t' or c ==' ' or c == ',' or c == ';'){
			if (data_found)
				ncols++;
			data_found = false; 

			// remove extra space or tab:
			while (c == '\t' or c ==' ')
				c = getc(f);
		}
		if (c == '\n' or c=='\r')
			EOL = true ;

    }
	if (data_found)
		ncols++;
	//rewind(f) ;
    cout << "ncols = " << ncols << endl ;
	    
	if (ncols ==2) {
	    // read from file:
	    nread = ncols ;
    	while ( nread == ncols) {
    		nread = fscanf(f,"%lf%lf", &wavelength, &nr) ;
    		if (!lineCountOnly) {
    			m->wavelength_r[m->nb_wavelength_r] = wavelength*mult ;
    			m->nr[m->nb_wavelength_r] = nr ;
    		}
    		m->nb_wavelength_r++ ;
    		//cout << wavelength << " "<< nr << endl;

    	}
    	
    	// removing possible empty lines
    	c = getc(f) ;
		while(c == '\t' or c ==' ' or c == '\n' or c == '\r')
			c = getc(f) ;

		// removing header
	    if (fscanf(f, "%*[^\n]")==0) {
	    	// second series
			nread = ncols ;
			while ( nread == ncols) {
				nread = fscanf(f,"%lf%lf", &wavelength, &ni) ;
	    		if (!lineCountOnly) {
					m->wavelength_i[m->nb_wavelength_i] = wavelength*mult ;
					m->ni[m->nb_wavelength_i] = ni ;
				}
				m->nb_wavelength_i++ ;
			}
		}
    } else if  (ncols ==3) {
    	nread = ncols ;
    	while ( nread == ncols) {
    		nread = fscanf(f,"%lf%lf%lf", &wavelength, &nr, &ni) ;
    		if (!lineCountOnly) {
    			m->wavelength_r[m->nb_wavelength_r] = wavelength*mult ;
    			m->nr[m->nb_wavelength_r] = nr ;
    			m->ni[m->nb_wavelength_r] = ni ;
    		}
    		m->nb_wavelength_r++ ;
    		//cout << wavelength << " "<< nr << endl;

    	}
    	
    } else {
    	throw parsefile_commandError("material: error material file format not"
    			"valid.\nOnly lambda n or lambda nr ni column format supported.");
    }

}

/* create_from_file called from command

filename: path to refractive index file to read
mult is a multiplicator factor to apply to the wavelength in the file to convert
	it to meter
name is the material name to be stored as a reference for further use.
*/
void material::create_from_file(char* filename, double mult, char* name) {
	
	material_data *m ;
	
	// check if name already exist. If yes find the index to overide
	int idx = get_material_id(name) ;
	if (idx == -1) {
		m = &(mat[number_of_material]) ;
		number_of_material ++; 
        if(number_of_material + 1 >= actual_size) {
            ensureSize(actual_size + INIT_SIZE);
    		m = &(mat[number_of_material-1]) ;
        }
	} else {
		cout << "Overiding material '" << name << "'" << endl;
		m = &(mat[idx]) ; 
	}
	cout << m->nb_wavelength_r << endl ;
	
	// get the number of wavelength for real and imaginary refractive index:
	read_from_file(filename, mult, true, m) ;

	// allocate the arrays to store data:
	if (m->nb_wavelength_r > 0) {
		m->wavelength_r = new double[m->nb_wavelength_r] ;
		m->nr = new double[m->nb_wavelength_r] ;
		if (m->nb_wavelength_i == 0){
			m->ni = new double[m->nb_wavelength_r] ;
			for (int i; i < m->nb_wavelength_r ; i++)
				m->ni[i] = 0 ;
		}
		m->name = new char[strlen(name)+2] ;
		strcpy(m->name, name) ;
	}
	if (m->nb_wavelength_i > 0) {
		m->wavelength_i = new double[m->nb_wavelength_i] ;
		m->ni = new double[m->nb_wavelength_i] ;
	}
	
	// read file and populate the arrays:
	read_from_file(filename, mult, false, m) ;
		
	double nr, ni ;
	get_refractive_index(name, 500e-9, &nr, &ni);
		
	
}
/*
	return the material id if the name is associated with a material
	Otherwise return -1 (meaning not found):
*/
int material::get_material_id(char* name) {

	for (int i=0 ; i < number_of_material ; i++) {
		if (strcmp(name, mat[i].name)==0)
			return i ;
	}
	return -1 ;
}
/** return refractive index for the given material at the given wavelenth
Parameters:
	name: name of the material
	wavelength: wavelength at which to extract the refractive index
Return:
	nr and ni are the returned refractive index;
	N.B. ni is returned as a negative value as it is the convention used in 
		the code
	
*/
void material::get_refractive_index(char *name, double wavelength, 
		double *nr, double *ni) {
	int idx = get_material_id(name) ;
	if (idx == -1) {
    	throw parsefile_commandError("material: unkown material name.");
	}
	material_data *m = &(mat[idx]) ;
	
	// look for wavelength:
	int i = m->nb_wavelength_r/2 ; 
	for (int step=m->nb_wavelength_r/4; step>0 ; step=step/2) {
		if (m->wavelength_r[i] >= wavelength)
			i -=step;
		else
			i+=step ;
	}
	// make sure that i is pointing to the wavelength just before the targeted
	// wavelength
	while (i<m->nb_wavelength_r && m->wavelength_r[i]<wavelength)
		i++;
	while (i>1 && m->wavelength_r[i]>wavelength)
		i--;
	
	if (m->wavelength_r[i]>wavelength || m->wavelength_r[i+1]<wavelength) {
		char str [VARDIM];
		sprintf(str, "material: don't find wavelength of %g in %s", wavelength, name);
		throw parsefile_commandError(str);
	}
		
	*nr = (m->nr[i+1]-m->nr[i])/(m->wavelength_r[i+1] - m->wavelength_r[i])*
		(wavelength - m->wavelength_r[i]) + m->nr[i] ;
	
	if (m->nb_wavelength_i == 0){
		*ni = -((m->ni[i+1]-m->ni[i])/(m->wavelength_r[i+1] - m->wavelength_r[i])*
		(wavelength - m->wavelength_r[i]) + m->ni[i]) ;
	} else {
		// look for wavelength in wavelength_i:
		i = m->nb_wavelength_i/2 ; 
		for (int step=m->nb_wavelength_i/4; step>0 ; step=step/2) {
			if (m->wavelength_i[i] >= wavelength)
				i -=step;
			else
				i+=step ;
		}
		// make sure that i is pointing to the wavelength just before the targeted
		// wavelength
		while (i<m->nb_wavelength_i && m->wavelength_i[i]<wavelength)
			i++;
		while (i>1 && m->wavelength_i[i]>wavelength)
			i--;
		
		if (m->wavelength_i[i]>wavelength || m->wavelength_i[i+1]<wavelength) {
			string stringvar= "material: don't find wavelength  of " + to_string(wavelength) +
				"in "+name +"\n";
			throw parsefile_commandError(stringvar.c_str());
		}
			
		*ni = -( (m->ni[i+1]-m->ni[i])/(m->wavelength_i[i+1] - m->wavelength_i[i])*
		(wavelength - m->wavelength_i[i]) + m->ni[i]) ;
	}
	
}

/* update variables that store refractive index (i.e. material name)

Parameters:
	p: pointer to structure 
	wavelength, wavelenth at which we compute the refractive index
Output:
	For each material create a variable with a array with real and imaginary
		refractive index
*/
void material::update_variables(structure *p, double wavelength) {
	
	double *n =new double[2];
	cout << "Wavelength = " << wavelength << endl;
	for (int i=0 ; i < number_of_material ; i++) {
		// Save refractive index values as an array into variable with name.
		get_refractive_index(mat[i].name, wavelength, &n[0], &n[1]) ;
		//cout << "pushing " << n[0] << " & " << n[1] << " into " << mat[i].name <<endl; 
		p->nP->insertArray(mat[i].name, 2, n);
    }
    delete[] n;
}

/* list the material already registered*/
void material::print_material(void) {
	for (int i=0 ; i < number_of_material ; i++) {
		cout << i << " " << mat[i].name <<  " " << mat[i].nb_wavelength_r << " nr points and "<<mat[i].nb_wavelength_i << " ni points" << endl; 
		
		/*
		for (int j=0; j < mat[i].nb_wavelength_r ; j++){
			if (mat[i].nb_wavelength_i == 0)
				cout << mat[i].wavelength_r[j] << " " << mat[i].nr[j] << " " <<  mat[i].ni[j]<< endl;
			else
				cout << mat[i].wavelength_r[j] << " " << mat[i].nr[j] << endl;
			}
		for (int j=0; j < mat[i].nb_wavelength_i ; j++){
			cout << mat[i].wavelength_i[j] << " " << mat[i].ni[j] << endl;
		}*/
	
	}
}



