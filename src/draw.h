/***************************************************************************
*     CLASS draw                                                      *
*	  Jérôme Michallon									   				*
*                                                                          *
****************************************************************************
*     Provides a class of methods for easily draw structure and 			*
*     generate normal field													*
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

    Davide Bucci, 2008-2026
    Jérôme Michallon, 2012-2014 2025-2026
*/

#ifndef DRAW_H
#define DRAW_H

#include <iostream>
#include "block_matrix.h"
#include "parsefile.h"

using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#define NBP 512	// multiple de 2
#define NBPnf 257	// puissance de 2 +1*/

#define PI 3.141592654
#define NB_LAYERS_MAX 100
#define NB_EDGES_MAX 100
#define precision 1e-8


class draw {
private:

	// parsing utilities:
	static void parse_input(int argc,char *argv[], int &readCount, bool interactive,
		bool ignoreOptional, const char* text, char *out) ;
	static void parse_input(int argc,char *argv[], int &readCount, bool interactive,
		bool ignoreOptional, const char* text, int &out);
	static void parse_input(int argc,char *argv[], int &readCount, bool interactive,
		bool ignoreOptional, const char* text, double &out);
	static void parse_skip_optional(int argc,char *argv[], int &readCount)	;
	// subfunction for perparing points
	void treat_input_point(int layerID, bool SortAngle) ;
	void treat_input_polygone_circle(int layerID, double alpha0);
	void treat_input_rectangle(int layerID) ;
	
	// subfunctions for normal field generation
	bool contained_in_polygone(int layerID, double x, double y) ;
	
	int interpolateNF ( int centerx[NB_LAYERS_MAX], int centery[NB_LAYERS_MAX],
		int Polygonex[NB_EDGES_MAX][NB_LAYERS_MAX] , int Polygoney[NB_EDGES_MAX][NB_LAYERS_MAX],
		int i1,int j1,int i2,int j2,int i3,int j3,int i4,int j4,int m,int n) ;

	// generation of index and normal field
	void normal_field_generation(void) ;
	void polygone(void ) ;


	// function to save to file
	void writeoutputfield(db_matrix &out, char* header, bool optiwaveformat, 
		FILE *f) ;
	void writeoutputfield(db_matrix &out, char* header, bool optiwaveformat,
		char *filename) ;
	void generate_command(char *command) ;
	


	// STORE ALL INTERNAL DATA:
	
	// x- and y-coordinate of the points of the polygone 
	double x_poly[NB_EDGES_MAX][NB_LAYERS_MAX] ;
	double y_poly[NB_EDGES_MAX][NB_LAYERS_MAX] ;
	
	// file name where to save files for standalone:
	char filename_ind[200];
	
	// number of layer used in the structure
	int nb_layers;
	
	// shape type for each layer:
	char shape_type[NB_LAYERS_MAX] ;
		
	// number of points that define the shape (3=triangle, 4=rectangle, 6=hexagone...)
	int nb_edges[NB_LAYERS_MAX];
		
	// x-and y-coordinate of the center of the polygone shape
	double origineX[NB_LAYERS_MAX];
	double origineY[NB_LAYERS_MAX];
	
	// size of the computation window in x and y directions or ratio of sy/sx with sx=1
	double sizeY;
	double sizeX ;
	
	//orientation of the normal field along X, Y direction (toward inside the structure or toward a direction)
	double OrientationX;
	double OrientationY;
	//x- and y-coordinate of the origin of the normal field
	double CentreX;
	double CentreY ;
	
	// radius or ratio radius / sx where sx and is the size of the computation window in x direction
	double r_a[NB_LAYERS_MAX];
	// radius or ratio radius / sy where sy and is the size of the computation window in y direction
	double r_b[NB_LAYERS_MAX];
	
	// width_x or ratio width_x /sx where sx and is the size of the computation window in x direction
	double wx[NB_LAYERS_MAX] ;
	// width_y or ratio width_y /sy where sy and is the size of the computation window in y direction
	double wy[NB_LAYERS_MAX] ;
	
	// origin angle for the first point of a shape (rotation of the shape)
	double alpha0[NB_LAYERS_MAX] ;
	
	// real and imaginary part of the refractive index of the shape
	double n_int[NB_LAYERS_MAX];
	double k_int[NB_LAYERS_MAX];
	
	// real and imaginary part of the refractive index of surronding medium
	double n_ext;
	double k_ext;
	
	// Matrices with the x- and y- component of the normal field
	db_matrix Nx;
	db_matrix Ny ; 
	
	// matrix in which we write the complex refractive indices
	db_matrix ind ;
	
	// Set to true when we only compute refractive index and not normal field
	bool OnlyIndex ;
	
protected:
	//void init(structure &s);
    
public:
    
    // Standard constructor
    draw();

    /* Here follows callable draw functions */
    // parsing functions
    int parse_legacy(int argc,char *argv[]);
    int parse(int argc,char *argv[], bool interactive);
    int parse_standalone(int argc,char *argv[]);

	// getter and setter for incorporating in crcwa:
	db_matrix &get_refractive_index(void) { return ind;	}
	void set_size(double sx, double sy) {sizeX=sx ; sizeY=sy ; }
	bool get_OnlyIndex(void) {return OnlyIndex ;};
	db_matrix &get_Normal_Field_x(void) { return Nx;	}
	db_matrix &get_Normal_Field_y(void) { return Ny;	}
	
	//void saveToFile(void) ;
	
	// information:
	void help(void) ;
};

#endif
