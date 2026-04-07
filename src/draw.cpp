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


#include "draw.h"
#include "parsefile.h"

using namespace std;

// Standard constructor: fabricate a null matrix in which the A pointer is not
// defined.
draw::draw(void)
{
	sizeX = 1 ;
	sizeY = 1 ;
	OnlyIndex = false ;
	strcpy(filename_ind, "filename.rid");
}

double abs2(double a)
{
	if (a<0)
	{
		return (-a);
	}
		return (a);
}



// check if the point x,y is contained in the polygone pointed by layerID
// layerID is the layer to where we want to know if the point x,y belong to it
// x an y are coordinate of the point
bool draw::contained_in_polygone(int layerID, double x, double y) {

	double point1[2], point2[2] ;
	
	double angles_sum = 0 ;
	for (int k=0; k<nb_edges[layerID] ; k++) {
		// 0.5 is added to NBP to get a rounding
		
		point1[0] = origineX[layerID]*NBP + x_poly[k][layerID] *(NBP+0.5) - x;
		point1[1] = origineY[layerID]*NBP + y_poly[k][layerID] *(NBP+0.5) - y;

		point2[0] = origineX[layerID]*NBP + x_poly[(k+1)%nb_edges[layerID]][layerID] *(NBP+0.5) - x ;
		point2[1] = origineY[layerID]*NBP + y_poly[(k+1)%nb_edges[layerID]][layerID] *(NBP+0.5) - y ;
		
		angles_sum += atan2( (point1[0]*point2[1] -point1[1]*point2[0]),
			(point1[0]*point2[0] + point1[1]*point2[1]) ) ;
	}
	return (abs2(angles_sum - 2*PI) <= precision ) ;
}	

// interpolateNF the normal field and update Nx and Ny, which are the normal 
// field in x and y direction
int draw::interpolateNF (int centerx[NB_LAYERS_MAX], int centery[NB_LAYERS_MAX], 
	int Polygonex[NB_EDGES_MAX][NB_LAYERS_MAX] , 
	int Polygoney[NB_EDGES_MAX][NB_LAYERS_MAX],
	int i1,int j1,int i2,int j2,int i3,int j3,int i4,int j4,
	int m,int n)
{

	int nbp, k,l,o,  i,j ;
	double alpha ;
	Nx(n,m) = 0;
	Ny(n,m) = 0;
	double norm ;
	double ax,bx,cx ;
	double ay,by,cy ;

	// weight the point from the radius:
	// for each shape
	for (l = 0 ; l < nb_layers ; l++ )
	{
		// for each side of the shape
		for (o = 0 ; o < nb_edges[l] ; o++)
		{

			// vector a along the side
			ax = Polygonex[(o+1)%nb_edges[l]][l] -Polygonex[o][l] ;
			ay = Polygoney[(o+1)%nb_edges[l]][l] -Polygoney[o][l] ;

			// vector b perpendicular pointing to outside of the calculation window
			if (ay !=0)
			{

				// a vector perpendicular to the side
				bx = 1 ;
				by = -ax/ay*bx ;
				
				// normalisation
				bx = bx / (sqrt ( bx*bx + by*by ) );				
				by = by / (sqrt ( bx*bx + by*by ) );
			}			
			else
			{
				bx = 0 ;
				by = 1 ;
			}

			if (OrientationX == 0 && OrientationY ==0 ){
	
				if (CentreX != -1) {
					// a vector pointing form the center to the first point of the shape
					cx = Polygonex[o][l] + centerx[l] - (CentreX+0.5)*NBPnf ; 
					cy = Polygoney[o][l] + centery[l] - (CentreY+0.5)*NBPnf ;
					
				} else {
					// a vector pointing form the centre of the polygone to the 
					//point of the shape
					cx = Polygonex[o][l] ; 
					cy = Polygoney[o][l] ;
				}
				// inverse the vector to make it pointing to the outside of the 
				// calculation window
				if ( cx*bx + by*cy < 0 )
				{
					bx = -bx ;
					by = -by ;
				}
			} else {
				// Make sure that the vector along the side is pointing in the same
				// direction than the Orientation vector:
				if (OrientationX*bx + OrientationY*by < 0 ) {
					bx *= -1 ;
					by *= -1 ;
				}
			}



			// number of points that will be calculated along the side
			nbp = abs2( ax );
			if ( abs2( ay) > nbp) 
				nbp = abs2( ay ) ;

			// for each point along the side: 
			for (k=0 ; k < nbp; k++ )
			{
				i = centerx[l] + Polygonex[o][l] + (int)((double)k/nbp*ax) ;
				j = centery[l] + Polygoney[o][l] + (int)((double)k/nbp*ay) ;

				if (i == m && j == n)
				{
					 Nx(n,m) = bx ;
					 Ny(n,m) = by ;
					 return 1 ;
				}
				else
				{
					Nx(n,m) = Nx(n,m) + bx/ abs( double(pow(i-m,2) + pow(j-n,2)) );                       
					Ny(n,m) = Ny(n,m) + by/ abs( double(pow(i-m,2) + pow(j-n,2)) );   
				} 
			}

		}

	}

	// weigt the 4 points :
	Nx(n,m) = Nx(n,m) + Nx(j1,i1)/ abs( double(pow(i1-m,2) + pow(j1-n,2)) ); 
	Ny(n,m) = Ny(n,m) + Ny(j1,i1)/ abs( double(pow(i1-m,2) + pow(j1-n,2)) ); 

	Nx(n,m) = Nx(n,m) + Nx(j2,i2)/ abs( double(pow(i2-m,2) + pow(j2-n,2)) );                       
	Ny(n,m) = Ny(n,m) + Ny(j2,i2)/ abs( double(pow(i2-m,2) + pow(j2-n,2)) );  

	Nx(n,m) = Nx(n,m) + Nx(j3,i3)/ abs( double(pow(i3-m,2) + pow(j3-n,2)) );                       
	Ny(n,m) = Ny(n,m) + Ny(j3,i3)/ abs( double(pow(i3-m,2) + pow(j3-n,2)) );  

	Nx(n,m) = Nx(n,m) + Nx(j4,i4)/ abs( double(pow(i4-m,2) + pow(j4-n,2)) );                       
	Ny(n,m) = Ny(n,m) + Ny(j4,i4)/ abs( double(pow(i4-m,2) + pow(j4-n,2)) ); 


	// normalize the vectors :
	norm = sqrt( pow(Nx(n,m).real(),2) + pow(Ny(n,m).real(),2) ) ;
	if (norm !=0)	
	{
		Nx(n,m)= Nx(n,m) /norm ;                       
		Ny(n,m) = Ny(n,m)/norm ;
	}
	return 1;
}




// compute the normal field and store to Nx and Ny matrices
void draw::normal_field_generation(void)
{
	int dx,dy ;
	int R[NB_LAYERS_MAX],Cx[NB_LAYERS_MAX],Cy[NB_LAYERS_MAX] ;
	int Polygonex[NB_EDGES_MAX][NB_LAYERS_MAX],Polygoney[NB_EDGES_MAX][NB_LAYERS_MAX] ;
	int i,j,k,l,m,n ;

	double alpha, nbp ;
	double norm ;
	double ax, ay;
	int xmin, xmax, ymin, ymax ;
	FILE *foutx,*fouty ;

	db_matrix O(NBPnf, NBPnf);
	Nx = O;
	Ny = O ;
	

 	//create the field at the boundaries of the windows
	// and adjust their direction accordingly to OrientationX and OrientationY:
	int sign=1 ;
    for (k=1; k<NBPnf ;k++)
	{
	    Nx(0,k)= 0.0;                       
        Ny(0,k)= -1.0*sign ; 
		if (OrientationY*Ny(0,k).real() < 0)
			Ny(0,k) *=-1*sign ;
		
        Nx(NBPnf-1,k)= 0.0 ;                       
        Ny(NBPnf-1,k)= 1.0 ;
		if (OrientationY*Ny(NBPnf-1,k).real() < 0)
			Ny(NBPnf-1,k) *=-1 ;
    }
	sign=1 ;
    for (k=1; k<NBPnf ;k++)
	{

	    Nx(k,0)= -1.0 ;                      
        Ny(k,0)= 0 ;
		if (OrientationX*Nx(k,0).real() < 0)
			Nx(k,0) *=-1 ;
		
        Nx(k,NBPnf-1)= 1.0*sign;                       
        Ny(k,NBPnf-1)= 0 ;
		if (OrientationX*Nx(k,NBPnf-1).real() < 0)
			Nx(k,NBPnf-1) *=-1*sign ;
    }
    Nx(0,0)= -1.0;                       
    Ny(0,0)= -1.0 ; 
	if (OrientationX*Nx(0,0).real() + OrientationY*Ny(0,0).real()  < 0 ) {
		Nx(0,0) *= -1 ;
		Ny(0,0) *= -1 ;
	}

    Nx(NBPnf-1,NBPnf-1)= 1.0 ;                       
    Ny(NBPnf-1,NBPnf-1)= 1.0; 
	if (OrientationX*Nx(NBPnf-1,NBPnf-1).real() + OrientationY*Ny(NBPnf-1,NBPnf-1).real()  < 0 ) {
		Nx(NBPnf-1,NBPnf-1) *= -1 ;
		Ny(NBPnf-1,NBPnf-1) *= -1 ;
	}
    
    Nx(0,NBPnf-1)= 1.0 ;                       
    Ny(0,NBPnf-1)= -1.0;
	if (OrientationX*Nx(0,NBPnf-1).real() + OrientationY*Ny(0,NBPnf-1).real()  < 0 ) {
		Nx(0,NBPnf-1) *= -1 ;
		Ny(0,NBPnf-1) *= -1 ;
	}
    
    Nx(NBPnf-1,0)= -1.0 ;                       
    Ny(NBPnf-1,0)= 1.0; 
	if (OrientationX*Nx(NBPnf-1,0).real() + OrientationY*Ny(NBPnf-1,0).real()  < 0 ) {
		Nx(NBPnf-1,0) *= -1 ;
		Ny(NBPnf-1,0) *= -1 ;
	}
 
    // set the value at the middle of the windows
    Nx(NBPnf/2,NBPnf/2)= OrientationX;                       
    Ny(NBPnf/2,NBPnf/2)= OrientationY;
    
    // set the value of the radius and circle center in term of pixel:
	for ( i = 0 ; i < nb_layers ; i++)
	{
		Cx[i] = origineX[i] * NBPnf + NBPnf/2 ;
		Cy[i] = origineY[i] * NBPnf + NBPnf/2 ;
		for (j=0 ; j < nb_edges[i] ; j++)
		{
			Polygonex[j][i] = (int)(x_poly[j][i]* (double)(NBPnf)) ;
			Polygoney[j][i] = (int)(y_poly[j][i]* (double)(NBPnf)) ;

		}
   	}
    //interpolate the field
    dx = NBPnf/2 ;
    dy = NBPnf/2 ;

    while (dx > 2 && dy > 2)
	{
        // first rectangle
        //initialize indexes
        i=0;
        k = dx ;
        m = dx/2;
        while ( k < NBPnf)
		{
            j = 0;
            l = dy ;
            n = dy/2 ;
            while (l < NBPnf)
			{
                interpolateNF(Cx,Cy, Polygonex,Polygoney,
                	i,j, i,l, k,j, k,l, m,n) ;
                //shift 
		        j = l ;
		        l = l + dy ;
                n = n + dy ;
            } 
            // shift:
		    i = k ;
		    k = k + dx ;
            m = m + dx ;
        }

        //tilted rectangle

        //initialize indexes
        i=0;
		k = 0 ;
        while ( i < NBPnf-2*dx)
		{
            j = 0;
			l = 0 ;
            while (j < NBPnf-2*dy)
			{
                
                interpolateNF(Cx,Cy, Polygonex,Polygoney,
                	k+dx/2,l+dy/2, i,l+dy, k+dx,l+dy, k+dx/2,l+dy+dy/2, k+dx/2,l+dy) ;
                interpolateNF(Cx,Cy, Polygonex,Polygoney, 
                	k+dx/2,l+dy/2, k+dx,j, k+dx,l+dy, k+dx/2+dx,l+dy/2, k+dx,l+dy/2) ;
                //shift 
                j = j + dy ;
				l = l + dy ;
            }
            // shift:
            i = i + dx ;
			k = k + dx ;
        }
        
        // refine the interpolation:
        dx = dx/2 ;
        dy = dy/2 ;
         
    }

    // for all remaining zero vectors comming from discretization,
    // interpolation is performed
    for (i=1; i < NBPnf-1 ; i++)
	{
		for (j=1; j < NBPnf-1 ; j++)
		{

            if (Nx(j,i).real() == 0 && Ny(j,i).real() == 0 )
			{

				interpolateNF(Cx,Cy ,Polygonex,Polygoney,
					i-1,j, i+1,j, i,j+1, i,j-1, i,j) ;
			
                // normalize
                norm = sqrt(Nx(j,i).real()*Nx(j,i).real() + 
                		Ny(j,i).real()*Ny(j,i).real() ) ;
                if (norm != 0 )
				{
                     Nx(j,i) = Nx(j,i) / norm ;
                     Ny(j,i) = Ny(j,i) / norm ;
                 }
            } 
        } 
    } 
    
	// norm the vectors and save
    for (j=0; j < NBPnf ; j++)
	{
		for (i=0; i < NBPnf ; i++)
		{

            norm = sqrt(Nx(j,i).real()*Nx(j,i).real() + 
					Ny(j,i).real()*Ny(j,i).real() ) ;
            if (norm != 0 )
			{
                Nx(j,i) = Nx(j,i) / norm ;
                Ny(j,i) = Ny(j,i) / norm ;
            }

        }
    }
	
}


// compute the refractive index for the draw shapes and store it in ind matrix
void draw::polygone(void)


{
	int i, j, k,l ;
	int flag;
	double trianglex[3], triangley[3] ;

	db_matrix O(NBP, NBP);
	ind = O;
	//return ;

	for (j = -NBP/2 ; j < NBP/2 ; j++)
	{
		for (i = -NBP/2 ; i < NBP/2 ; i++)
		{
			flag = nb_layers ;
			for (l=nb_layers-1 ; l >=0 ; l--)		// parcours les couches de la couche externe vers la couche interne
			{
//				if ( (i-origineX[l]*NBP)/(r_a[l]*NBP)*(i-origineX[l]*NBP)/(r_a[l]*NBP)+
//					(j-origineY[l]*NBP)/(r_b[l]*NBP)*(j-origineY[l]*NBP)/(r_b[l]*NBP) <= 1 )
//				{
					// au coeur
					// verifie si on est dans le polygone ou non:
					if (contained_in_polygone(l, i, j)) {
						flag = l ;
					}

//				}
			}
			if (flag ==nb_layers){		// on n'est dans aucune des couches:
				ind(NBP/2+j,NBP/2+i)=complex<double>(n_ext,k_ext) ;
			} else {			// apartient  au polygone
				ind(NBP/2+j,NBP/2+i)=complex<double>(n_int[flag],k_int[flag]) ;
			}
		}
	}
}


/** Write on the given file the refractive index or the normal field.
    @param out is the matrix to write on the file
    @param header is a one line string to be written at the begining of the file
    	as a header.
    @param f the file where the field will be written.

*/
void draw::writeoutputfield(db_matrix &out, char* header,
	bool optiwaveformat, FILE *f)
{
    double x;
    double y;

    double dx = sizeX/(out.getNcol()-1) ;
    double dy = sizeY/(out.getNrow()-1);

    if(out.isEmpty()) {
         throw parsefile_commandError("The out matrix shouldn't be empty! "
            "Programming error!");
    }

	// write header
	fprintf(f, "%s\n", header );
	if(optiwaveformat) {
		// write number of points in x and y directions
        fprintf(f, "%d %d\n", out.getNcol(), out.getNrow() );
        // write the size in x and y directions
        fprintf(f, "%le %le %le %le\n", 0.0, sizeX, 0.0, sizeY );
	}

    for(int i=0; i<out.getNrow(); ++i) {
        for(int j=0; j<out.getNcol(); ++j) {
            // calculate the current position (we are in curvilinear
            // coordinates).
            x = j*dx-(out.getNcol()-1)/2.0*dx;
            y = i*dy-(out.getNrow()-1)/2.0*dy;

            // Output of the fields.
            if(optiwaveformat) {
                fprintf(f, "%le, %le \n", out(i,j).real(),out(i,j).imag());
            } else {
                fprintf(f, "%le %le %le %le \n",
                    x, y, out(i,j).real(),out(i,j).imag());
            }
        }
        if(!optiwaveformat) {
        	fprintf(f, "\n");
        }
    }
}

/** Write on the given file the refractive index or the normal field.
    @param out is the matrix to write on the file
    @param header is a one line string to be written at the begining of the file
    	as a header.
    @param filename is the path+filename of the file to write.
*/
void draw::writeoutputfield(db_matrix &out, char* header,
	bool optiwaveformat, char * filename)
{
		
		FILE * fout = fopen(filename,"w");
		if (fout == NULL)
		{
			printf("draw: not able ot open file: %s\n",filename);
	        throw parsefile_commandError("Error in writting file.");
		}
		writeoutputfield(out, header, true, fout) ;
		fclose(fout);
}

/*  TREAT_INPUT_POINT
	function that is called to compute x_poly, y_poly, origineX, origineY, r_a, r_b
	based on point directly input through x_poly, y_poly
	parameters:
		layerID: index of the layer where we want to process data 
		SortAngle: if true, will sort point with increase angle.
*/
void draw::treat_input_point(int layerID, bool SortAngle)
{
	int i ;
	bool change;
	double angle[NB_EDGES_MAX][NB_LAYERS_MAX] ;
	double R ;
	double alpha ;
	
	// compute the center through the average of point coordinates
	origineX[layerID] = 0 ;
	origineY[layerID] = 0 ;
	for (i = 0 ; i < nb_edges[layerID] ; i++) {
		origineX[layerID] += x_poly[i][layerID]/nb_edges[layerID] ;
		origineY[layerID] += y_poly[i][layerID]/nb_edges[layerID] ;
	}

	r_a[layerID] = -1 ;
	for (i = 0 ; i < nb_edges[layerID] ; i++) {
		// refers polygon points to shape origin:
		x_poly[i][layerID] = x_poly[i][layerID] - origineX[layerID] ;
		y_poly[i][layerID] = y_poly[i][layerID] - origineY[layerID] ;

		R = sqrt( pow(x_poly[i][layerID],2) + pow(y_poly[i][layerID],2) ) ;

		// compute r_a as the smallest distance from center
		if (r_a[layerID] < R  )
			r_a[layerID] = R ;
		
		// compute angle with respect of x axis:
		angle[i][layerID] = acos(x_poly[i][layerID]/R) ;
		// if y is negative, the angle is negative :
		if ( y_poly[i][layerID] < 0) 
			angle[i][layerID] = 2*PI - angle[i][layerID] ;			
	}
	r_b[layerID] = r_a[layerID] ;

	if (SortAngle) {

		// sort points such as angle is increasing:
		change=true ;
		while (change) {
			change=false ;
			for (i = 0 ; i < nb_edges[layerID]-1 ; i++) {
				if (angle[i][layerID] > angle[i+1][layerID]) {
					// inverse both values
					alpha = angle[i][layerID] ;
					angle[i][layerID] = angle[i+1][layerID] ;
					angle[i+1][layerID] = alpha;

					R = x_poly[i][layerID] ;
					x_poly[i][layerID] = x_poly[i+1][layerID] ;
					x_poly[i+1][layerID] = R ;

					R = y_poly[i][layerID] ;
					y_poly[i][layerID] = y_poly[i+1][layerID] ;
					y_poly[i+1][layerID] = R ;
					change=true ;
				}
			}
		}
	} else {
		// actually check if the points are ordered in anti-clockwise direction
		// if not, change their order:
		int clockwise = 0 ;
		int anticlockwise = 0 ;
		for (i = 1 ; i < nb_edges[layerID] ; i++) {
			if (angle[i][layerID] > angle[i-1][layerID])
				anticlockwise++ ;
			else
				clockwise++ ;
		}
	
		if (clockwise>anticlockwise) {
			// reversing the order to obtain anti-clockwise order (increase of
			// angle):
			for (i = 0 ; i < nb_edges[layerID]/2 ; i++) {
				printf("%d\n",i);
				// first and last values
				alpha = angle[i][layerID] ;
				angle[i][layerID] = angle[nb_edges[layerID]-1-i][layerID] ;
				angle[nb_edges[layerID]-1-i][layerID] = alpha;

				R = x_poly[i][layerID] ;
				x_poly[i][layerID] = x_poly[nb_edges[layerID]-1-i][layerID] ;
				x_poly[nb_edges[layerID]-1-i][layerID]= R ;

				R = y_poly[i][layerID] ;
				y_poly[i][layerID] = y_poly[nb_edges[layerID]-1-i][layerID] ;
				y_poly[nb_edges[layerID]-1-i][layerID] = R ;
			}

		} 
	}
	
}

/*  TREAT_INPUT_RECTANGLE
	function that is called to compute x_poly, y_poly, r_a, r_b
	based on rectangle directly input through wx, wy, origineX, origineY
	parameters:
		layerID: index of the layer where we want to process data 
		SortAngle: if true, will sort point with increase angle.
*/
void draw::treat_input_rectangle(int layerID)
{
	int i ;
	bool change;
	double angle[NB_EDGES_MAX][NB_LAYERS_MAX] ;
	double R ;
	double alpha ;
	
	nb_edges[layerID] = 4 ;
	i=0 ;
	x_poly[i][layerID] = wx[layerID]/2.0 ;
	y_poly[i++][layerID] = wy[layerID]/2.0 ;
	
	x_poly[i][layerID] = -wx[layerID]/2.0 ;
	y_poly[i++][layerID] = wy[layerID]/2.0 ;
	
	x_poly[i][layerID] = -wx[layerID]/2.0 ;
	y_poly[i++][layerID] = -wy[layerID]/2.0 ;
	
	x_poly[i][layerID] = wx[layerID]/2.0 ;
	y_poly[i++][layerID] = -wy[layerID]/2.0 ;
	
	r_a[layerID] =sqrt(2)*wx[layerID]/2.0 ;
	r_b[layerID] =sqrt(2)*wy[layerID]/2.0 ;


}

/*  TREAT_INPUT_POLYGONE_CIRCLE
	function that is called to compute x_poly, y_poly
	based on point directly input through nb_edges, alpha0, r_a and r_b
	parameters:
		layerID: index of the layer where we want to process data 
		alpha0: angle of origin for the point (rotation of the shape).

*/
void draw::treat_input_polygone_circle(int layerID, double alpha0)
{
	int i ;
	double alpha ;
	alpha = 2*PI/nb_edges[layerID] ;
	for (i = 0 ; i < nb_edges[layerID] ; i++)
	{
		// coordonnées des points du polygone  entre -0.5 et 0.5:
		// dessine le polygone centré en zéro
		x_poly[i][layerID] = r_a[layerID]* cos (alpha0+alpha*i) ;
		y_poly[i][layerID] = r_b[layerID]/(sizeY/sizeX) * sin (alpha0+alpha*i) ;
	}

}
int draw::parse_legacy(int argc,char *argv[])
{

	double a=1 ;
	n_ext = 1;
	k_ext = 0 ;
	int i, j, k;

	char filename_nfx[250] ;
	char filename_nfy[250] ;
	char string[350] ;
		
	int nb_edges_max ;
	sizeY=1.0 ;
	double sizeX=1.0 ;
	double R ;

	double alphar[NB_LAYERS_MAX] ;


	bool change=true ;
	bool SpecifyPoint=false ;
	bool SortAngle=false ;
	bool Orientation=false ;
	bool Centre=false ;
	OnlyIndex=false ;
	OrientationX=0;
	OrientationY=0 ;
	CentreX=-1;
	CentreY=-1 ;
	int nb_option=0 ;
	j = 1 ;
	if (SpecifyPoint == true)
		j++ ;

	// lit les options
	while (argc>j && change==true) {
		change=false ;
		if (strncmp(argv[j], "-pr",3)==0){
			SortAngle = true ;
			SpecifyPoint = true ;
			change=true ;
			j++ ;
			nb_option+=1;
			printf("_pr\n");
		}
		if (argc>j) {
			if (strncmp(argv[j], "-p",2)==0) {
				SpecifyPoint = true ;
				change=true ;
				j++ ;
				nb_option++;
				printf("_p\n");
			}
		}
		if (argc>j) {
			if (strncmp(argv[j], "-o",2)==0){
				Orientation = true ;
				change=true ;
				j++ ;
				nb_option++;
			}
		}
		if (argc>j) {
			if (strncmp(argv[j], "-c",2)==0){
				Centre = true ;
				change=true ;
				j++ ;
				nb_option++;
			}
		}
		if (argc>j) {
			if (strncmp(argv[j], "-i",2)==0){
				OnlyIndex = true ;
				change=true ;
				j++ ;
				nb_option++;
			}
		}
	}

	if (argc==1+nb_option  )
	{
		//tout est tapé à la main
		if (SpecifyPoint==false) {
			printf("NOTA BENE: l'option '-p' permet de specifier directement les points du polynome.\n");
			printf("           l'option '-pr' permet de ranger les points du polynome dans le sens trigonometrique. Si vous n'utilisez pas cette fonction. Rangez les vous même dans cet ordre.\n");
			printf("           l'option '-o' permet de choisir l'orientation du champ normal (par defaut, il est radial.\n");
			printf("           l'option '-c' permet de choisir le centre du champ radial (par defaut, il est au milieu de chaque polygone.\n");
			printf("           l'option '-i' permet de calculer seulement l'indice optique (pas le champ normal).\n");
		}

		if (Centre==true) {
			printf("Coordonnée selon x du centre pour le champ radial ?\n");
			if (scanf("%lf", &CentreX) !=1)
				printf("ERROR: input CentreX\n");
			printf("Coordonnée selon y du centre pour le champ radial ?\n");
			if (scanf("%lf", &CentreY)!=1)
				printf("ERROR: input CentreY\n");
		} else if (Orientation==true) {
			printf("composante selon x (nf_x) de l'orientation principale du champ normal (nf_x, nf_y) ?\n");
			if (scanf("%lf", &OrientationX)!=1)
				printf("ERROR: input OrientationX\n");
			printf("composante selon x (nf_x) de l'orientation principale du champ normal (nf_x, nf_y) ?\n");
			if (scanf("%lf", &OrientationY)!=1)
				printf("ERROR: input OrientationY\n");
		}

		printf("Nom du fichier à écrire\n");
		if (scanf("%s", filename_ind)!=1)
			printf("ERROR: input filename_ind\n");
		//printf("%s", filename_ind);


		printf("taille de la fenetre en y (dimension réduite sizeY/sizeX) ?\n");
		if (scanf("%lf", &sizeY)!=1)
			printf("ERROR: input sizeY\n");
		printf("nombre de couches ?\n");
		if (scanf("%d", &nb_layers)!=1)
			printf("ERROR: input nb_layers\n");

		printf("\ndéfinition des couches, la couche 1 étant la couche centrale !\n");

		for (i = 0 ; i < nb_layers ; i++)
		{
			printf("\nCOUCHE%d\n",i+1);
			if (SpecifyPoint==false) {
				printf("rapport (r/a = rayon / taille de la maille x )\n");
				if (scanf("%lf", &r_a[i])!=1)
					printf("ERROR: input r/a\n");	
				r_b[i] = r_a[i]/sizeY;
			}
			printf("nombre de cotes \n");
			if (scanf("%d", &nb_edges[i])!=1)
				printf("ERROR: input nb_edges\n");
			if (SpecifyPoint==false) {
				printf("origine x (entre 0 et 1) \n");
				if (scanf("%lf", &origineX[i])!=1)
					printf("ERROR: input origineX\n");
				printf("origine y (entre 0 et %g) \n", sizeY);
				if (scanf("%lf", &origineY[i])!=1)
					printf("ERROR: input origineY\n");
				printf("angle à l'origine (deg)\n");
				if (scanf("%lf", &alpha0[i])!=1)
					printf("ERROR: input alpha0\n");	
				alpha0[i] = alpha0[i]/180*PI ;	//rad
			} else {
				for (k = 0 ; k < nb_edges[i] ; k++) {
					printf("x du point %d (entre 0 et 1) \n", k+1);
					if (scanf("%lf", &x_poly[k][i])!=1)
						printf("ERROR: input x\n");
					printf("y du point %d (entre 0 et 1) \n", k+1);
					if (scanf("%lf", &y_poly[k][i])!=1)
						printf("ERROR: input y\n");
				}
			}
			printf("indice réel interne\n");
			if (scanf("%lf", &n_int[i])!=1)
				printf("ERROR: input n_int\n");	
			printf("indice imaginaire interne\n");
			if (scanf("%lf", &k_int[i])!=1)
				printf("ERROR: input k_int\n");		
		}

		printf("indice réel externe\n");
		if (scanf("%lf", &n_ext)!=1)
				printf("ERROR: input n_ext\n");	
		printf("indice imaginaire externe\n");
		if (scanf("%lf", &k_ext)!=1)
				printf("ERROR: input k_ext\n");
	}
	else
	{
		if (Centre==true) {
			CentreX=atof(argv[j++]);
			CentreY=atof(argv[j++]);
		} else if (Orientation==true) {
			OrientationX=atof(argv[j++]);
			OrientationY=atof(argv[j++]);
		}

		sscanf(argv[j++],"%s", filename_ind);

		printf("\n fichier %s  ",filename_ind);

		sizeY = atof(argv[j++]) ;
		nb_layers=atoi(argv[j++]);
		for (i = 0 ; i < nb_layers ; i++)
		{
			if (SpecifyPoint==false) {
				r_a[i]=atof(argv[j++]);
				r_b[i] = r_a[i];
			}
			nb_edges[i]=atoi(argv[j++]);
			if (SpecifyPoint==false) {
				origineX[i]=atof(argv[j++]);
				origineY[i]=atof(argv[j++]);
				alpha0[i] = atof(argv[j++]) ;
				alpha0[i] = alpha0[i]/180*PI ;	//rad
			} else {
				for (k = 0 ; k < nb_edges[i] ; k++) {
					x_poly[k][i] = atof(argv[j++]);
					y_poly[k][i] = atof(argv[j++]);
				}
			}
			n_int[i]=atof(argv[j++]);
			k_int[i]=atof(argv[j++]);

		}
		n_ext=atof(argv[j++]);
		k_ext=atof(argv[j++]);

	}

	// apply sizeY
	for (i = 0 ; i < nb_layers ; i++) {
		origineY[i] = origineY[i]/sizeY;
	}
	
	// convertie origine de -0.5 à 0.5 :
	for (j=0 ; j < nb_layers ; j++)
	{

		if (SpecifyPoint==false) {
			shape_type[j] = 'c' ;
			origineX[j] = (origineX[j]-0.5) ;
			origineY[j] = (origineY[j]-0.5) ;
		} else {
			shape_type[j] = 'p' ;
			for (i = 0 ; i < nb_edges[j] ; i++) {
				// convertie les points du polygone entre -0.5 et 0.5
				x_poly[i][j] = x_poly[i][j]/sizeX -0.5 ;
				y_poly[i][j] = y_poly[i][j]/sizeY -0.5 ;
			}
			treat_input_point(j, SortAngle);

//		printf("origine = %g %g \n",origineX[j],origineY[j]);
		printf("{ %d-cote r/a = %e ;x=%g y=%g ; n_int=%g ; k_int=%g } ",nb_edges[j],r_a[j],origineX[j],origineY[j],n_int[j],k_int[j]);
		}
	}
			

	printf(" ; n_ext=%g ; k_ext=%g \n",n_ext,k_ext);

	if (SpecifyPoint==false) {
		// genere les points du polygone :
		for (j=0 ; j < nb_layers ; j++)
		{
			treat_input_polygone_circle(j, alpha0[j]) ;
		}
	}
	if (Centre==true) {
		CentreX = CentreX - 0.5 ;
		CentreY =CentreY/sizeY -0.5 ;
		printf("Le centre du champ normal est pris à (%g, %g)\n", CentreX, CentreY);
	}
	// normalise le vecteur d'orientation principal:
	if (Orientation==true) {
		R = sqrt(OrientationX*OrientationX + OrientationY*OrientationY) ;
		OrientationX /= R ;
		OrientationY /= R ;
		printf("L'orientation du champ normal est pris selon (%g, %g)\n",OrientationX, OrientationY);
	}
	// calcul le champ d'indice :

	
	printf("creation du fichier d'indice\n");
	polygone() ;
	printf("Writing in %s\n",filename_ind);
	
	generate_command(string);
		
	writeoutputfield(ind, string, true, filename_ind) ;
		
	// champ normal :
	if (!OnlyIndex) {
		printf("creation du fichier de champ normal\n");
		sprintf(filename_nfx,"%s_nvf_x", filename_ind);
		sprintf(filename_nfy,"%s_nvf_y", filename_ind);

		normal_field_generation() ;
		writeoutputfield(Nx, string, true, filename_nfx) ;
		writeoutputfield(Ny, string, true, filename_nfy) ;

	}

	printf("Executed command: %s\n", string);
	

	return 0;
}
/*
	parse_input
	function that is either parsing value from argv or asking from user
	in case of interactive
	
	parameters:
	argc: number of arguments
	argv: vector with text arguments
	readCount: pointer to number of already read from argv
	interactive = true, will ask user if needed
	ignoreOptional = true, skip all parameters inside brakets [...]
	text is the text of the request to user
	out is a pointer where input value will be returned. out can be char*, 
		double or int (i.e. 3 prototpe functions)
*/
void draw::parse_input(int argc,char *argv[], int &readCount, bool interactive,
		bool ignoreOptional, const char* text, char *out)  {
	
	char string[100] ;
	
	if (readCount>argc-1) {
		printf("%s? ", text);
		if (interactive) {
			if (scanf("%s", out)!=1){
				sprintf(string, "ERROR: getting %s\n", text) ;
				throw parsefile_commandError(string);
			 }
		} else {
			sprintf(string, "draw: parameter '%s' missing\n", text) ;
			throw parsefile_commandError(string);
		}
	} else {

		if (argv[readCount][0] == '[') {
			if (ignoreOptional){
				parse_skip_optional(argc,argv,readCount);
				strcpy(out, argv[readCount++]);
			}else
				strcpy(out, &argv[readCount++][1]);
		} else
			strcpy(out, argv[readCount++]);
			
		if (out[strlen(out)-1] == ']')
			out[strlen(out)-1] = '\0';
			
		if (interactive) {
			printf("%s = %s\n",text, out); 
		 }
	}
}

void draw::parse_input(int argc,char *argv[], int &readCount, bool interactive,
		bool ignoreOptional, const char* text, int &out)  {

	char string[100] ;
	if (readCount>argc-1) {
		printf("%s? ", text);
		if (interactive) {
			if (scanf("%d", &out)!=1){
				sprintf(string, "ERROR: getting %s\n", text) ;
				throw parsefile_commandError(string);
			 }
		} else {
			sprintf(string, "draw: parameter '%s' missing\n", text) ;
			throw parsefile_commandError(string);
		}
	} else {
		if (argv[readCount][0] == '[') {
			if (ignoreOptional){
				parse_skip_optional(argc,argv,readCount);
				out=atoi(argv[readCount++]);
			}else
				out=atoi(&argv[readCount++][1]);
		} else
			out=atoi(argv[readCount++]);
			
		if (interactive) {
			printf("%s = %d\n",text, out); 
		 }
	}
}

void draw::parse_input(int argc,char *argv[], int &readCount, bool interactive,
		bool ignoreOptional, const char* text, double &out)  {

	char string[100] ;
	if (readCount>argc-1) {
		printf("%s? ", text);
		if (interactive) {
			if (scanf("%lf", &out) !=1){
				sprintf(string, "ERROR: getting %s\n", text) ;
				throw parsefile_commandError(string);
			 }
		} else {
			sprintf(string, "draw: parameter '%s' missing\n", text) ;
			throw parsefile_commandError(string);
		}
	} else {
		if (argv[readCount][0] == '[') {
			if (ignoreOptional){
				parse_skip_optional(argc,argv,readCount);
				out=atof(argv[readCount++]);
			} else
				out=atof(&argv[readCount++][1]);
		} else
			out=atof(argv[readCount++]);
		if (interactive) {
			printf("%s = %g\n",text, out); 
		 }
	}
}
/*
	parse_skip_optional
	function that will continue reading parameters when parameters are whithin
	square braket [...]
	
	parameters:
	argc: number of arguments
	argv: vector with text arguments
	readCount: pointer to number of already read from argv
*/
void draw::parse_skip_optional(int argc,char *argv[], int &readCount){

	if (argv[readCount][0] == '['){
		while (argv[readCount][strlen(argv[readCount])-1] != ']' 
			and readCount < argc) {
			readCount++;
		}

		if (argv[readCount][strlen(argv[readCount])-1] == ']')
			readCount++;
		else
	         throw parsefile_commandError("draw command without closing braket ]!");
	}
}
/*
	PARSE_STANDALONE the following arguments. See help for more details
	
		draw filename sizeX sizeY
		nb_layer 
  		shape_type parameters 
  		shape_type parameters 
  		...  				
  		n_ext k_ext 
		OPTIONAL_PARAMETERS		

	NB: save results directly to file
*/

int draw::parse_standalone(int argc,char *argv[]) 
{
	int readCount=1 ;
	bool interactive = true; 
	bool skipOptional=false ;
	char filename_nfx[250] ;
	char filename_nfy[250] ;
	char header[300] ;
		
	sizeX = 1 ;
		
	// filename
	parse_input(argc,argv, readCount, interactive, skipOptional, "filename", filename_ind);
	sprintf(filename_nfx,"%s_nvf_x", filename_ind);
	sprintf(filename_nfy,"%s_nvf_y", filename_ind);

	// sizeX:
	parse_input(argc,argv, readCount, interactive, skipOptional, 
		"size of computation window in x", sizeX);

	// sizeY:
	parse_input(argc,argv, readCount, interactive, skipOptional,
		"size of computation window in y", sizeY);
	
	// parse the rest of the parameters
	parse(argc-readCount+1,&argv[readCount-1], interactive);
	
	
	// write the index and normal field to files

	printf("Writing in %s\n",filename_ind);
	
	generate_command(header);
	writeoutputfield(ind, header, true, filename_ind) ;
		
	// champ normal :
	if (!OnlyIndex) {
		sprintf(filename_nfx,"%s_nvf_x", filename_ind);
		sprintf(filename_nfy,"%s_nvf_y", filename_ind);
		
		printf("writing normal field files: %s and %s\n", filename_nfx, filename_nfy);		
		writeoutputfield(Nx, header, true, filename_nfx) ;
		writeoutputfield(Ny, header, true, filename_nfy) ;

	}
	return 0 ;
	
}
    
/*
	PARSE the following arguments. See help for more details
	
		draw nb_layer 
  		shape_type parameters 
  		shape_type parameters 
  		...  				
  		n_ext k_ext 
		OPTIONAL_PARAMETERS		



*/
int draw::parse( int argc,char *argv[], bool interactive)
{
	int i, readCount=1 , k;
	char string[350] ;

	bool skipOptional = true ;
	// centre =-1 = in the center of the shape
	CentreX=-1;
	CentreY=-1 ;
	
	// no orientation set:
	OrientationX=0;
	OrientationY=0 ;
	
	
	OnlyIndex = false;
	bool orderPoints = false;
	

	// first argument is nb_layer
	parse_input(argc,argv, readCount, interactive, skipOptional,
		 "number of layers", nb_layers);
	
	for (i = 0 ; i < nb_layers ; i++) {
		// parse shape_type
		parse_input(argc,argv, readCount, interactive, skipOptional,
			"shape type (c, e, r or p)", string);
		shape_type[i] = string[0] ;
	
		//n_int
		parse_input(argc,argv, readCount, interactive, skipOptional,
			"Internal real index", n_int[i]);

		//k_int
		parse_input(argc,argv, readCount, interactive, skipOptional,
			"Internal imaginary index", k_int[i]);
			
		if (shape_type[i] == 'c' or shape_type[i] == 'e') {
			//rx
			sprintf(string, "radius along x (between 0 and %g)",sizeX/2.0); 
			parse_input(argc,argv, readCount, interactive, skipOptional,
				string, r_a[i]);

			// ry
			if (shape_type[i] == 'e') {
				sprintf(string, "radius along y (between 0 and %g)",sizeY/2.0);
				parse_input(argc,argv, readCount, interactive, skipOptional,
					string, r_b[i]);	
			} else
				r_b[i] = r_a[i];
			
			// nb_edge:
			parse_input(argc,argv, readCount, interactive, skipOptional, 
				"Number of edges", nb_edges[i]);

			
		} else if (shape_type[i] == 'r') {
			// wx
			sprintf(string, "width of the rectangle along x (between 0 and %g)"
				,sizeX);
			parse_input(argc,argv, readCount, interactive, skipOptional,
				string, wx[i]);
			
			// wy
			sprintf(string, "width of the rectangle along y (between 0 and %g)"
				,sizeY);
			parse_input(argc,argv, readCount, interactive, skipOptional,
				string, wy[i]);
		} else if (shape_type[i] == 'p') {
			// nb_edge:
			parse_input(argc,argv, readCount, interactive, skipOptional,
				"Number of edges", nb_edges[i]);
			for (k = 0 ; k < nb_edges[i] ; k++) {
				
				//x_k
				sprintf(string, "x for point %d (between -%g and %g)", k,
					sizeX/2.0, sizeX/2.0); 
				parse_input(argc,argv, readCount, interactive, skipOptional, 
					string, x_poly[k][i]);
				
				// y_k
				sprintf(string, "y for point %d (between -%g and %g)", k,
					sizeY/2.0, sizeY/2.0); 
				parse_input(argc,argv, readCount, interactive, skipOptional,
					string, y_poly[k][i]);

			}
				
		} else {
			sprintf(string, "ERROR parsing shape_type = '%c' and should be {c|e|r|p}\n"
				,shape_type[i]);
			throw parsefile_commandError(string);
		}
		
		if (shape_type[i] != 'p') {
			// px
			sprintf(string, "Origine X (between -%g and %g)",sizeX/2.0, sizeX/2.0); 
			parse_input(argc,argv, readCount, interactive, skipOptional, 
				string, origineX[i]);
			
			//py
			sprintf(string, "Origine Y (between -%g and %g)",sizeY/2.0, sizeY/2.0); 
			parse_input(argc,argv, readCount, interactive, skipOptional,
				string, origineY[i]);
		}
		
		// shift_angle
		if (shape_type[i] == 'c' or shape_type[i] == 'e') {
			parse_input(argc,argv, readCount, interactive, skipOptional,
				"Shift angle (deg)", alpha0[i]);
			alpha0[i] +=360.0/(2.0*nb_edges[i]) ;
			alpha0[i] = alpha0[i]/180*PI ;	//rad
		}
										
	}
	// parsing n_ext k_ext
	parse_input(argc,argv, readCount, interactive, skipOptional,
		 "External real index", n_ext);
	parse_input(argc,argv, readCount, interactive, skipOptional,
		"External imaginary index", k_ext);

	// parsing optional_parameters
	while (argc>readCount) {

		if (strncmp(argv[readCount], "-order-points",13)==0 ){
			orderPoints = true ;
			readCount++ ;
			if (interactive)
				printf("Points of the shapes will ordered in anticlockwise direction\n");
		} else if (strncmp(argv[readCount], "-o",2)==0 ){
			// OrientationX 
			parse_input(argc,argv, ++readCount, interactive, skipOptional,
				"Normal field orientation X", OrientationX);
			
			//OrientationY
			parse_input(argc,argv, readCount, interactive, skipOptional,
				"Normal field orientation Y", OrientationY);
			if (interactive)
				printf("orientation set to (%g, %g)\n",OrientationX, OrientationY);		
		} else if (strncmp(argv[readCount], "-c",2)==0 ){
			// CentreX
			sprintf(string, "Normal field center X (between -%g and %g)",
				sizeX/2.0, sizeX/2.0); 
			parse_input(argc,argv, ++readCount, interactive, skipOptional,
				 string, CentreX);
			
			//CentreY
			sprintf(string, "Normal field center Y (between -%g and %g))",
				sizeY/2.0, sizeY/2.0); 
			parse_input(argc,argv, readCount, interactive, skipOptional,
				string, CentreY);
			if (interactive)
				printf("Normal field center set to (%g, %g)\n",CentreX, CentreY);
		} else if (strncmp(argv[readCount], "-i",2)==0 or 
				strncmp(argv[readCount], "-no-nf",6)==0){
			OnlyIndex = true ;
			readCount++ ;
			if (interactive)
				printf("Only computing refractive index (without normal field )\n");
		} else {
			sprintf(string, "ERROR: unkown %s optional parameter\n", argv[readCount]);
			throw parsefile_commandError(string);
		}

	}	
	
	// reduce variables
	if (CentreX != -1 and CentreY!=-1) {
		CentreX = CentreX/sizeX;
		CentreY = CentreY/sizeY;
	}
	for (i = 0 ; i < nb_layers ; i++) {

		origineX[i] = origineX[i]/sizeX ;
		origineY[i] = origineY[i]/sizeY ;

		// generate poly points
		if (shape_type[i] == 'c' or shape_type[i] == 'e') {
			r_a[i] = r_a[i]/sizeX ;
			r_b[i] = r_b[i]/sizeX ;
			treat_input_polygone_circle(i, alpha0[i]);
		} else if (shape_type[i] == 'r') {
			wx[i] = wx[i]/sizeX ;
			wy[i] = wy[i]/sizeY ;
			treat_input_rectangle(i) ;
		} else if (shape_type[i] == 'p') {
			for (k = 0 ; k < nb_edges[i] ; k++) {
				x_poly[k][i] = x_poly[k][i] /sizeX ;
				y_poly[k][i] = y_poly[k][i] /sizeY ;
			}
			treat_input_point(i, orderPoints) ;
		}else {
			sprintf(string, "ERROR parsing shape_type = '%c' and should be {c|e|r|p}\n"
				,shape_type[i]);
			throw parsefile_commandError(string);
		}
	}
	

	generate_command(string);
	printf("Command executed: %s", string);
	
	// generate the index and normal fields	
	polygone() ;
	printf("Refractive index map generated\n");
	
	// Normal field :
	if (!OnlyIndex) {
		normal_field_generation() ;
		printf("Normal field maps generated\n");
	}
	
	return 0;
}

    
/*
	GENERATE COMMAND 
	Generate command that can then be input to the parser to generate again
	the structure. See help for more details
	
		draw nb_layer 
  		shape_type parameters 
  		shape_type parameters 
  		...  				
  		n_ext k_ext 
		OPTIONAL_PARAMETERS		

	output:
	command generated


*/
void draw::generate_command(char *command)
{
	command[0] = '\0';
	
	int i,k ;
	char string[300];
	
	sprintf(string, "draw [%s %g %g] %d", filename_ind, sizeX, sizeY, nb_layers );
	strcat(command, string);
		
	for (i = 0 ; i < nb_layers ; i++) {
		sprintf(string, " %c %g %g", shape_type[i], n_int[i], k_int[i] );
		strcat(command, string);
		if (shape_type[i] == 'c') {
			sprintf(string, " %g %d", r_a[i]*sizeX, nb_edges[i]);
			strcat(command, string);
		}else if (shape_type[i] == 'e') {
			sprintf(string, " %g %g %d", r_a[i]*sizeX, r_b[i]*sizeY, nb_edges[i]);
			strcat(command, string);
		} else if (shape_type[i] == 'r') {
			sprintf(string, " %g %g", wx[i]*sizeX, wy[i]*sizeY);
			strcat(command, string);
		} else if (shape_type[i] == 'p') {
			sprintf(string, " %d", nb_edges[i]);
			strcat(command, string);
			for (k = 0 ; k < nb_edges[i] ; k++) {
				sprintf(string, " %g %g", (x_poly[k][i]+ origineX[i])*sizeX,
						(y_poly[k][i] + origineY[i])*sizeY);
				strcat(command, string);
			}		
		} else {
			sprintf(string, "ERROR programming error: shape_type = '%c' and should be {c|e|r|p}\n"
				,shape_type[i]);
			throw parsefile_commandError(string);
		}
		
		if (shape_type[i] != 'p') {
			sprintf(string, " %g %g", origineX[i]*sizeX, origineY[i]*sizeY);
			strcat(command, string);
		}
		
		// shift_angle
		if (shape_type[i] == 'c' or shape_type[i] == 'e') {
			sprintf(string, " %g", (alpha0[i]/PI - 1.0/nb_edges[i])*180.0);
			strcat(command, string); 
		}	
	}
	sprintf(string, " %g %g ", n_ext, k_ext); 
	strcat(command, string);
	if (OrientationX!=0 or OrientationY!=0) {
		sprintf(string, " -o %g %g ", OrientationX, OrientationY);
		strcat(command, string);
	}
	if ( CentreX!=-1 or CentreY!=-1) {
		sprintf(string, " -c %g %g ", CentreX, CentreY);
		strcat(command, string);
	}
	if (OnlyIndex) {
		strcat(command, " -no-nf");
	}
		
	strcat(command, "\n");
}

void draw::help(void)
{
const char *helptext =
"USAGE AS cRCWA COMMAND		 													\n"														
"draw [...] nb_layer ...														\n"
"  		shape_type parameters ... 												\n"
"  		shape_type parameters ...												\n"
"  		...  																	\n"
"  		n_ext k_ext ...															\n"
"		OPTIONAL_PARAMETERS														\n"	
"																				\n"
"USAGE AS STANDALONE TOOL: \\ 													\n"
"./draw test.ind sizeX sizeY \\													\n"
" 		nb_layer \\ 															\n"
"  		shape_type parameters \\ 												\n"
"  		shape_type parameters \\ 												\n"
"  		... \\ 																	\n"
"  		n_ext k_ext \\ 															\n"
"		OPTIONAL_PARAMETERS														\n"
"																				\n"
"	shape_type = c or e or r or p												\n"
"		- c: polygone generated from CIRCLE: 									\n"
"			all points are taken regularly on a circle							\n"
"			* parameters =  n_int k_int rx nb_edge px py shift_angle			\n"
"		- e: polygone generated with all points are taken regularly on an ELIPSE\n"
"			* parameters = n_int k_int rx ry nb_edge px py shift_angle 			\n"
"		- r: RECTANGLE generated												\n"
"			* parameters = n_int k_int wx wy px py 								\n"
"		- p: shape generated with specifying its POINTS							\n"
"			* parameters = n_int k_int nb_edge x0 y0 x1 y1 ... xN yN			\n"
"																				\n"
"	OPTIONAL_PARAMETERS are:													\n"
"		- '-o' ox oy : set the orientation of the normal field. 				\n"
"					ox and oy are the x- y- direction of the normal field		\n"
"				Default is radial normal fileld									\n"
"		- '-c' cx, cy : set the cordinate of the point where the normal field is\n"
"				pointing to. For example, when the point is inside a polygone,	\n" 
"				normal field will be close to radial with a vortex inside the	\n"
"				shape. When the point is pointing far away from the shape, the 	\n"
"				normal field parallel to each other pointing toward the specified\n"
"				point 															\n"	
"				cx and cy are the x- y- coordinate of the center point			\n"
"				By default it is centered on the polygone						\n"
"		- '-i' or '-no-nf' only generate the refractive index, not the  		\n"
"				normal field. When used inside crwca, previously computed 		\n"
"				nf will be used.												\n" 
"		- '-order-points' when considering input points, order them		  		\n"
"				   																\n" 
"	General parameter description:												\n"
"		* sizeX and sizeY are respectively the size of the computation window	\n"
"		* n_int, k_int are respectively the real and imaginary index of the shape\n"
"		* wx, wy are the x- and y- width of the shape							\n"
"		* rx, ry are the radius of the circle or elipse along x and y directions\n"
"		* nb_edge is the number of edge of the shape							\n"
"		* shift_angle is an added angle to rotate the shape. The origin of angle\n"
" 			is taken as 360/(2*nb_edges) so that the side of the polygone is 	\n"
"			always perpendicular to the x axis.									\n"
"		* px, py, x- and y-position of the center of the shape. 				\n"
"			px=py=0 means that the shape is centered in the middle of the 		\n"
"			computational window.												\n"
"		* x0,y0..xN,yN are points to specify the position of the edges			\n"
"		* nf_x and nf_y are the x and y component of the normal field			\n"
"																				\n"
"																				\n"
"N.B. [...]: parameters inside square brakets are ignored when called from crcwa\n"
"   and read when using command line. It can be advantageously used to input    \n"
"   [test.ind sizeX sizeY] parameters											\n"
"N.B. N+1 layer is below the layer N. Layer 1, is therefore on top of all other	\n"	
" 	layers.																		\n"
"																				\n"
"USAGE IN LEGACY MODE															\n"
"./draw --legacy 																\n"
" or																			\n"
"./draw --legacy test.ind sizeY \\												\n"
" 		nb_layer \\																\n"
"  		rx nb_edge px py angle_shift n_int k_int \\								\n"
"  		rx nb_edge px py angle_shift n_int k_int \\								\n"
"  		n_ext k_ext																\n" ;

printf("%s", helptext);
}
