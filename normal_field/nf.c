#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#define NBP 1024    // multiple de 2
#define NBPnf 256   // multiple de 2

#define PI 3.141592654

/* Normal field calculation based on the work of Götz et al.
    BibTeX reference:


 @article{gotz2008normal,
  title={Normal vector method for the RCWA with automated vector field generation},
  author={G{\"o}tz, P. and Schuster, T. and Frenner, K. and Rafler, S. and Osten, W.},
  journal={Optics express},
  volume={16},
  number={22},
  pages={17295--17301},
  year={2008},
  publisher={Optical Society of America}
 }

*/


typedef enum tag_tc {autre, disc, rect, coeurcoquille} tc;


int interpolateRect (double Nx[NBPnf][NBPnf],double Ny[NBPnf][NBPnf], 
    int centerx, int centery,int width, int height, // Carac du rectangle
    int i1,int j1,int i2, int j2,   // Points 1 et 2
    int i3,int j3,int i4,int j4,    // Points 3 et 4
    int m,int n)                    // Center
{

    int nbp, k, i,j, c;
    double alpha;
    Nx[n][m] = 0;
    Ny[n][m] = 0;
    double norm;
    
    double cx, cy;

    nbp = 2*(width+height); 
    
    for (c=0; c<4; ++c) {
        for (k=0 ; k < nbp; ++k) {

            if(c==0) {
                // Vertical gauche
                i=centerx-width/2;
                j=centery-(height/2)+k*height/nbp;
            } else if(c==1) {
                // Vertical droit
                i=centerx+width/2;
                j=centery-(height/2)+k*height/nbp;
            } else if(c==2) {
                // Horisontal haut
                i=centerx-width/2+k*width/nbp;
                j=centery+height/2;
            } else if(c==3) {
                // Horisontal bas
                i=centerx-width/2+k*width/nbp;
                j=centery-height/2;
            }
            
            if(c==0) {
                cx=-1;
                cy=0;
            } else if(c==1) {
                cx=1;
                cy=0;
            } else if(c==2) {
                cx=0;
                cy=1;
            } else if(c==3) {
                cx=0;
                cy=-1;
            }
            
            if (i == m && j == n) { 
                Nx[n][m]=cx;
                Ny[n][m]=cy;
                return 1 ;
            } else {
                Nx[n][m] = Nx[n][m] + cx/fabs((double)(pow(i-m,2) +
                    pow(j-n,2)));                       
                Ny[n][m] = Ny[n][m] + cy/fabs((double)(pow(i-m,2) + 
                    pow(j-n,2)));     
            } 
        }
    }
    
    Nx[n][m] = Nx[n][m] + Nx[j1][i1]/ fabs((double)(pow(i1-m,2) + pow(j1-n,2))); 
    Ny[n][m] = Ny[n][m] + Ny[j1][i1]/ fabs((double)(pow(i1-m,2) + pow(j1-n,2))); 

    Nx[n][m] = Nx[n][m] + Nx[j2][i2]/ fabs((double)(pow(i2-m,2) + pow(j2-n,2)));                       
    Ny[n][m] = Ny[n][m] + Ny[j2][i2]/ fabs((double)(pow(i2-m,2) + pow(j2-n,2)));  

    Nx[n][m] = Nx[n][m] + Nx[j3][i3]/ fabs((double)(pow(i3-m,2) + pow(j3-n,2)));                       
    Ny[n][m] = Ny[n][m] + Ny[j3][i3]/ fabs((double)(pow(i3-m,2) + pow(j3-n,2)));  

    Nx[n][m] = Nx[n][m] + Nx[j4][i4]/ fabs((double)(pow(i4-m,2) + pow(j4-n,2)));                       
    Ny[n][m] = Ny[n][m] + Ny[j4][i4]/ fabs((double)(pow(i4-m,2) + pow(j4-n,2))); 

    // normalize the vectors :
    norm = sqrt(pow(Nx[n][m],2) + pow(Ny[n][m],2));
    if (norm !=0)   
    {
        Nx[n][m] = Nx[n][m]/norm;                       
        Ny[n][m] = Ny[n][m]/norm;
    }
    return 1;
}

int interpolateDisc (double Nx[NBPnf][NBPnf],double Ny[NBPnf][NBPnf], int centerx, int centery,int radiusx, int radiusy, int i1,int j1,int i2,int j2,int i3,int j3,int i4,int j4,int m,int n)
{

    int nbp, k, i,j ;
    double alpha ;
    Nx[n][m] = 0;
    Ny[n][m] = 0;
    double norm ;
    // weight the point from the radius:
    nbp = 7*radiusx ;   // 2*pi radius
    for (k=0 ; k < nbp; k++ )
    {
        alpha = -PI + (double)k/nbp*2*PI ;

        i =centerx + (int)((double)(radiusx)* cos(alpha));
        j =centery + (int)((double)(radiusy)* sin(alpha));
        if (i == m && j == n)
        {
             Nx[n][m] = cos(alpha) ;
             Ny[n][m] = sin(alpha) ;
             return 1 ;
        }
        else
        {
            Nx[n][m] = Nx[n][m] + cos(alpha)/ fabs((double)(pow(i-m,2) +
                pow(j-n,2)));                       
            Ny[n][m] = Ny[n][m] + sin(alpha)/ fabs((double)(pow(i-m,2) + 
                pow(j-n,2)));     
        } 
        /*Nx[j][i] = cos(alpha) ;
        Ny[j][i] = sin(alpha) ;*/
    }
//printf("%d %d Nx = %g , Ny = %g \n",n,m,Nx[n][m], Ny[n][m] );
    // weigt the 4 points :
    Nx[n][m] = Nx[n][m] + Nx[j1][i1]/ fabs((double)(pow(i1-m,2) + pow(j1-n,2))); 
    Ny[n][m] = Ny[n][m] + Ny[j1][i1]/ fabs((double)(pow(i1-m,2) + pow(j1-n,2))); 

    Nx[n][m] = Nx[n][m] + Nx[j2][i2]/ fabs((double)(pow(i2-m,2) + pow(j2-n,2)));                       
    Ny[n][m] = Ny[n][m] + Ny[j2][i2]/ fabs((double)(pow(i2-m,2) + pow(j2-n,2)));  

    Nx[n][m] = Nx[n][m] + Nx[j3][i3]/ fabs((double)(pow(i3-m,2) + pow(j3-n,2)));                       
    Ny[n][m] = Ny[n][m] + Ny[j3][i3]/ fabs((double)(pow(i3-m,2) + pow(j3-n,2)));  

    Nx[n][m] = Nx[n][m] + Nx[j4][i4]/ fabs((double)(pow(i4-m,2) + pow(j4-n,2)));                       
    Ny[n][m] = Ny[n][m] + Ny[j4][i4]/ fabs((double)(pow(i4-m,2) + pow(j4-n,2))); 

    // normalize the vectors :
    norm = sqrt( pow(Nx[n][m],2) + pow(Ny[n][m],2) ) ;
    if (norm !=0)   
    {
        Nx[n][m] = Nx[n][m]/norm;                       
        Ny[n][m] = Ny[n][m]/norm;
    }
    return 1;
}
int normal_field_generation(const char *filename_x,const char *filename_y,
    tc typeint, double r_a, double r_b, double Tx, double Ty)
{
    int dx,dy;
    int Rx,Ry,Cx,Cy;
    int i,j,k,l,m,n;
    double Nf_x[NBPnf][NBPnf] ,Nf_y[NBPnf][NBPnf];
    double alpha, nbp ;
    double norm ;
    FILE *foutx,*fouty ;
    
    if(typeint==coeurcoquille) typeint=disc;

    // open the files :

    foutx = fopen(filename_x,"w");
    if (foutx == NULL)
    {
        printf("impossible d'ouvrir le fichier : %s\n",filename_x);
        exit(0);
    }

    fouty = fopen(filename_y,"w");
    if (fouty == NULL)
    {
        printf("impossible d'ouvrir le fichier : %s\n",filename_y);
        exit(0);
    }

    // write the header
    fprintf(foutx,"normal vector field x : r/a = %e \n",r_a);
    fprintf(foutx,"%d %d\n",NBPnf, NBPnf);
    fprintf(foutx,"0 %f 0 %f\n", Tx, Ty);

    fprintf(fouty,"normal vector field y : r/a = %e \n",r_a);
    fprintf(fouty,"%d %d\n",NBPnf, NBPnf);
    fprintf(fouty,"0 %f 0 %f\n", Tx, Ty);

    // initialize the array :
    for (i=1; i<NBPnf-1 ;i++)
    {
        for (j=1; j<NBPnf-1 ;j++)
        {
            Nf_x[j][i]= 0.0;                       
            Nf_y[j][i]= 0.0 ; 
        }
    }
    //create the field at the boundaries of the windows
    for (k=1; k<NBPnf ;k++)
    {
        Nf_x[0][k]= 0.0;                       
        Nf_y[0][k]= -1.0 ; 
        
        Nf_x[NBPnf-1][k]= 0.0 ;                       
        Nf_y[NBPnf-1][k]= 1.0 ;
    }
    for (k=1; k<NBPnf ;k++)
    {
        Nf_x[k][0]= -1.0 ;                      
        Nf_y[k][0]= 0 ;
        
        Nf_x[k][NBPnf-1]= 1.0;                       
        Nf_y[k][NBPnf-1]= 0 ;
    }
    Nf_x[0][0]= -1.0;                       
    Nf_y[0][0]= -1.0 ; 
        
    Nf_x[NBPnf-1][NBPnf-1]= 1.0 ;                       
    Nf_y[NBPnf-1][NBPnf-1]= 1.0; 
    
    Nf_x[0][NBPnf-1]= 1.0 ;                       
    Nf_y[0][NBPnf-1]= -1.0;
    
    Nf_x[NBPnf-1][0]= -1.0 ;                       
    Nf_y[NBPnf-1][0]= 1.0;  
    // set the value at the middle of the windows
    Nf_x[NBPnf/2][NBPnf/2]= 0.0;                       
    Nf_y[NBPnf/2][NBPnf/2]= 0.0;
    
    // set the size of the structure in pixels:
    Rx = (int)(r_a* (double)(NBPnf));
    Ry = (int)(r_b* (double)(NBPnf));
    Cx = NBPnf/2 ;
    Cy = NBPnf/2 ;
   
    //interpolate the field
    dx = NBPnf/2 ;
    dy = NBPnf/2 ;

    while (dx > 2 && dy > 2)
    {
        // first rectangle
        //initialize indexes
        i=0;
        k = dx-1 ;
        m = dx/2-1;
        while ( k <= NBPnf)
        {
            j = 0;
            l = dy-1 ;
            n = dy/2-1 ;
            while (l <= NBPnf)
            {
                if(typeint==disc) {
                    interpolateDisc(Nf_x,Nf_y,Cx,Cy,Rx,Rx,i,j,i,l,k,j, k,l,m,n);
                } else {
                    interpolateRect(Nf_x,Nf_y,Cx,Cy,Rx,Ry,i,j,i,l,k,j, k,l,m,n);
                }
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
        k = -1 ;
        while (i <= NBPnf-2*dx)
        {
            j = 0;
            l = -1 ;
            while (j <= NBPnf-2*dy) {
                
                if(typeint==disc) {
                    interpolateDisc(Nf_x,Nf_y, Cx,Cy ,Rx,Rx, k+dx/2,l+dy/2, i,l+dy, k+dx,l+dy, k+dx/2,l+dy+dy/2, k+dx/2,l+dy);
                    interpolateDisc(Nf_x,Nf_y, Cx,Cy, Rx,Rx, k+dx/2,l+dy/2, k+dx,j, k+dx,l+dy, k+dx/2+dx,l+dy/2, k+dx,l+dy/2);
                } else {
                    interpolateRect(Nf_x,Nf_y, Cx,Cy ,Rx,Ry, k+dx/2,l+dy/2, i,l+dy, k+dx,l+dy, k+dx/2,l+dy+dy/2, k+dx/2,l+dy);
                    interpolateRect(Nf_x,Nf_y, Cx,Cy, Rx,Ry, k+dx/2,l+dy/2, k+dx,j, k+dx,l+dy, k+dx/2+dx,l+dy/2, k+dx,l+dy/2);
                }
                
                //shift 
                j = j + dy;
                l = l + dy;
            }
            // shift:
            i = i + dx;
            k = k + dx;
        }
        
        // refine the interpolation:
        dx = dx/2;
        dy = dy/2;
    }

    // for all remaining zero vectors comming from discretization,
    // interpolation is performed
    for (i=1; i < NBPnf-1 ; i++)
    {
        for (j=1; j < NBPnf-1 ; j++)
        {
            if (Nf_x[j][i] == 0 && Nf_y[j][i] == 0 )
            {
                if(typeint==disc) {
                    interpolateDisc(Nf_x,Nf_y, Cx,Cy ,Rx,Rx, i-1,j, i+1,j, i,j+1, i,j-1, i,j) ;
                } else {
                    interpolateRect(Nf_x,Nf_y, Cx,Cy ,Rx,Ry, i-1,j, i+1,j, i,j+1, i,j-1, i,j) ;
                }
            
                // normalize
                norm = sqrt(Nf_x[j][i]*Nf_x[j][i] + Nf_y[j][i]*Nf_y[j][i] ) ;
                if (norm != 0 )
                {
                     Nf_x[j][i] = Nf_x[j][i] / norm;
                     Nf_y[j][i] = Nf_y[j][i] / norm;
                 }
            } 
        } 
    } 

    // norm the vectors and save
    for (j=0; j < NBPnf ; j++)
    {
        for (i=0; i < NBPnf ; i++)
        {

           norm = sqrt(Nf_x[j][i]*Nf_x[j][i] + Nf_y[j][i]*Nf_y[j][i] ) ;
            if (norm != 0 )
            {
                Nf_x[j][i] = Nf_x[j][i] / norm ;
                Nf_y[j][i] = Nf_y[j][i] / norm ;
            }
            // save

            fprintf(foutx,"%g, %g\n",Nf_x[j][i],0.0);
            fprintf(fouty,"%g, %g\n",Nf_y[j][i],0.0);
        }
    }

    // close the files :
    fclose (foutx);
    fclose (fouty);
    return 0;
}

int disque(const char * filename, double r_a, double n_int, double k_int, 
    double n_ext, double k_ext, double Tx, double Ty)
{
    int i,j ;
    FILE *fout ;



    fout = fopen(filename,"w");
    if (fout == NULL)
    {
        printf("impossible d'ouvrir le fichier : %s\n",filename);
        exit(0);
    }
    // ecrit l'entete:
    fprintf(fout,"disque r/a = %e ; n_int=%g ; k_int=%g ; n_ext=%g ; k_ext=%g \n",r_a,n_int,k_int,n_ext,k_ext);
    fprintf(fout,"%d %d\n",NBP, NBP);
    fprintf(fout,"0 %f 0 %f\n", Tx, Ty);

    for (i = -NBP/2 ; i < NBP/2 ; i++)
    {
        for (j = -NBP/2 ; j < NBP/2 ; j++)
        {
            if (i*i+j*j < (r_a*NBP)*(r_a*NBP) )
            {
                // au coeur
                fprintf(fout,"%e, %e\n",n_int,k_int);
            }
            else
            {
                // autour
                fprintf(fout,"%e, %e\n",n_ext,k_ext);
            }
        }
    }
    fclose (fout);
    return 0;
}

int indcoeurcoquille(const char * filename, double r_a, double n_int, double k_int,
    double r1, double n1_int, double k1_int,
    double n_ext, double k_ext, double Tx, double Ty)
{
    int i,j ;
    FILE *fout ;



    fout = fopen(filename,"w");
    if (fout == NULL)
    {
        printf("impossible d'ouvrir le fichier : %s\n",filename);
        exit(0);
    }
    // ecrit l'entete:
    fprintf(fout,"coeur coquille r/a = %e ; n_int=%g ; k_int=%g ; r1/a = %e ; n1_int=%g ; k1_int=%g ; n_ext=%g ; k_ext=%g \n",r_a,n_int,k_int,r1,n1_int,k1_int,n_ext,k_ext);
    fprintf(fout,"%d %d\n",NBP, NBP);
    fprintf(fout,"0 %f 0 %f\n", Tx, Ty);

    for (i = -NBP/2 ; i < NBP/2 ; i++) {
        for (j = -NBP/2 ; j < NBP/2 ; j++) {
            if (i*i+j*j < (r_a*NBP)*(r_a*NBP)) {
                // au coeur
                fprintf(fout,"%e, %e\n",n_int,k_int);
            } else if (i*i+j*j < (r1*NBP)*(r1*NBP)) {
                // coquille
                fprintf(fout,"%e, %e\n",n1_int,k1_int);
            } else {
                // autour
                fprintf(fout,"%e, %e\n",n_ext,k_ext);
            }
        }
    }
    fclose (fout);
    return 0;
}

int rectangle(const char * filename, double r_a, double r_b, double n_int, 
    double k_int, double n_ext, double k_ext, double Tx, double Ty)
{
    int i,j ;
    FILE *fout;

    fout = fopen(filename,"w");
    if (fout == NULL)
    {
        printf("impossible d'ouvrir le fichier : %s\n",filename);
        exit(0);
    }
    // ecrit l'entete:
    fprintf(fout,"rectangle r/a = %e ; r/b = %e ; n_int=%g ; k_int=%g ; n_ext=%g ; k_ext=%g \n",r_a,r_b,n_int,k_int,n_ext,k_ext);
    fprintf(fout,"%d %d\n",NBP, NBP);
    fprintf(fout,"0 %f 0 %f\n", Tx, Ty);
    
    double r_a2=r_a/2.0*NBP;
    double r_b2=r_b/2.0*NBP;
    

    for (i = -NBP/2 ; i < NBP/2 ; ++i) {
        for (j = -NBP/2 ; j < NBP/2 ; ++j) {
            if (i<r_b2 && i>-r_b2 && j<r_a2 && j>-r_a2)
            {
                // au coeur
                fprintf(fout,"%e, %e\n",n_int,k_int);
            } else {
                // autour
                fprintf(fout,"%e, %e\n",n_ext,k_ext);
            }
        }
    }
    fclose (fout);
    return 0;
}


int main(int argc, char *argv[])
{
    double a=1.0, r_a=1.0, r_b=1.0 ;
    double r_1=1.0;
    double n_int = 3 , k_int = 0.02 ;
    double n1_int = 3 , k1_int = 0.02 ;
    double n_ext = 1 , k_ext = 0 ;
    double alpha ;
    FILE *fout , *fout2,  *fout3 ;
    int i, j ;
    char fichier_a_ecrire[200];
    char fichier_a_ecrire2[200];
    char fichier_a_ecrire3[200];
    char disquerect[200];
    double Tx, Ty;
    
    tc typecalc;
    
    // This should never happen
    if(argc<1)
        return 1;
    
    printf("------------------------------------------------------------------------------\n");
    printf("   Normal field calculation (disc or rectangle)\n");
    printf("   Jérôme Michallon, Davide Bucci, october 2013\n");
    printf("       type %s -h to get an help\n", argv[0]);
    printf("------------------------------------------------------------------------------\n\n");
    
    if(argc==2 && strcmp(argv[1],"-h")==0) {
        printf("   This is a program to calculate some refractive index distibutions,\n"
               "   as well as the normal field vectors, by using the strategy published\n"
               "   in Götz et al., Optics express vol. 16, no. 22, pp. 17295-17301, 2008\n\n"
               "   Several geometrical shapes are available: disc, rectangle, core/shell.\n"
               "   You might use the software in the interactive mode, not specifying any\n"
               "   command line parameter. You can also specify everything in the command\n"
               "   line. This might be useful if you need it in a script.\n"
               "   Three files are genrated: one with the refractive index structure.\n"
               "   Two with the x and y component of the normal vector field (their name\n"
               "   ends with '_nvf_x' and '_nvf_y' respectively).\n"
               "   All files are in the OptiWave RID format.\n\n"
               "   The command line options begins as follows:\n"
               "   %s Tx Ty type filename ...\n"
               "   Where options following the filename depend on the geometry specified. \n\n"
               "   (d) type=d (disk), the command line is as follows:\n"
               "      %s Tx Ty d filename r_x n_int k_int n_ext k_ext\n"
               "      where:  Tx is the horizontal size of the calculation window\n"
               "              Ty is the vertical size of the calculation window\n"
               "              d  specifies a disk structure.\n"
               "              filename is the name of the files (do not use an extension).\n"
               "              r_x is the radius over Tx ratio.\n"
               "              n_int and k_int are the real and imaginary part of refractive\n"
               "                  index to be used INSIDE the disk\n"
               "              n_ext and k_ext are the real and imaginary part of refractive\n"
               "                  index to be used OUTSIDE the disk\n\n"
               "   (r) type=r (rectangle), the command line is as follows:\n"
               "      %s Tx Ty r filename r_x r_y n_int k_int n_ext k_ext\n"
               "      where:  Tx is the horizontal size of the calculation window\n"
               "              Ty is the vertical size of the calculation window\n"
               "              r  specifies a rectangle structure.\n"
               "              filename is the name of the files (do not use an extension).\n"
               "              r_x is the width over Tx ratio.\n"
               "              r_y is the height over Ty ratio.\n"
               "              n_int and k_int are the real and imaginary part of refractive\n"
               "                  index to be used INSIDE the rectangle\n"
               "              n_ext and k_ext are the real and imaginary part of refractive\n"
               "                  index to be used OUTSIDE the rectangle\n\n"
               "   (r) type=c (core/shell), the command line is as follows:\n"
               "      %s Tx Ty r filename r_int n_int k_int r1 n1_int k1_int n_ext k_ext\n"
               "      where:  Tx is the horizontal size of the calculation window\n"
               "              Ty is the vertical size of the calculation window\n"
               "              c  specifies a core/shell structure.\n"
               "              filename is the name of the files (do not use an extension).\n"
               "              r_int is the core radius over Tx ratio.\n"
               "              n_int and k_int are the real and imaginary part of refractive\n"
               "                  index to be used INSIDE the core\n"
               "              r1 is the shell over Tx ratio.\n"
               "              n1_int and k1_int are the real and imaginary part of refractive\n"
               "                  index to be used INSIDE the shell\n"
               "              n_ext and k_ext are the real and imaginary part of refractive\n"
               "                  index to be used OUTSIDE the shell\n\n"
               "   NOTE: for a lossy structure, all the imaginary parts of refractive indices\n"
               "   should be negative or non zero.\n\n"
               "   OUTPUT: three output files will be genreated. The first one contains the\n"
               "   refractive index distribution and it is called with the same name specified\n"
               "   as an input parameter. The second one is the x component of the normal \n"
               "   vector field and has '_nvf_x' added at the end of the filename. The third \n"
               "   is the y component of the vector field and ends with '_nvf_y'\n"
               , argv[0],argv[0],argv[0],argv[0]);
               return 0;
    }
               
    if (argc<10) {
        printf("Horizontal size of the calculation window? ");
        scanf("%lf", &Tx);
        printf("Vertical size of the calculation window? ");
        scanf("%lf", &Ty);
        
        /*tout est tap    la main*/
        printf("[d]isk, [r]ectangle, [c]core/shell structure? ");
        scanf("%s", disquerect);
        if(tolower(disquerect[0])=='d') {
            typecalc=disc;
        } else if(tolower(disquerect[0])=='r') {
            typecalc=rect;
        } else if(tolower(disquerect[0])=='c') {
            typecalc=coeurcoquille;
        } else {
            fprintf(stderr, "Uh? I do not know how to deal with that...");
            return 1;
        }

        printf("Name of the output file? ");
        scanf("%s", fichier_a_ecrire);

        if (typecalc==disc) {
            printf("Ratio (r/a = radius/calculation window width)? ");
            scanf("%lf", &r_a); 
        } else if (typecalc==coeurcoquille) {
            printf("Ratio (r/a = core radius/calculation window width)? ");
            scanf("%lf", &r_a); 
            printf("Ratio (r1/a = shell radius/calculation window width)?");
            scanf("%lf", &r_1);
            printf("Real part of the shell index? ");
            scanf("%lf", &n1_int);  
            printf("Imaginary part of the shell index? ");
            scanf("%lf", &k1_int);      
        } else {
            printf("Ratio (r/a = width/calculation window width)? ");
            scanf("%lf", &r_a);     
            printf("Ratio (r/b = height/calculation window height)? ");
            scanf("%lf", &r_b);     
        }

        printf("Real part of the internal refractive index? ");
        scanf("%lf", &n_int);   
        printf("Imaginary part of the internal refractive index? ");
        scanf("%lf", &k_int);       


        printf("Real part of the external refractive index? ");
        scanf("%lf", &n_ext);   
        printf("Imaginary part of the external refractive index? ");
        scanf("%lf", &k_ext);
    } else {
        /*tout est sp cifi  au programme*/

        int count = 1;
        
        Tx=atof(argv[count++]);
        Ty=atof(argv[count++]);
        
        if(tolower(argv[count][0])=='d') {
            typecalc=disc;
            if (argc!=10) {
                fprintf(stderr, "Invalid number of parameters\n");
                return 1;
            }
        } else if(tolower(argv[count][0])=='r') {
            typecalc=rect;
            if (argc!=11) {
                fprintf(stderr, "Invalid number of parameters\n");
                return 1;
            }
        } else if(tolower(argv[count][0])=='c') {
            typecalc=coeurcoquille;
            if (argc!=13) {
                fprintf(stderr, "Invalid number of parameters\n");
                return 1;
            }
        }   else {
            fprintf(stderr, "Unrecognized type %s\n", argv[count]);
            return 1;
        }
        ++count;
        sscanf(argv[count++],"%s", fichier_a_ecrire);
        r_a=atof(argv[count++]);
        if (typecalc==rect) {
            printf("Rectangle: ");
            r_b=atof(argv[count++]);
        } else {
            printf("Disc: ");
        }
        
        n_int=atof(argv[count++]);
        k_int=atof(argv[count++]);
        
        if (typecalc==coeurcoquille) {
            r_1=atof(argv[count++]);
            n1_int=atof(argv[count++]);
            k1_int=atof(argv[count++]);
            printf("r_1=%g, n1_int=%g, k1_int=%g\n", r_1, n1_int, k1_int);
        }
        n_ext=atof(argv[count++]);
        k_ext=atof(argv[count++]);

        printf("Tx=%g, Ty=%g, Type: %s\n", Tx, Ty, 
            (typecalc==rect)?"rectangle":
            (typecalc==disc)?"disc":"core/shell");
        printf("Fichier %s, r/a =%g, n=%g, k=%g, n=%g, k=%g\n",
            fichier_a_ecrire, r_a, n_int, k_int, n_ext, k_ext);
    
    }


    // calcul le champ d'indice :
    printf("Creation of refractive index file\n");
    
    if (typecalc==disc) {
        disque(fichier_a_ecrire, r_a,  n_int,  k_int,  n_ext,  k_ext, Tx, Ty);
    } else if(typecalc==coeurcoquille) {
        indcoeurcoquille(fichier_a_ecrire, r_a, n_int, k_int, r_1, n1_int, k1_int,
            n_ext, k_ext, Tx, Ty);
    } else {
        rectangle(fichier_a_ecrire, r_a, r_b,  n_int,  k_int,  n_ext,  
            k_ext, Tx, Ty);
    }

    // champ normal :
    printf("Creation of normal field fiels\n");
    snprintf(fichier_a_ecrire2,199,"%s_nvf_x", fichier_a_ecrire);
    snprintf(fichier_a_ecrire3,199,"%s_nvf_y", fichier_a_ecrire);
     
    normal_field_generation(fichier_a_ecrire2,fichier_a_ecrire3, 
        typecalc, r_a, r_b, Tx, Ty);


    return 0;
}

