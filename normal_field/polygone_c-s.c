#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#define NBP 1024    // multiple de 2
#define NBPnf 257   // puissance de 2 +1
#define PI 3.141592654
#define NB_COUCHES_MAX 30
#define NB_COTES_MAX 100
#define precision 1e-5

////////////DOIT ETRE COMPILE AVEC -O3 (sous redhat)/////////////////////////

// dessine la structure voulue.
// Permet de dessiner toutes les formes géometriques classiques en calculant
// un polygone inscrit dans une cercle.


/* Normal field calculation based on the work of Götz et al.
    BibTeX reference:


 @article{gotz2008normal,
  title={Normal vector method for the RCWA with automated vector field
  generation},
  author={G{\"o}tz, P. and Schuster, T. and Frenner, K. and Rafler, S. and
  Osten, W.},
  journal={Optics express},
  volume={16},
  number={22},
  pages={17295--17301},
  year={2008},
  publisher={Optical Society of America}
 }

*/


double abs2(double a)
{
    if (a<0)
    {
        return (-a);
    }
        return (a);
}

int appartient_triangle(double triangle_x[3],double triangle_y[3],double x,
    double y)
{
    double aire1, aire2;
    int i;

    aire1 = (triangle_x[1]-triangle_x[0])*(triangle_y[2]-triangle_y[0]) -
        (triangle_y[1]-triangle_y[0])*(triangle_x[2]-triangle_x[0]);

    aire2 = 0;
    for (i=0; i<3; i++)
        aire2 += abs2((triangle_x[i]-x)*(triangle_y[(i+1)%3]-y) -
        (triangle_y[i]-y)*(triangle_x[(i+1)%3]-x));
            // somme les aires des triangles


    if (aire1+ precision > aire2 && aire1 - precision < aire2)
        return 1;
    return 0;
}

int interpolate (double Nx[NBPnf][NBPnf],double Ny[NBPnf][NBPnf],
    int centerx[NB_COUCHES_MAX], int centery[NB_COUCHES_MAX],
    int Polygonex[NB_COTES_MAX][NB_COUCHES_MAX] ,
    int Polygoney[NB_COTES_MAX][NB_COUCHES_MAX], int nb_cote[NB_COUCHES_MAX] ,
    int nb_couches, int i1,int j1,int i2,int j2,int i3,int j3,int i4,int j4,
    int m,int n)
{

    int nbp, k,l,o,  i,j;
    double alpha;
    Nx[n][m] = 0;
    Ny[n][m] = 0;
    double norm;
    double ax,bx,cx;
    double ay,by,cy;
    // weight the point from the radius:
    // for each shape
    for (l = 0; l < nb_couches; l++)
    {
        // for each side of the shape
        for (o = 0; o < nb_cote[l]; o++)
        {

            // vector a along the side
            ax = Polygonex[(o+1)%nb_cote[l]][l] -Polygonex[o][l];
            ay = Polygoney[(o+1)%nb_cote[l]][l] -Polygoney[o][l];

            // vector b perpendicular pointing to outside of the
            // calculation window
            if (ay !=0)
            {

                // a vector perpendicular to the side
                bx = 1;
                by = -ax/ay*bx;

                // normalisation
                bx = bx / (sqrt (bx*bx + by*by));
                by = by / (sqrt (bx*bx + by*by));
            }
            else
            {
                bx = 0;
                by = 1;
            }

            // a vector pointing form the origin to the firt point of the shape
            cx = Polygonex[o][l] + centerx[l] - NBPnf/2;
            cy = Polygoney[o][l] + centery[l] - NBPnf/2;


            // inverse the vector to make it pointing to the outside of the
            // calculation window
            if (cx*bx + by*cy < 0)
            {
                bx = -bx;
                by = -by;
            }

            // number of points that will be calculated along the side
            nbp = abs2(ax);
            if (abs2(ay) > nbp)
                nbp = abs2(ay);

            // for each point along the side:
            for (k=0; k < nbp; k++)
            {
                i = centerx[l] + Polygonex[o][l] + (int)((double)k/nbp*ax);
                j = centery[l] + Polygoney[o][l] + (int)((double)k/nbp*ay);

                if (i == m && j == n)
                {
                     Nx[n][m] = bx;
                     Ny[n][m] = by;
                     return 1;
                }
                else
                {
                    Nx[n][m] = Nx[n][m] + bx/ abs(double(pow(i-m,2) +
                        pow(j-n,2)));
                    Ny[n][m] = Ny[n][m] + by/ abs(double(pow(i-m,2) +
                        pow(j-n,2)));
                }
            }

        }

    }
    // weigt the 4 points :
    Nx[n][m] = Nx[n][m] + Nx[j1][i1]/ abs(double(pow(i1-m,2) + pow(j1-n,2)));
    Ny[n][m] = Ny[n][m] + Ny[j1][i1]/ abs(double(pow(i1-m,2) + pow(j1-n,2)));

    Nx[n][m] = Nx[n][m] + Nx[j2][i2]/ abs(double(pow(i2-m,2) + pow(j2-n,2)));
    Ny[n][m] = Ny[n][m] + Ny[j2][i2]/ abs(double(pow(i2-m,2) + pow(j2-n,2)));

    Nx[n][m] = Nx[n][m] + Nx[j3][i3]/ abs(double(pow(i3-m,2) + pow(j3-n,2)));
    Ny[n][m] = Ny[n][m] + Ny[j3][i3]/ abs(double(pow(i3-m,2) + pow(j3-n,2)));

    Nx[n][m] = Nx[n][m] + Nx[j4][i4]/ abs(double(pow(i4-m,2) + pow(j4-n,2)));
    Ny[n][m] = Ny[n][m] + Ny[j4][i4]/ abs(double(pow(i4-m,2) + pow(j4-n,2)));


    // normalize the vectors :
    norm = sqrt(pow(Nx[n][m],2) + pow(Ny[n][m],2));
    if (norm !=0)
    {
        Nx[n][m]= Nx[n][m] /norm;
        Ny[n][m] = Ny[n][m]/norm;
    }
    return 1;
}
int normal_field_generation(const char * filename_x,const char * filename_y, double x_poly[NB_COTES_MAX][NB_COUCHES_MAX] , double y_poly[NB_COTES_MAX][NB_COUCHES_MAX], int nb_cote[NB_COUCHES_MAX], double origineX[NB_COUCHES_MAX], double origineY[NB_COUCHES_MAX], int nb_couches, double tailleY)
{
    int dx,dy;
    int R[NB_COUCHES_MAX],Cx[NB_COUCHES_MAX],Cy[NB_COUCHES_MAX];
    int Polygonex[NB_COTES_MAX][NB_COUCHES_MAX],
        Polygoney[NB_COTES_MAX][NB_COUCHES_MAX];
    int i,j,k,l,m,n;
    double Nf_x[NBPnf][NBPnf] ,Nf_y[NBPnf][NBPnf];
    double alpha, nbp;
    double norm;
    FILE *foutx,*fouty;

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
    fprintf(foutx,"normal field x ");
    for (j=0; j < nb_couches; j++)
        fprintf(foutx,"{ %d-cote; x=%g; y=%g } ",nb_cote[j] ,
            origineX[j],origineY[j]);
    fprintf(foutx,"\n");

    fprintf(foutx,"%d %d\n",NBPnf, NBPnf);
    fprintf(foutx,"0 1 0 1\n");

    fprintf(fouty,"normal field y ");
    for (j=0; j < nb_couches; j++)
        fprintf(fouty,"{ %d-cote; x=%g; y=%g } ",
            nb_cote[j],origineX[j],origineY[j]);
    fprintf(fouty,"\n");

    fprintf(fouty,"%d %d\n",NBPnf, NBPnf);
    fprintf(fouty,"0 1 0 %g\n", tailleY);

    // initialize the array :
    for (i=1; i<NBPnf-1;i++)
    {
        for (j=1; j<NBPnf-1;j++)
        {
            Nf_x[j][i]= 0.0;
            Nf_y[j][i]= 0.0;
        }
    }
    //create the field at the boundaries of the windows
    for (k=1; k<NBPnf;k++)
    {
        Nf_x[0][k]= 0.0;
        Nf_y[0][k]= -1.0;

        Nf_x[NBPnf-1][k]= 0.0;
        Nf_y[NBPnf-1][k]= 1.0;
    }
    for (k=1; k<NBPnf;k++)
    {
        Nf_x[k][0]= -1.0;
        Nf_y[k][0]= 0;

        Nf_x[k][NBPnf-1]= 1.0;
        Nf_y[k][NBPnf-1]= 0;
    }
    Nf_x[0][0]= -1.0;
    Nf_y[0][0]= -1.0;

    Nf_x[NBPnf-1][NBPnf-1]= 1.0;
    Nf_y[NBPnf-1][NBPnf-1]= 1.0;

    Nf_x[0][NBPnf-1]= 1.0;
    Nf_y[0][NBPnf-1]= -1.0;

    Nf_x[NBPnf-1][0]= -1.0;
    Nf_y[NBPnf-1][0]= 1.0;
    // set the value at the middle of the windows
    Nf_x[NBPnf/2][NBPnf/2]= 0.0;
    Nf_y[NBPnf/2][NBPnf/2]= 0.0;

    // set the value of the radius and circle center in term of pixel:
    for (i = 0; i < nb_couches; i++)
    {
        //R[i] = (int)(r_a[i]* (double)(NBPnf));
        Cx[i] = origineX[i] * NBPnf + NBPnf/2;
        Cy[i] = origineY[i] * NBPnf + NBPnf/2;
        for (j=0; j < nb_cote[i]; j++)
        {
            Polygonex[j][i] = (int)(x_poly[j][i]* (double)(NBPnf));
            Polygoney[j][i] = (int)(y_poly[j][i]* (double)(NBPnf));

        }
    }
    //interpolate the field
    dx = NBPnf/2;
    dy = NBPnf/2;

    while (dx > 2 && dy > 2)
    {
        // first rectangle
        //initialize indexes
        i=0;
        k = dx;
        m = dx/2;
        while (k < NBPnf)
        {
            j = 0;
            l = dy;
            n = dy/2;
            while (l < NBPnf)
            {
                interpolate(Nf_x,Nf_y, Cx,Cy, Polygonex,Polygoney,nb_cote,
                    nb_couches, i,j, i,l, k,j, k,l, m,n);
                //shift
                j = l;
                l = l + dy;
                n = n + dy;
            }
            // shift:
            i = k;
            k = k + dx;
            m = m + dx;
        }

        //tilted rectangle

        //initialize indexes
        i=0;
        k = 0;
        while (i < NBPnf-2*dx)
        {
            j = 0;
            l = 0;
            while (j < NBPnf-2*dy)
            {

                interpolate(Nf_x,Nf_y, Cx,Cy, Polygonex,Polygoney,nb_cote,
                    nb_couches, k+dx/2,l+dy/2, i,l+dy, k+dx,l+dy,
                    k+dx/2,l+dy+dy/2, k+dx/2,l+dy);
                interpolate(Nf_x,Nf_y, Cx,Cy, Polygonex,Polygoney,nb_cote,
                    nb_couches, k+dx/2,l+dy/2, k+dx,j, k+dx,l+dy,
                    k+dx/2+dx,l+dy/2, k+dx,l+dy/2);
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
    for (i=1; i < NBPnf-1; i++)
    {
        for (j=1; j < NBPnf-1; j++)
        {

            if (Nf_x[j][i] == 0 && Nf_y[j][i] == 0)
            {

                interpolate(Nf_x,Nf_y, Cx,Cy ,Polygonex,Polygoney,nb_cote,
                    nb_couches, i-1,j, i+1,j, i,j+1, i,j-1, i,j);

                // normalize
                norm = sqrt(Nf_x[j][i]*Nf_x[j][i] + Nf_y[j][i]*Nf_y[j][i]);
                if (norm != 0)
                {
                     Nf_x[j][i] = Nf_x[j][i] / norm;
                     Nf_y[j][i] = Nf_y[j][i] / norm;
                 }
            }
        }
    }
    // create onto the controur the value of the normal field.
    /*
    nbp = 7*R; // 2*pi radius
    for (k=0; k < nbp; k++)
    {
        alpha = -PI + (double)k/nbp*2*PI;

        i =Cx + int(double(R)* cos(alpha));
        j =Cy + int(double(R)* sin(alpha));

        Nf_x[j][i] = cos(alpha);
        Nf_y[j][i] = sin(alpha);
    }
   //create the field at the boundaries of the windows
    for (k=1; k<NBPnf;k++)
    {
        Nf_x[0][k]= 0.0;
        Nf_y[0][k]= -1.0;

        Nf_x[NBPnf-1][k]= 0.0;
        Nf_y[NBPnf-1][k]= 1.0;
    }
    for (k=1; k<NBPnf;k++)
    {
        Nf_x[k][0]= -1.0;
        Nf_y[k][0]= 0;

        Nf_x[k][NBPnf-1]= 1.0;
        Nf_y[k][NBPnf-1]= 0;
    }
    Nf_x[0][0]= -1.0;
    Nf_y[0][0]= -1.0;

    Nf_x[NBPnf-1][NBPnf-1]= 1.0;
    Nf_y[NBPnf-1][NBPnf-1]= 1.0;

    Nf_x[0][NBPnf-1]= 1.0;
    Nf_y[0][NBPnf-1]= -1.0;

    Nf_x[NBPnf-1][0]= -1.0;
    Nf_y[NBPnf-1][0]= 1.0;  */

    // norm the vectors and save
    for (j=0; j < NBPnf; j++)
    {
        for (i=0; i < NBPnf; i++)
        {

           norm = sqrt(Nf_x[j][i]*Nf_x[j][i] + Nf_y[j][i]*Nf_y[j][i]);
            if (norm != 0)
            {
                Nf_x[j][i] = Nf_x[j][i] / norm;
                Nf_y[j][i] = Nf_y[j][i] / norm;
            }
            // save

            fprintf(foutx,"%g, %g\n",Nf_x[j][i],0.0);
            fprintf(fouty,"%g, %g\n",Nf_y[j][i],0.0);
        }
    }

    // close the files :
    fclose (foutx);
    fclose (fouty);
}

int polygone(const char * filename, int nb_couches, double origineX[NB_COUCHES_MAX], double origineY[NB_COUCHES_MAX],double r_a[NB_COUCHES_MAX],double r_b[NB_COUCHES_MAX],double n_int[NB_COUCHES_MAX] , double k_int[NB_COUCHES_MAX], double n_ext, double k_ext, double x_poly[NB_COTES_MAX][NB_COUCHES_MAX] , double y_poly[NB_COTES_MAX][NB_COUCHES_MAX], int nb_cote[NB_COUCHES_MAX], double tailleY)
{
    FILE *fout;
    int i, j, k,l;
    int flag;
    double trianglex[3], triangley[3];

    fout = fopen(filename,"w");
    if (fout == NULL)
    {
        printf("impossible d'ouvrir le fichier : %s\n",filename);
        exit(0);
    }


     // genere les polygone
    fprintf(fout,"polygone_c/s ");
    for (j=0; j < nb_couches; j++)
        fprintf(fout,"{ %d-cote r/a = %e; x=%g y=%g; n_int=%g; k_int=%g } ",nb_cote[j],r_a[j],origineX[j],origineY[j],n_int[j],k_int[j]);
    fprintf(fout,"; n_ext=%g; k_ext=%g \n",n_ext,k_ext);

    fprintf(fout,"%d %d\n",NBP, NBP);
    fprintf(fout,"0 1 0 %g\n", tailleY);

    for (j = -NBP/2; j < NBP/2; j++)
    {
        for (i = -NBP/2; i < NBP/2; i++)
        {
            flag = nb_couches;
            for (l=nb_couches-1; l >=0; l--)      // parcours les couches de la couche externe vers la couche interne
            {
//              if ((i-origineX[l]*NBP)*(i-origineX[l]*NBP)+(j-origineY[l]*NBP)*(j-origineY[l]*NBP) < (r_a[l]*NBP)*(r_a[l]*NBP))
                if ((i-origineX[l]*NBP)/(r_a[l]*NBP)*(i-origineX[l]*NBP)/(r_a[l]*NBP)+(j-origineY[l]*NBP)/(r_b[l]*NBP)*(j-origineY[l]*NBP)/(r_b[l]*NBP) < 1)
                {
                    // au coeur
                    // verifie si on est dans le polygone ou non:
                    for (k=0; k<nb_cote[l]; k++)
                    {
                        trianglex[0] = origineX[l]*NBP + 0;
                        triangley[0] = origineY[l]*NBP + 0;

                        trianglex[1] = origineX[l]*NBP + x_poly[k][l] *NBP;
                        triangley[1] = origineY[l]*NBP + y_poly[k][l] *NBP;

                        trianglex[2] = origineX[l]*NBP + x_poly[(k+1)%nb_cote[l]][l] *NBP;
                        triangley[2] = origineY[l]*NBP + y_poly[(k+1)%nb_cote[l]][l] *NBP;

                        if (appartient_triangle(trianglex,triangley,i,j))
                        {
                            k = nb_cote[l]+2; // fin
                            flag = l;
                        }

                    }

                }
            }
            if (flag ==nb_couches)      // on n'est dans aucune des couches:
                fprintf(fout,"%e, %e\n",n_ext,k_ext);
            else            // apartient  au polygone
                fprintf(fout,"%e, %e\n",n_int[flag],k_int[flag]);
        }
    }
    fclose (fout);
    return 0;
}

int main(int argc, char *argv[])
{
    double a=1;
    double n_ext = 1 , k_ext = 0;
    int i, j;


    char fichier_a_ecrire[200];
    char fichier_a_ecrire2[200];
    char fichier_a_ecrire3[200];
    double alpha;
    int nb_cote_max;
    int nb_couches= 1;
    double tailleY=1.0;
    double tailleX=1.0;

    double r_a[NB_COUCHES_MAX];
    double r_b[NB_COUCHES_MAX];
    double n_int[NB_COUCHES_MAX] , k_int[NB_COUCHES_MAX], alpha0[NB_COUCHES_MAX];
    double alphar[NB_COUCHES_MAX], origineX[NB_COUCHES_MAX],origineY[NB_COUCHES_MAX];
    int nb_cote[NB_COUCHES_MAX];


    printf("------------------------------------------------------------------\n");
    printf("   Normal field calculation (regular shapes)\n");
    printf("   Jérôme Michallon, Davide Bucci, April 2013\n");
    printf("------------------------------------------------------------------\n\n");


    if (argc==1)
    {
        /*tout est tapé à la main*/

        printf("Nom du fichier à écrire\n");
        scanf("%s", fichier_a_ecrire);
        //printf("%s", fichier_a_ecrire);

        /*printf("taille de la cellule (a)\n");
        scanf("%e", a);*/
        printf("taille de la fenetre en y (dimension réduite tailleY/tailleX) ?\n");
        scanf("%lf", &tailleY);
        printf("nombre de couches ?\n");
        scanf("%d", &nb_couches);

        /*double r_a[nb_couches];
        double n_int[nb_couches] , k_int[nb_couches], alpha0[nb_couches];
        int nb_cote[nb_couches];*/
        printf("\ndéfinition des couches, la couche 1 étant la couche centrale !\n");

        for (i = 0; i < nb_couches; i++)
        {
            printf("\nCOUCHE%d\n",i+1);
            printf("rapport (r/a = rayon / taille de la maille x)\n");
            scanf("%lf", &r_a[i]);
            r_b[i] = r_a[i];
            printf("nombre de cotes \n");
            scanf("%d", &nb_cote[i]);
            printf("origine x (entre 0 et 1) \n");
            scanf("%lf", &origineX[i]);
            printf("origine y (entre 0 et %g) \n", tailleY);
            scanf("%lf", &origineY[i]);
            printf("angle à l'origine (deg)\n");
            scanf("%lf", &alpha0[i]);
            alpha0[i] = alpha0[i]/180*PI;  //rad

            printf("indice réel interne\n");
            scanf("%lf", &n_int[i]);
            printf("indice imaginaire interne\n");
            scanf("%lf", &k_int[i]);
        }

        printf("indice réel externe\n");
        scanf("%lf", &n_ext);
        printf("indice imaginaire externe\n");
        scanf("%lf", &k_ext);
    }
    else
    {
        j = 1;
        sscanf(argv[j++],"%s", fichier_a_ecrire);

        printf("\n fichier %s  ",fichier_a_ecrire);

        tailleY = atof(argv[j++]);
        nb_couches=atoi(argv[j++]);
        for (i = 0; i < nb_couches; i++)
        {
            r_a[i]=atof(argv[j++]);
            r_b[i] = r_a[i];
            nb_cote[i]=atoi(argv[j++]);
            origineX[i]=atof(argv[j++]);
            origineY[i]=atof(argv[j++]);
            alpha0[i] = atof(argv[j++]);
            alpha0[i] = alpha0[i]/180*PI;  //rad
            n_int[i]=atof(argv[j++]);
            k_int[i]=atof(argv[j++]);
            printf("{ %d-cote r/a = %e;x=%g y=%g; n_int=%g; k_int=%g } ",nb_cote[i],r_a[i],origineX[i],origineY[i],n_int[i],k_int[i]);
        }
        n_ext=atof(argv[j++]);
        k_ext=atof(argv[j++]);

        printf("; n_ext=%g; k_ext=%g \n",n_ext,k_ext);

    }

        // convertie origine de -1 à 1 :
    for (j=0; j < nb_couches; j++)
    {
        origineX[j] = (origineX[j]-0.5);
        origineY[j] = (origineY[j]-0.5);
        printf("origine = %g %g \n",origineX[j],origineY[j]);
    }

    // coordonnées des points du polygone  entre -0.5 et 0.5:

    double x_poly[NB_COTES_MAX][NB_COUCHES_MAX] ,y_poly[NB_COTES_MAX][NB_COUCHES_MAX];

    // genere les points du polygone :
    for (j=0; j < nb_couches; j++)
    {


        alpha = 2*PI/nb_cote[j];
        for (i = 0; i < nb_cote[j]; i++)
        {


            // dessine le polygone centré en zéro
            x_poly[i][j] = r_a[j]* cos (alpha0[j]+alpha*i);
            y_poly[i][j] = r_b[j]/tailleY * sin (alpha0[j]+alpha*i);


        }


    }
    // calcul le champ d'indice :
    printf("creation du fichier d'indice\n");
    polygone(fichier_a_ecrire, nb_couches, origineX, origineY,r_a,r_b, n_int ,k_int, n_ext, k_ext, x_poly, y_poly, nb_cote, tailleY);

    // champ normal :
    printf("creation du fichier de champ normal\n");
    sprintf(fichier_a_ecrire2,"%s_nvf_x", fichier_a_ecrire);
    sprintf(fichier_a_ecrire3,"%s_nvf_y", fichier_a_ecrire);

    normal_field_generation(fichier_a_ecrire2,fichier_a_ecrire3, x_poly, y_poly, nb_cote, origineX, origineY, nb_couches, tailleY);

    return 0;
}

