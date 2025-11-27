function [Nx Ny] = interpolate (Nx,Ny,centerx, centery,radiusx, radiusy, i1,j1,i2,j2,i3,j3,i4,j4,m,n)

Nx(n,m) = 0;
Ny(n,m) = 0;

%tolerance = 0.05 ;
% weight the circle and the four point described by i and j to get the
% normal field

nbp = 7*radiusx ;
for k=1:nbp
    alpha = -pi + k/nbp*2*pi ;

    i =centerx + floor(radiusx* cos(alpha)) ;
    j =centery + floor(radiusy* sin(alpha)) ;
    if (i == 0)
        i = 1 ;
    end
    if (j == 0)
        j = 1 ;
    end

    if (i == m && j == n)
         Nx(n,m) = cos(alpha) ;
         Ny(n,m) = sin(alpha) ;
         break ;
    else
        Nx(n,m) = Nx(n,m) + cos(alpha)/ abs( (i-m)^2 + (j-n)^2 );                       
        Ny(n,m) = Ny(n,m) + sin(alpha)/ abs( (i-m)^2 + (j-n)^2 );     
    end 

end

% weigt the 4 points :
Nx(n,m) = Nx(n,m) + Nx(j1,i1)/ abs( (i1-m)^2 + (j1-n)^2 ); 
Ny(n,m) = Ny(n,m) + Ny(j1,i1)/ abs( (i1-m)^2 + (j1-n)^2 ); 

Nx(n,m) = Nx(n,m) + Nx(j2,i2)/ abs( (i2-m)^2 + (j2-n)^2 );                       
Ny(n,m) = Ny(n,m) + Ny(j2,i2)/ abs( (i2-m)^2 + (j2-n)^2 );  

Nx(n,m) = Nx(n,m) + Nx(j3,i3)/ abs( (i3-m)^2 + (j3-n)^2 );                       
Ny(n,m) = Ny(n,m) + Ny(j3,i3)/ abs( (i3-m)^2 + (j3-n)^2 );  

Nx(n,m) = Nx(n,m) + Nx(j4,i4)/ abs( (i4-m)^2 + (j4-n)^2 );                       
Ny(n,m) = Ny(n,m) + Ny(j4,i4)/ abs( (i4-m)^2 + (j4-n)^2 ); 


% normalize the vectors :
norm = sqrt(Nx(n,m)^2 + Ny(n,m)^2) ;
if (norm ~=0)
    Nx(n,m) = Nx(n,m) /norm ;                       
    Ny(n,m) = Ny(n,m)/norm ;
end

