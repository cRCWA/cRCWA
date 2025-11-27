function create_nf_halfsphere_jm_interpolated

    clc ;
	% all sizes are given in micrometers

    % number of point must be a power of 2
	totx=32
	toty=32
 
  
    n_sub = 1 
	n_ring = 3.9058
    k_ring = -0.022
    n_air = 1 
	
   
 	tot_x_size = 350
 	tot_y_size = 350
    
    R_x = 315/2
    R_y = 315/2

    % decalage=R_x-tot_x_size/2;
    decalage=0;
    
    % create the cylinder
	dx=tot_x_size/totx;
 	dy=tot_y_size/toty;
	A=n_sub*ones(totx, toty);       % real part of the refractive index
    B =zeros(totx, toty) ;          % imaginary part of the refractive index
	for i=1:totx
		for j=1:toty
			x1=(i-totx/2)*dx;
			y=(j-toty/2)*dy;
			
            if ((y/R_y)^2+(x1/R_x)^2 <=1)  %R�alisation d'une ellipse
				A(j,i)=n_ring;
                B(j,i)= k_ring;
            end;
        end
    end

    % create the field at the boundaries of the windows
    for k=2:totx-1
	    Nf_x(1,k)= 0.0;                       
        Nf_y(1,k)= -1.0 ; %-1
		
        Nf_x(toty,k)= 0.0 ;                       
        Nf_y(toty,k)= 1.0; %1
    end
    for k=2:toty-1
	    Nf_x(k,1)= -1.0 ;    %-1                   
        Nf_y(k,1)= 0 ;
		
        Nf_x(k,totx)= 1.0;                       
        Nf_y(k,totx)= 0 ;
    end
    Nf_x(1,1)= -1.0;                       
    Nf_y(1,1)= -1.0 ; 
		
    Nf_x(toty,totx)= 1.0 ;                       
    Nf_y(toty,totx)= 1.0; 
    
    Nf_x(1,totx)= 1.0 ;                       
    Nf_y(1,totx)= -1.0;
    
    Nf_x(toty,1)= -1.0 ;                       
    Nf_y(toty,1)= 1.0;  
    % set the value at the middle of the windows
    Nf_x(toty/2,totx/2)= 0.0;                       
    Nf_y(toty/2,totx/2)= 0.0;
    
    
    % set the value of the radius and circle center in term of pixel:
    Rx = R_x / tot_x_size * totx 
    Ry = R_y / tot_y_size * totx 
    Cx = totx/2 ;
    Cy = toty/2 ;
    
    %interpolate the field
    dx = totx/2 ;
    dy = toty/2 ;
    
    while (dx > 2 && dy > 2)
        % first rectangle
        %initialize indexes
        i=1;
        k = dx ;
        m = dx/2;
        while ( k <= totx)
            j = 1;
            l = dy ;
            n = dy/2 ;
            while (l <= toty)
                [Nf_x Nf_y] = interpolate(Nf_x,Nf_y, Cx,Cy, Rx,Ry, i,j, i,l, k,j, k,l, m,n) ;
                %shift 
                j = l ;
                l = l + dy ;
                n = n + dy ;
            end 
            % shift:
            i = k ;
            k = k + dx ;
            m = m + dx ;
        end

        %tilted rectangle
        if (dx >=4 || dy >= 4)
            %initialize indexes
            i=0;
            while ( i < totx-2*dx)
                j = 0;
                while (j < toty-2*dy)
                    if (i == 0)
                        k = 1;
                    else
                        k = i ;
                    end
                    if (j == 0)
                        l = 1;
                    else
                        l = j ;
                    end
                    [Nf_x Nf_y] = interpolate(Nf_x,Nf_y, Cx,Cy ,Rx,Ry, i+dx/2,j+dy/2, k,j+dy, i+dx,j+dy, i+dx/2,j+dy+dy/2, i+dx/2,j+dy) ;
                    [Nf_x Nf_y] = interpolate(Nf_x,Nf_y, Cx,Cy, Rx,Ry, i+dx/2,j+dy/2, i+dx,l, i+dx,j+dy, i+dx/2+dx,j+dy/2, i+dx,j+dy/2) ;
                    %shift 
                    j = j + dy ;
                end 
                % shift:
                i = i + dx ;
            end
        end
        
        % refine the interpolation:
        dx = dx/2 ;
        dy = dy/2 ;
         
    end

    % for all remaining zero vectors comming from discretization,
    % interpolation is performed
    for i=2:totx
		for j=2:toty
            if (Nf_x(j,i) == 0 && Nf_y(j,i) == 0 )
                Nf_x(j,i)= (Nf_x(j,i-1) +Nf_x(j,i+1) +Nf_x(j+1,i) +Nf_x(j-1,i) ) /4.0 ;
                Nf_y(j,i)= (Nf_y(j,i-1) +Nf_y(j,i+1) +Nf_y(j+1,i) +Nf_y(j-1,i) ) /4.0 ;
                % normalize
                norm = sqrt(Nf_x(j,i)*Nf_x(j,i) + Nf_y(j,i)*Nf_y(j,i) ) ;
                 if (norm ~= 0 )
                     Nf_x(j,i) = Nf_x(j,i) / norm ;
                     Nf_y(j,i) = Nf_y(j,i) / norm ;
                 end
            end 
        end 
    end 
    
    % create onto the controur the value of the normal field.
    conversionx = tot_x_size/totx ;
    conversiony = tot_y_size/toty ;
    for k=1:totx*4
        alpha = -pi + k/totx*2*pi ;

        i =floor(totx/2) + floor(R_x/conversionx * cos(alpha )) ;
        j =floor(toty/2) + floor(R_y/conversiony * sin(alpha )) ;
        if (i == 0)
            i = 1 ;
        end
        if (j == 0)
            j = 1 ;
        end
	    Nf_x(j,i)=cos(alpha);                       
        Nf_y(j,i)=sin(alpha);
	
    end
        % create the field at the boundaries of the windows
    for k=2:totx-1
	    Nf_x(1,k)= 0.0;                       
        Nf_y(1,k)= -1.0 ; %-1
		
        Nf_x(toty,k)= 0.0 ;                       
        Nf_y(toty,k)= 1.0; %1
    end
    for k=2:toty-1
	    Nf_x(k,1)= -1.0 ;    %-1                   
        Nf_y(k,1)= 0 ;
		
        Nf_x(k,totx)= 1.0;                       
        Nf_y(k,totx)= 0 ;
    end
    Nf_x(1,1)= -1.0;                       
    Nf_y(1,1)= -1.0 ; 
		
    Nf_x(toty,totx)= 1.0 ;                       
    Nf_y(toty,totx)= 1.0; 
    
    Nf_x(1,totx)= 1.0 ;                       
    Nf_y(1,totx)= -1.0;
    
    Nf_x(toty,1)= -1.0 ;                       
    Nf_y(toty,1)= 1.0;  
    
	% norm the vectors:
    for i=1:totx
		for j=1:toty
            norm = sqrt(Nf_x(j,i)*Nf_x(j,i) + Nf_y(j,i)*Nf_y(j,i) ) ;
            if (norm ~= 0 )
                Nf_x(j,i) = Nf_x(j,i) / norm ;
                Nf_y(j,i) = Nf_y(j,i) / norm ;
            end
        end
    end
	filename = 'boule_r6_0um_idx.rid';
    fid = fopen (filename, 'w');
    fprintf(fid, 'UPI3DRI 3.0\n');
    fprintf(fid, '%d %d\n', totx, toty);
    fprintf(fid, '0 %f 0 %f\n', tot_x_size, tot_y_size);

    for j=1:toty
		for i=1:totx
    		fprintf(fid, '%f , %f\n', A(j,i), B(j,i));
		end
	end
	
	filename = 'NVF_x_interp.f3d';
    fid = fopen (filename, 'w');
    fprintf(fid, 'UPI3DRI 3.0\n');
    fprintf(fid, '%d %d\n', totx, toty);
    fprintf(fid, '0 %f 0 %f\n', tot_x_size, tot_y_size);

    for j=1:toty
		for i=1:totx
    		fprintf(fid, '%f,0\n', Nf_x(j,i));
		end
	end
    
    fclose (fid);
    
    filename = 'NVF_y_interp.f3d';
    fid = fopen (filename, 'w');
    fprintf(fid, 'UPI3DRI 3.0\n');
    fprintf(fid, '%d %d\n', totx, toty);
    fprintf(fid, '0 %f 0 %f\n', tot_x_size, tot_y_size);

    for j=1:toty
		for i=1:totx
    		fprintf(fid, '%f,0\n', Nf_y(j,i));
		end
	end
    
    fclose (fid);
    
    
    figure(1);
	subplot(2,1,1);
    pcolor(A);
    shading interp
	subplot(2,1,2);
    quiver(Nf_x,Nf_y);
    
    figure(2);
	subplot(3,1,1);
    pcolor(Nf_x);
    shading interp
	subplot(3,1,2);
    pcolor(Nf_y);
    shading interp
    subplot(3,1,3);
    pcolor(log(sqrt(Nf_y.*Nf_y +Nf_x.*Nf_x)) );
    shading interp
    
    
    figure(3);
	subplot(3,1,1);
    pcolor(Nf_x.*Nf_x);
    shading interp
	subplot(3,1,2);
    pcolor(Nf_x.*Nf_y);
    shading interp
 	subplot(3,1,3);
    pcolor(Nf_y.*Nf_y);
    shading interp
    
    figure(4);
    quiver(Nf_x,Nf_y);
end