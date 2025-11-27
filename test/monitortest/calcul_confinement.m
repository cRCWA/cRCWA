%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Grenoble INP - Phelma
%%            TD2 cours optique guidée PNS 3A - OM
%%
%%            Script pour Gnu Octave (probablement, il tourne aussi sous
%%            Matlab, peut-être en rajustant un poil certaines fonctions
%%            d'affichage).
%%
%%            Version 1.1, Davide Bucci, novembre 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DONNEES DU CALCUL

% Longueur d'onde en micron
lambda0 = 1.55

% Indices des couches du guide
% Core
n_c = 1.57
% Superstrat
n_sup =1.331
% Substrat
n_sub = 1.51
% Epaisseur en micron
d=1
% Ordre du mode à chercher
m=0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Calcul du vecteur d'onde
k0=2*pi/lambda0;


% Recherche de l'indice effectif par bisection
min = n_sub;
max = n_c;

toll = 1e-6;

for i=1:100
	n_eff=(max+min)/2;
	v=k0*d*sqrt(n_c^2-n_eff^2)-atan(sqrt(n_eff^2-n_sub^2)/sqrt(n_c^2-n_eff^2))-atan(sqrt(n_eff^2-n_sup^2)/sqrt(n_c^2-n_eff^2));
	
	if v>m*pi
		min = n_eff;
		
	else
		max = n_eff;
	end
	
	if max-min<toll 
		break;
	end
end

% Vérification si le mode est à la coupure
if n_eff-n_sub<toll
	printf("Mode non supporté par la structure !\n");
	break;
end

% Impression de l'indice effectif trouvé
printf("Indice effectif : %f\n", n_eff);

% Calcul de plusieurs paramètres pour la forme du mode
k_core = k0*sqrt(n_c^2-n_eff^2)
alpha_sup = k0*sqrt(n_eff^2-n_sup^2)
alpha_sub = k0*sqrt(n_eff^2-n_sub^2)

% Trois calculs différents pour psi
% Attention! Les trois résultats dépendent beaucoup de neff. Pour obtenir
% le même résultat, il faut déterminer neff avec beaucoup de précision.
psi = .5 * (atan(alpha_sup/k_core)-atan(alpha_sub/k_core))+m*pi/2
% psi = atan(alpha_sup/k_core)-d/2*k_core+m*pi
% psi = -atan(alpha_sub/k_core)+d/2*k_core

A=cos(k_core*d/2+psi)
B=cos(-k_core*d/2+psi)

% Dessin du profil modal. Au même temps, on en profite pour calculer
% numériquement les intégrales nécessaires au calcul du taux de confinement.
% Les intégrales sont calculées avec une méthode très basique, cela ne donne
% pas une grande précision, mais elle nous permet de vérifier nos calculs.
dy=.001;
y=[-10:dy:10];
phi=y;
norm=0;
norm_sub=0;
norm_core=0;
norm_sup=0;
tot_y=20;   % Largeur totale montrée (en micron).
ind=y;

for i=1:(tot_y/dy)+1
	if y(i)<-d/2    % Substrat
		phi(i)=B*exp(alpha_sub*(y(i)+d/2));
		norm_sub=norm_sub+phi(i)^2*dy;
		ind(i)=n_sub;
	elseif y(i)>d/2 % Superstrat
		phi(i)=A*exp(-alpha_sup*(y(i)-d/2));
		norm_sup=norm_sup+phi(i)^2*dy;
		ind(i)=n_sup;
	else            % Coeur du guide
		phi(i)=cos(k_core*y(i)+psi);
		norm_core=norm_core+phi(i)^2*dy;
		ind(i)=n_c;
	end	
	norm=norm+phi(i)^2*dy;
end

% Normes calculées numériquement
% Enlever le ";" pour que Octave les montre.
norm;
norm_sub;
norm_core;
norm_sup;

Gamma=norm_sup/norm;
printf("Calcul numérique, coeff. d'interaction : Gamma=%f %%\n", Gamma*100);


% Normes calculées analytiquement
norm_core=1/2*(cos(2*psi)*sin(k_core*d)+k_core*d)/k_core;
norm_sub=B^2/2/alpha_sub;
norm_sup=A^2/2/alpha_sup;

% Si tout se passe bien, norm_p ne doit pas être trop loin de norm
norm_p=norm_sub+norm_core+norm_sup;

% Représentation. Montre d'abord l'indice de réfraction
clf
subplot(2,1,1)
plot(y,ind);
xlabel ("y/µm")
ylabel ("Refractive index");
subplot(2,1,2)
% Montre la forme modale. On applique une normalisation, d'où le sqrt(norm_p)
plot(y,phi/sqrt(norm_p));
ylabel ("Modal field/a.u.")
xlabel ("y/µm")
% Calcul du coefficient d'intéraction 
Gamma=norm_sup/norm;
printf("Calcul analytique, coeff. d'interaction : Gamma=%f %%\n", Gamma*100);

% Attenuations
c1=1e-2
alpha_1=exp(-20*c1*1*Gamma)

10*log(alpha_1)/log(10)

c2=1e-1
alpha_2=exp(-20*c2*1*Gamma)

10*log(alpha_2)/log(10)

c3=1
alpha_3=exp(-20*c3*1*Gamma)
10*log(alpha_3)/log(10)
