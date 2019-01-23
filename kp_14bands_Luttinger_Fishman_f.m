function[E]=kp_14bands_Luttinger_Fishman_f(k_list, Eg8c, Eg7c, Eg6c, Dso, EP, EP14, F, gc123, g123)

% Guy Fishman
% "Semi-Conducteurs: les Bases de la Theorie k.p " (2010)
% Structure de bande: II, 4.B Appendice : La matrice 14x14 ou hamiltonien H14, page 274
% https://www.amazon.fr/Semi-Conducteurs-Bases-Theorie-K-P-Fishman/dp/2730214976/ref=sr_1_fkmr1_1?ie=UTF8&qid=1548234034&sr=8-1-fkmr1&keywords=guy+fishman+kp
% https://www.abebooks.fr/semi-conducteurs-bases-th%C3%A9orie-k.p-Fishman-ECOLE/30091636895/bd
% https://www.decitre.fr/livres/semi-conducteurs-les-bases-de-la-theorie-k-p-9782730214971.html
% https://www.unitheque.com/Livre/ecole_polytechnique/Semi_conducteurs_les_bases_de_la_theorie_K.p-35055.html
% https://www.eyrolles.com/Sciences/Livre/semi-conducteurs-les-bases-de-la-theorie-k-p-9782730214971/

% Moustafa El Kurdi
% PhD thesis: "Dispositifs a ilots de GeSi pour la microphotonique proche infrarouge sur silicium"
% Chapter 1, page 31
% https://www.theses.fr/2003PA112148

% Said Ridene et al., PRB 64, 085329 (2001)
% "Infrared absorption in Si/Si(1-x)Ge(x)/Si quantum wells"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
e=1.602176487E-19;              %% charge de l electron [Coulomb]
m0=9.10938188E-31;              %% electron mass [kg]
H0=hbar^2/(2*m0) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eg8c  =  Eg8c*e;
Eg7c  =  Eg7c*e;
Eg6c  =  Eg6c*e;
Eg8v  =  0;
Eg7v  = -Dso*e;
%Dso  =  Dso*e;

EP  = EP *e;
EPx = EP14(1)*e;
EPp = EP14(2)*e;

P =  sqrt(EP *hbar^2/(2*m0)) ;
Px=  sqrt(EPx*hbar^2/(2*m0)) ;
Pp=  sqrt(EPp*hbar^2/(2*m0)) ;

%gc= 1+2*F + EP*(Eg6c-2*Eg7v/3) / (Eg6c*(Eg6c-Eg7v)) ;   % =1/mc  electron in CB eff mass

% renormalization of all the gamma for the 8bands Hamiltonian
g1=g123(1); g2=g123(2); g3=g123(3);
%gc = gc - EP/3*( 2/Eg6c + 1/(Eg6c+Dso) ) + EPp/3*( 2/(Eg8c-Eg6c) + 1/(Eg7c-Eg6c) );
gc = 1+2*F + EPp/3*( 2/(Eg8c-Eg6c) + 1/(Eg7c-Eg6c) );
g1 = g1 - EP/(3*Eg6c) - EPx/3*(1/(Eg7c-Eg8v) + 1/(Eg8c-Eg8v) )  ;
g2 = g2 - EP/(6*Eg6c) + EPx/6/(Eg7c-Eg8v);
g3 = g3 - EP/(6*Eg6c) - EPx/6/(Eg7c-Eg8v);

gD1=g1;gD2=g2;gD3=g3;
% gD1 = g1 - EP/3/Eg6c - EPx/3*(1/Eg7c + 1/Eg8c - 2/(Eg8c-Eg7v)) ;
% gD2 = g2 - EP/6/Eg6c - EPx/12*(1/Eg8c + 1/(Eg8c-Eg7v) -2/Eg7c) ;
% gD3 = g3 - EP/6/Eg6c + EPx/12*(1/Eg8c + 1/(Eg8c-Eg7v) -2/Eg7c) ;

gc1=gc123(1); gc2=gc123(2); gc3=gc123(3);
%gc123=[0 0 0]; %% flat band
%gc1 = gc123(1) + EPp/(3*(Eg8c-Eg6c)) + EPx/3*(1/(Eg8c-Eg8v) + 1/(Eg8c-Eg7v) )  ;
%gc2 = gc123(2) + EPp/(6*(Eg8c-Eg6c)) - EPx/6/(Eg8c-Eg7v);
%gc3 = gc123(3) + EPp/(6*(Eg8c-Eg6c)) + EPx/6/(Eg8c-Eg7v);

% if Silicon or Ge
%gc1=-1;
%gc2=0;
%gc3=0;

%gcD1 = gc1 + EPx/3 *(1/(Eg8c-Eg8v) + 1/(Eg8c-Eg7v) - 2/Eg7c       ) ; 
%gcD2 = gc2 + EPx/12*(1/(Eg8c-Eg8v) + 1/Eg7c        - 2/(Eg8c-Eg7v)) ; 
%gcD3 = gc3 - EPx/12*(1/(Eg8c-Eg8v) + 1/Eg7c        - 2/(Eg8c-Eg7v)) ; 

gcD1=gc1;gcD2=gc2;gcD3=gc3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%k+ =  kx + 1i*ky
%k- =  kx - 1i*ky

for i=1:length(k_list(:,1))
  
kx = k_list(i,1);
ky = k_list(i,2);
kz = k_list(i,3);

k=sqrt(kx.^2 + ky.^2 + kz.^2);
kpp =  kx + 1i*ky;
kmm =  kx - 1i*ky;

AA  = H0*g2 *( 2*kz^2 - kx.^2 - ky.^2 );
AAc = H0*gc2*( 2*kz^2 - kx.^2 - ky.^2 );
BB  = H0*2*sqrt(3)*g3 *kz*(kx - 1i*ky) ;
BBc = H0*2*sqrt(3)*gc3*kz*(kx - 1i*ky) ;
CC  = H0*sqrt(3)*(g2 *(kx^2-ky^2)-2i*g3 *kx*ky);
CCc = H0*sqrt(3)*(gc2*(kx^2-ky^2)-2i*gc3*kx*ky);

AA_D  = H0*gD2 *( 2*kz^2 - kx.^2 - ky.^2 );
AAc_D = H0*gcD2*( 2*kz^2 - kx.^2 - ky.^2 );

BB_D  = H0*2*sqrt(3)*gD3 *kz*(kx - 1i*ky) ;
BBc_D = H0*2*sqrt(3)*gcD3*kz*(kx - 1i*ky) ;

CC_D  = H0*sqrt(3)*(gD2 *(kx^2-ky^2)-2i*gD3 *kx*ky);
CCc_D = H0*sqrt(3)*(gcD2*(kx^2-ky^2)-2i*gcD3*kx*ky);

EH8c = Eg8c - gc1 *H0*k^2 + AAc;
EL8c = Eg8c - gc1 *H0*k^2 - AAc;
E7c  = Eg7c - gcD1*H0*k^2      ;
E6c  = Eg6c + gc  *H0*k^2      ;

EH8v =  0  - g1  *H0*k^2 + AA ;
EL8v =  0  - g1  *H0*k^2 - AA ;
E7v  = Eg7v- gD1 *H0*k^2      ;

Hdiag = [ EH8c EL8c EL8c EH8c E7c E7c E6c E6c EH8v EL8v EL8v EH8v E7v E7v];

% Ec  Ec      HH                  LH                LH               HH             SO               SO

H8x8=[
  0   0  -sqrt(1/2)*P*kpp   sqrt(2/3)*P*kz    sqrt(1/6)*P*kmm         0         sqrt(1/3)*P*kz    sqrt(1/3)*P*kmm  % Ec
  0   0        0           -sqrt(1/6)*P*kpp   sqrt(2/3)*P*kz   sqrt(1/2)*P*kmm  sqrt(1/3)*P*kpp  -sqrt(1/3)*P*kz   % Ec
  0   0        0                  BB                CC                0         sqrt(1/2)*BB_D    sqrt(2)  *CC_D   % HH
  0   0        0                   0                 0               CC        -sqrt(2)  *AA_D   -sqrt(3/2)*BB_D   % LH
  0   0        0                   0                 0              -BB        -sqrt(3/2)*BB_D'   sqrt(2)  *AA_D   % LH
  0   0        0                   0                 0                0        -sqrt(2)  *CC_D'   sqrt(1/2)*BB_D'  % HH
  0   0        0                   0                 0                0              0                0            % SO
  0   0        0                   0                 0                0              0                0            % SO
];

% HHc  LHc  LHc  HHc            SOc                SOc

H6x6=[
   0   BBc  CCc   0     sqrt(1/2)*BBc_D    sqrt(2)  *CCc_D  % HHc
   0    0    0   CCc   -sqrt(2)  *AAc_D   -sqrt(3/2)*BBc_D  % LHc
   0    0    0  -BBc   -sqrt(3/2)*BBc_D'   sqrt(2)  *AAc_D  % LHc
   0    0    0    0    -sqrt(2)  *CCc_D'   sqrt(1/2)*BBc_D' % HHc
   0    0    0    0              0                  0       % SOc
   0    0    0    0              0                  0       % SOc
];

Deltap=0;%Dso; % unknown parameter...

%        Ec                 Ec                HH              LH                 LH                HH              SO                 SO

H6x8=[
-sqrt(1/2)*Pp*kmm           0              1/3*Deltap  sqrt(1/3)*Px*kpp  sqrt(1/3)*Px*kz            0       sqrt(1/6)*Px*kpp   sqrt(2/3)*Px*kz  % HHc
 sqrt(2/3)*Pp*kz  -sqrt(1/6)*Pp*kmm -sqrt(1/3)*Px*kmm        1/3*Deltap           0        sqrt(1/3)*Px*kz          0        -sqrt(1/2)*Px*kpp  % LHc
 sqrt(1/6)*Pp*kpp  sqrt(2/3)*Pp*kz  -sqrt(1/3)*Px*kz           0               1/3*Deltap -sqrt(1/3)*Px*kpp sqrt(1/2)*Px*kmm           0        % LHc
          0        sqrt(1/2)*Pp*kpp           0       -sqrt(1/3)*Px*kz   sqrt(1/3)*Px*kmm        1/3*Deltap sqrt(2/3)*Px*kz  -sqrt(1/6)*Px*kmm  % HHc
 sqrt(1/3)*Pp*kz   sqrt(1/3)*Pp*kmm -sqrt(1/6)*Px*kmm          0        -sqrt(1/2)*Px*kpp -sqrt(2/3)*Px*kz       -2/3*Deltap           0        % SOc
 sqrt(1/3)*Pp*kpp -sqrt(1/3)*Pp*kz  -sqrt(2/3)*Px*kz   sqrt(1/2)*Px*kmm           0        sqrt(1/6)*Px*kpp          0             -2/3*Deltap  % SOc
];


H=[H6x6 H6x8 ; zeros(8,6) H8x8];

H14 = H' + H + diag(Hdiag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(H14)/e;

end

end