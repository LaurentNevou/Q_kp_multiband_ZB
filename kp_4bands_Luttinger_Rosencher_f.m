function[E]=kp_4bands_Luttinger_Rosencher_f(k_list, g1, g2, g3)

% Emmanuel Rosencher
% "Optoelectronics"
% Complement to Chapter 5C: "Kane’s k · p method", page 239
% https://www.cambridge.org/core/books/optoelectronics/86B6621671230A798D5BFBE24266EE3F
% https://www.amazon.fr/Optoelectronics-Emmanuel-Rosencher/dp/0521778131

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
e=1.602176487E-19;              %% charge de l electron [Coulomb]
m0=9.10938188E-31;              %% electron mass [kg]
H0=hbar^2/(2*m0) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%k+ = 1/sqrt(2)*( kx + 1i*ky )
%k- = 1/sqrt(2)*( kx - 1i*ky )

for i=1:length(k_list(:,1))
  
kx=k_list(i,1);
ky=k_list(i,2);
kz=k_list(i,3);

Hh = -H0*(g1+g2)*(kx^2+ky^2) -H0*(g1-2*g2)*kz^2;
Hl = -H0*(g1-g2)*(kx^2+ky^2) -H0*(g1+2*g2)*kz^2;
b  =  H0*2*sqrt(3)*g3*(kx-1i*ky)*kz;
c  =  H0*  sqrt(3)* (g2*(kx^2-ky^2)-2i*g3*kx*ky);

%   HH    LH    LH    HH
H=[
    Hh    c     b     0  % HH
    c'    Hl    0    -b  % LH
    b'    0     Hl    c  % LH
    0    -b'    c'    Hh % HH
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(H)/e;

end

end