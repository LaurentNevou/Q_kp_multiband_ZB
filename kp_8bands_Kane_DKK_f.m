function[E]=kp_8bands_Kane_DKK_f(k_list, Eg, EP, Dso)

% DKK model: Dresselhaus, Kip and Kittel

% Calin Galeriu
% PhD thesis: "k.p Theory of semiconductor nanostructures" (2005)
% Chapter 3, page 26
% Download:
% https://web.wpi.edu/Pubs/ETD/Available/etd-120905-095359/unrestricted/cgaleriu.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
e=1.602176487E-19;              %% electron charge [Coulomb]
m0=9.10938188E-31;              %% electron mass [kg]
H0=hbar^2/(2*m0) ;

Dso = Dso*e;
Eg  = Eg*e;
EP  = EP*e;
P   = sqrt(EP*hbar^2/(2*m0)) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(k_list(:,1))

kx = k_list(i,1);
ky = k_list(i,2);
kz = k_list(i,3);

k=sqrt(kx.^2 + ky.^2 + kz.^2);
Hdiag = H0*k^2  +  [Eg -Dso/3 -Dso/3 -Dso/3 ];

H4=[

  0     1i*P*kx   1i*P*ky   1i*P*kz
  0        0         0         0
  0        0         0         0
  0        0         0         0

];

HH4 = H4' + H4 + diag(Hdiag); 

Hso=[

 0   0   0   0   0   0   0   0
 0   0   1   0   0   0   0   1i
 0  -1   0   0   0   0   0   1
 0   0   0   0   0  -1i -1   0
 0   0   0   0   0   0   0   0
 0   0   0  -1i  0   0  -1   0
 0   0   0   1   0   1   0   0
 0   1i -1   0   0   0   0   0

];

H=[HH4  zeros(4,4) ; zeros(4,4)  HH4] + Hso*Dso/(3i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(H)/e ;

end

end