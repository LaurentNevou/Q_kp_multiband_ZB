function[E]=kp_6bands_Luttinger_DKK_f(k_list, Dso, g1, g2, g3)

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
e=1.602176487E-19;              %% charge de l electron [Coulomb]
m0=9.10938188E-31;              %% electron mass [kg]
H0=hbar^2/(2*m0) ;
Dso=Dso*e;

L =  H0*(-1-g1-4*g2);
M =  H0*(-1-g1+2*g2);
N = -H0*6*g3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(k_list(:,1))
  
kx = k_list(i,1);
ky = k_list(i,2);
kz = k_list(i,3);

k=sqrt(kx.^2 + ky.^2 + kz.^2);

Hdiag = H0*k^2  - Dso/3*ones(1,6);

H_DKK=[

L*kx^2+M*(ky^2+kz^2)           N*kx*ky                 N*kx*kz
      N*kx*ky           L*ky^2+M*(kx^2+kz^2)           N*ky*kz
      N*kx*kz                  N*ky*kz           L*kz^2+M*(kx^2+ky^2)
];

Hso=[

 0   1   0   0   0   1i
-1   0   0   0   0   1
 0   0   0  -1i -1   0
 0   0  -1i  0  -1   0
 0   0   1   1   0   0
 1i -1   0   0   0   0

];

H = diag(Hdiag) + [H_DKK  zeros(3,3) ; zeros(3,3)  H_DKK] + Hso*Dso/(3i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(H)/e ;

end

end