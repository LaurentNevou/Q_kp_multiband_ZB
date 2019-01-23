function[E]=kp_8bands_Kane_Bastard_f(k_list, Eg, EP, Dso)

% Gerard Bastard
% "Wave mechanics applied to semiconductor heterostructures" (1992)
% page 43

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

%k+ = 1/sqrt(2)*( kx + 1i*ky )
%k- = 1/sqrt(2)*( kx - 1i*ky )

for i=1:length(k_list(:,1))
  
kx=k_list(i,1);
ky=k_list(i,2);
kz=k_list(i,3);

k=sqrt(kx.^2 + ky.^2 + kz.^2);
kpp = sqrt(1/2)*( kx + 1i*ky );
kmm = sqrt(1/2)*( kx - 1i*ky );

%%%%%%%%%%%%%%%%%%%%%% Bastard filling method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=zeros(8,8);

Hdiag = [ Eg+H0*k^2  H0*k^2  H0*k^2  H0*k^2-Dso    Eg+H0*k^2  H0*k^2  H0*k^2  H0*k^2-Dso ];

%  Ec       LH             HH            SO              EC            LH            HH        SO

H=[
   0  -sqrt(2/3)*P*kz     P*kpp     sqrt(1/3)*P*kz       0        -sqrt(1/3)*P*kmm   0   -sqrt(2/3)*P*kmm  % EC
   0        0               0            0        sqrt(1/3)*P*kmm       0            0         0           % LH
   0        0               0            0               0              0            0         0           % HH
   0        0               0            0        sqrt(2/3)*P*kmm       0            0         0           % SO
   0        0               0            0               0        -sqrt(2/3)*P*kz  P*kmm  sqrt(1/3)*P*kz   % EC
   0        0               0            0               0              0            0         0           % LH
   0        0               0            0               0              0            0         0           % HH
   0        0               0            0               0              0            0         0           % SO
];

H=H'+H+diag(Hdiag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(H)/e;

end

end