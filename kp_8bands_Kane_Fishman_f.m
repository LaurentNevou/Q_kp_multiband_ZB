function[E]=kp_8bands_Kane_Fishman_f(k_list, Eg, EP, Dso)

% Guy Fishman
% "Semi-Conducteurs: les Bases de la Theorie k.p " (2010)
% 3.3.2/ L’hamiltonien à l’intérieur de {G6 ; G8 ; G7} ou hamiltonien de Kane
% page 152
% https://www.amazon.fr/Semi-Conducteurs-Bases-Theorie-K-P-Fishman/dp/2730214976/ref=sr_1_fkmr1_1?ie=UTF8&qid=1548234034&sr=8-1-fkmr1&keywords=guy+fishman+kp
% https://www.abebooks.fr/semi-conducteurs-bases-th%C3%A9orie-k.p-Fishman-ECOLE/30091636895/bd
% https://www.decitre.fr/livres/semi-conducteurs-les-bases-de-la-theorie-k-p-9782730214971.html
% https://www.unitheque.com/Livre/ecole_polytechnique/Semi_conducteurs_les_bases_de_la_theorie_K.p-35055.html
% https://www.eyrolles.com/Sciences/Livre/semi-conducteurs-les-bases-de-la-theorie-k-p-9782730214971/

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

%k+ =  kx + 1i*ky 
%k- =  kx - 1i*ky 

for i=1:length(k_list(:,1))
  
kx=k_list(i,1);
ky=k_list(i,2);
kz=k_list(i,3);

k=sqrt(kx.^2 + ky.^2 + kz.^2);
kpp =  kx + 1i*ky;
kmm =  kx - 1i*ky;

%%%%%%%%%%%%%%%%%%%%%% Fishman filling method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=zeros(8,8);

Hdiag = H0*k^2*ones(1,8) + [ Eg Eg 0 0 0 0 -Dso -Dso ];

%  Ec      EC       HH                 LH                LH               HH             SO               SO

H=[
   0       0  -sqrt(1/2)*P*kpp   sqrt(2/3)*P*kz    sqrt(1/6)*P*kmm         0         sqrt(1/3)*P*kz    sqrt(1/3)*P*kmm % EC
   0       0        0           -sqrt(1/6)*P*kpp   sqrt(2/3)*P*kz   sqrt(1/2)*P*kmm  sqrt(1/3)*P*kpp  -sqrt(1/3)*P*kz  % EC
   0       0        0                   0                 0                0              0                0           % HH
   0       0        0                   0                 0                0              0                0           % LH
   0       0        0                   0                 0                0              0                0           % LH
   0       0        0                   0                 0                0              0                0           % HH
   0       0        0                   0                 0                0              0                0           % SO
   0       0        0                   0                 0                0              0                0           % SO
];

H=H'+H+diag(Hdiag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(H)/e; 

end

end