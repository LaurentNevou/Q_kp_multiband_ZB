function[E]=kp_6bands_Luttinger_Bastard_f(k_list, Dso, g1, g2, g3)

% Gerard Bastard
% "Wave mechanics applied to semiconductor heterostructures" (1992)
% page 52

% Bastard shows a 8bands matrix but the CB is not couple to the other bands, p52
% Therefore, I just take into account the 6bands
% It should be the Hamiltonian of Bir and Pikus
% Moreover, I saw some mistakes in 
% H = -D*kz*( kx - 1i*ky) instead of 
% H = -D*kz*( kx + 1i*ky)
% I = -sqrt(3)/2*B*(kx^2-ky^2) + 1i*D*kx*ky  instead of 
% I = +sqrt(3)/2*B*(kx^2-ky^2) - 1i*D*kx*ky

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

A = (L+2*M)/3;
B = (L-M)/3;
D = N/sqrt(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(k_list(:,1))
  
kx = k_list(i,1);
ky = k_list(i,2);
kz = k_list(i,3);

k=sqrt(kx.^2 + ky.^2 + kz.^2);

F = A*k^2+B/2*(k^2-3*kz^2);
G = A*k^2-B/2*(k^2-3*kz^2);
H = -D*kz*(kx-1i*ky);                    % Here, there is a mistake in the book
I = -sqrt(3)/2*B*(kx^2-ky^2)+1i*D*kx*ky; % Here, there is a mistake in the book

Hdiag = H0*k^2*[1 1 1] + [ G  F  1/2*(F+G)-Dso ];
Hdiag = [Hdiag Hdiag];

%  LH    HH       SO              LH              HH            SO

%HH=[
%    0    -H'  sqrt(1/2)*(G-F)     0              I       sqrt(3/2)*H    % LH
%    0     0   -sqrt(1/2)*H        I              0       sqrt(2)*I      % HH
%    0     0       0           sqrt(3/2)*H   -sqrt(2)*I       0          % SO
%    0     0       0               0              H     -sqrt(1/2)*(G-F) % LH
%    0     0       0               0              0     -sqrt(1/2)*H'    % HH
%    0     0       0               0              0           0          % SO
%];

HH=[
    0   1i*H'  -sqrt(1/2)*(G-F)     0              I      1i*sqrt(3/2)*H    % LH
    0     0     1i*sqrt(1/2)*H     -I              0       -sqrt(2)*I       % HH
    0     0        0           -1i*sqrt(3/2)*H  sqrt(2)*I       0           % SO
    0     0        0                0           -1i*H      -sqrt(1/2)*(G-F) % LH
    0     0        0                0              0       -1i*sqrt(1/2)*H' % HH
    0     0        0                0              0            0           % SO
];

HH=HH'+HH+diag(Hdiag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(HH)/e;

end

end