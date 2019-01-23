function[E]=kp_6bands_Luttinger_Cardona_f(k_list, Dso, g1, g2, g3)

% Yu and Cardona
% "Fundamentals of Semiconductors"
% p78, 2.6/The k·p Method of Band-Structure Calculations
% Take care, there is a lot of mistake in the book of Cardona p78...
% Therefore, the values of H are sightly different from the book to make it work

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

MM=M;
LL=L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(k_list(:,1))
  
kx = k_list(i,1);
ky = k_list(i,2);
kz = k_list(i,3);

k=sqrt(kx.^2 + ky.^2 + kz.^2);

Hdiag(1) =  (H0+(2*M+L)/3)*k^2  -  (M-L)/6*(kx^2 + ky^2 - 2*kz^2);
Hdiag(2) =  H0*k^2 + 1/3*(M+2*L)*k^2 - 1/2*(L-M)*(kx^2+ky^2);
Hdiag(3) =  Hdiag(2);
Hdiag(4) =  Hdiag(1);
Hdiag(5) =  H0*k^2 + 1/3*(2*MM+LL)*k^2 - Dso;
Hdiag(6) =  Hdiag(5);

H=zeros(6,6);

H(1,2) =  sqrt(1/3)*N*(kx-1i*ky)*kz;
H(1,3) = -1/2*sqrt(1/3)*( (L-M)*(kx^2-ky^2)-2i*N*kx*ky   );
H(1,5) =  sqrt(1/2)*H(1,2);
H(1,6) =  sqrt(2)*H(1,3);

H(2,4) =  H(1,3);
H(2,5) =  sqrt(1/2)*(Hdiag(2)-Hdiag(1));
H(2,6) = -sqrt(3/2)*H(1,2);

H(3,4) = -H(1,2);
H(3,5) =  H(2,6)';
H(3,6) = -H(2,5);

H(4,5) = -sqrt(2)*H(1,3)';
H(4,6) =  H(1,5)';

H=H'+H+diag(Hdiag);

%     HH  LH  LH  HH  SO  SO
% HH   .   .   .   0   .   .
% LH   .   .   0   .   .   .
% LH   .   0   .   .   .   .
% HH   0   .   .   .   .   .
% SO   .   .   .   .   .   0
% SO   .   .   .   .   0   .

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(H)/e;

end

end