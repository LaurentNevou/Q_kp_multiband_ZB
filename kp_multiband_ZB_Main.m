%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% last update 23Jan2019, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Here, you have to choose your material among the following %%%%%%%%%% 
%%%%%% Take care, not all the models are available for all the materials! %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Material='AlAs';
Material='GaAs';
%Material='InAs';
%Material='AlSb';
%Material='GaSb';
%Material='InSb';
%Material='AlP';
%Material='GaP';
%Material='InP';
%Material='C';
%Material='Si';
%Material='Ge';
%Material='Sn';
%Material='3C-SiC';
%Material='ZnSe';
%Material='ZnTe';
%Material='CdS';
%Material='ZnS';
%Material='CdSe';
%Material='CdTe';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=300;                  % Temperature [Kelvin]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Library
ExtractParameters


Nk=100;                   %% number of k-points for the dispersion
[k_list,k]=kZB_f(Nk,a);   %% function to compute the k-vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Activate the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kp_3x3_Luttinger   = 0;  % HH, LH
kp_4x4_Luttinger   = 0;  % HH, LH
kp_6x6_Luttinger   = 1;  % HH, LH, SO
kp_8x8_Kane        = 0;  % EC, HH, LH, SO (only EC is good, HH has even an opposite mass)
kp_8x8_Luttinger   = 1;  % EC, HH, LH, SO (The most used model)
kp_14x14_Luttinger = 0;  % GaAs, InAs, Si and Ge only
kp_16x16_Luttinger = 1;  % GaAs, InAs, Si and Ge only

plot_Kane_parabole=0;
plot_Luttinger_parabole=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=0;FS=20;
c=[
1 0 0
0 1 0
0 0 1
1 0 1
0 1 1
1 1 0
];

if kp_3x3_Luttinger == 1;
  i=i+1;
  E{i}=kp_3bands_DKK_f(k_list, g1, g2, g3);
  s{i}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{',num2str(c(i,:)),'}k.p 3x3: Luttinger');
end
if kp_4x4_Luttinger == 1;
  i=i+1;
  E{i}=kp_4bands_Luttinger_Fishman_f(k_list, g1, g2, g3);
  %E{i}=kp_4bands_Luttinger_Rosencher_f(k_list, g1, g2, g3);
  s{i}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{',num2str(c(i,:)),'}k.p 4x4: Luttinger');
end
if kp_6x6_Luttinger == 1;
  i=i+1;
  E{i}=kp_6bands_Luttinger_DKK_f(k_list, Dso, g1, g2, g3);
  %E{i}=kp_6bands_Luttinger_Kohn_f(k_list, Dso, g1, g2, g3);
  %E{i}=kp_6bands_Luttinger_Bastard_f(k_list, Dso, g1, g2, g3);
  %E{i}=kp_6bands_Luttinger_Cardona_f(k_list, Dso, g1, g2, g3);
  s{i}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{',num2str(c(i,:)),'}k.p 6x6: Luttinger');
end
if kp_8x8_Kane == 1;
  i=i+1;
  E{i}=kp_8bands_Kane_DKK_f(k_list, Eg, EP_K, Dso);
  %E{i}=kp_8bands_Kane_Bastard_f(k_list, Eg, EP_K, Dso);
  %E{i}=kp_8bands_Kane_Fishman_f(k_list, Eg, EP_K, Dso);
  s{i}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{',num2str(c(i,:)),'}k.p 8x8: Kane');
end
if kp_8x8_Luttinger == 1;
  i=i+1;
  E{i}=kp_8bands_Luttinger_DKK_f(k_list, Eg, EP_L, Dso, F, g1, g2, g3);
  %E{i}=kp_8bands_Luttinger_Pistol_f(k_list, Eg, EP_L, Dso, F, g1, g2, g3);
  %E{i}=kp_8bands_Luttinger_Fishman_f(k_list, Eg, EP_L, Dso, F, g1, g2, g3, gD1, gD2, gD3);
  s{i}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{',num2str(c(i,:)),'}k.p 8x8: Luttinger');
end
if kp_14x14_Luttinger == 1;
  i=i+1;
  if isnan(EP14(1))
    display('ERROR: Take care, the kp14x14 bands is only available for: GaAs, InAs, Si and Ge...')
    break  
  end
  E{i}=kp_14bands_Luttinger_Fishman_f(k_list, Eg8c, Eg7c, Eg, Dso, EP14, EPx14, F, gc123_14, g123);
  s{i}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{',num2str(c(i,:)),'}k.p 14x14: Luttinger');
end
if kp_16x16_Luttinger == 1;
  i=i+1;
  if isnan(EPs16)
    display('ERROR: Take care, the kp16x16 bands is only available for: GaAs, InAs, Si and Ge...')
    break  
  end
  E{i}=kp_16bands_Luttinger_Fishman_f(k_list, Eg8c, Eg7c, Eg, Dso, Eg6v, EP14, EPx14, EPs16, F, gc123_14, g123);
  s{i}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{',num2str(c(i,:)),'}k.p 16x16: Luttinger');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_Kane_parabole == 1;
  i=i+1;
  s{i}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{0 0 0}Kane parabole');
end
if plot_Luttinger_parabole == 1;
  i=i+1;
  s{i}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{0 0 0}Luttinger parabole');
end
if isempty(E)
  display('Choose at least one of the k.p model...')
  break
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[50 200 1000 900]);

FS=20;
LW=2;

xscale=[-5 5];
yscale=[-2 4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,1,1,'fontsize',FS)
hold on;grid on;

xlabel('Wavevector (nm-1)')
ylabel('Energy (eV)')
xlim(xscale)
ylim(yscale)

text(0.7*(xscale(2)-xscale(1))+xscale(1),0.5,s);
text(0*(xscale(2)-xscale(1))+xscale(1),-0.05*(yscale(2)-yscale(1))+yscale(1),'[111]','fontsize',FS);
text(0.95*(xscale(2)-xscale(1))+xscale(1),-0.05*(yscale(2)-yscale(1))+yscale(1),'[100]','fontsize',FS);
title(strcat(M{1},' bandstructure @T=',num2str(T),'K'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:length(E)
  
  plot(k*1e-9,E{j},'g-', 'linewidth',LW,'color',c(j,:))
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_Kane_parabole==1
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  h=6.62606896E-34;               %% Planck constant [J.s]
  hbar=h/(2*pi);
  e=1.602176487E-19;              %% electron charge [Coulomb]
  m0=9.10938188E-31;              %% electron mass [kg]
  H0=hbar^2/(2*m0) ;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  mc_K = 1 / (  1+ EP_K/3*(2/Eg + 1/(Eg+Dso)));
  ml_K  =-1 / (1-2*EP_K/(3*Eg));
  mso_K =-1 / (1-EP_K/(3*(Eg+Dso)));  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  f = 1/e*H0*k(1:Nk).^2;
  plot(k(1:Nk)*1e-9 ,  f/mc_K+Eg   ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(1:Nk)*1e-9 , -f/ml_K      ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(1:Nk)*1e-9 , -f/mso_K-Dso ,'linewidth',LW,'color','k','linestyle','--')
  
  f = 1/e*H0*k(Nk:end).^2;
  plot(k(Nk:end)*1e-9 ,  f/mc_K+Eg   ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(Nk:end)*1e-9 , -f/ml_K      ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(Nk:end)*1e-9 , -f/mso_K-Dso ,'linewidth',LW,'color','k','linestyle','--')
  
%  f = 1/e*H0*(k).^2;
%  plot(k*1e-9 ,  f/mc_K+Eg   ,'linewidth',LW,'color','k','linestyle','--')
%  plot(k*1e-9 , -f/ml_K      ,'linewidth',LW,'color','k','linestyle','--')
%  plot(k*1e-9 , -f/mso_K-Dso ,'linewidth',LW,'color','k','linestyle','--')

end

if plot_Luttinger_parabole==1
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  h=6.62606896E-34;               %% Planck constant [J.s]
  hbar=h/(2*pi);
  e=1.602176487E-19;              %% electron charge [Coulomb]
  m0=9.10938188E-31;              %% electron mass [kg]
  H0=hbar^2/(2*m0) ;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  mc_L = 1 / (  1+2*F + EP_L*(Eg+2*Dso/3) / (Eg*(Eg+Dso)) );
  
  mh100 = 1 / (g1-2*g2);
  mh110 = 1 / (g1-sqrt(g2^2+3*g3^2));
  mh111 = 1 / (g1-2*g3);
  
  ml100 = 1 / (g1+2*g2);
  ml110 = 1 / (g1+sqrt(g2^2+3*g3^2));
  ml111 = 1 / (g1+2*g3);
  
  mso_L = 1 / (g1-(EP_L*Dso/(3*Eg*(Eg+Dso))) );
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  f = 1/e*H0*k(1:Nk).^2;
  plot(k(1:Nk)*1e-9 ,  f/mc_L+Eg   ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(1:Nk)*1e-9 , -f/mh111     ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(1:Nk)*1e-9 , -f/ml111     ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(1:Nk)*1e-9 , -f/mso_L-Dso ,'linewidth',LW,'color','k','linestyle','--')
  
  f = 1/e*H0*k(Nk:end).^2;
  plot(k(Nk:end)*1e-9 ,  f/mc_L+Eg   ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(Nk:end)*1e-9 , -f/mh100     ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(Nk:end)*1e-9 , -f/ml100     ,'linewidth',LW,'color','k','linestyle','--')
  plot(k(Nk:end)*1e-9 , -f/mso_L-Dso ,'linewidth',LW,'color','k','linestyle','--')
  
  %f = 1/e*H0*k.^2;
  %plot(k*1e-9 ,  f/mc_L+Eg   ,'linewidth',LW,'color','k','linestyle','--')
  %plot(k*1e-9 , -f/mh110     ,'linewidth',LW,'color','k','linestyle','--')
  %plot(k*1e-9 , -f/ml110     ,'linewidth',LW,'color','k','linestyle','--')
  %plot(k*1e-9 , -f/mso_L-Dso ,'linewidth',LW,'color','k','linestyle','--')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%