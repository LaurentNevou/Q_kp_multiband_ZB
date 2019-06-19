%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Extract general parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lattice parameter

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'a');
  if idx==1
    a=M{2}(i-1)*1e-10;
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EP from the Kane model

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'EP_Kane');
  if idx==1
    EP_K = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EP from the Luttinger model

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'EP_Luttinger');
  if idx==1
    EP_L = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy level of the band 7c

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'Eg7c');
  if idx==1
    Eg7c = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy level of the band 8c

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'Eg8c');
  if idx==1
    Eg8c = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy level of the band 6v

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'Eg6v');
  if idx==1
    Eg6v = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy level of the band 6c (the band gap in a direct band gap semiconductor)

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'Eg6cG');
  if idx==1
    Eg6c = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bandgap Temperature dependency parameter "alpha" at the Gamma point

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'alphaG');
  if idx==1
    alphaG = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bandgap Temperature dependency parameter "beta" at the Gamma point

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'betaG');
  if idx==1
    betaG = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eg   = Eg6c - (alphaG*T^2) ./ (T+betaG);   %Eg = Eg0 - (a*T.^2)./(T + b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% split-off band energy Delta-so

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'Dso');
  if idx==1
    Dso = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter for the electrons

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'F');
  if idx==1
    F = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter gamma1 

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'gamma1');
  if idx==1
    g1 = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter gamma2

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'gamma2');
  if idx==1
    g2 = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter gamma3

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'gamma3');
  if idx==1
    g3 = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g123 = [g1 g2 g3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% here are the parameters for the 14x14 and 16x16  bands %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EP14

for i=1:length(DB_kp16.textdata(:,1))
  idx=strcmp(DB_kp16.textdata{i,1},'EP_14');
  if idx==1
    EP14 = M{3}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPx14

for i=1:length(DB_kp16.textdata(:,1))
  idx=strcmp(DB_kp16.textdata{i,1},'EPx_14');
  if idx==1
    EPx14 = [ M{3}(i-1) M{3}(i)];
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPs16

for i=1:length(DB_kp16.textdata(:,1))
  idx=strcmp(DB_kp16.textdata{i,1},'EPs_16');
  if idx==1
    EPs16 = [ M{3}(i-1) M{3}(i)];
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conduction band GammaC

for i=1:length(DB_kp16.textdata(:,1))
  idx=strcmp(DB_kp16.textdata{i,1},'gc1_14');
  if idx==1
    gc123_14 = [ M{3}(i-1) M{3}(i) M{3}(i+1)];
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%