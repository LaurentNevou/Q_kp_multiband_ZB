%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Library load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DB      = importdata('materialDB.csv'           ,',');
DB_kp16 = importdata('materialDB_kp16bands.csv' ,',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M{1}=Material;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Here is a patch to be able to load the tables in Matlab AND Octave %%%%%%
% Matlab see the header in multiple cells while Octave see the header in one cell only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M{1} ',']);
    idxM=find(idxM==idx);
    
    M{2} = DB.data(:,idxM);
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M{1});
      if idx==1
        M{2}=DB.data(:,i-1);
        break  
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(DB_kp16.textdata(1,:))==1  %% Octave data load

    DB_kp16.textdata{1,1}=[DB_kp16.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB_kp16.textdata{1,1},',');
    idx=strfind(DB_kp16.textdata{1,1},[',' M{1} ',']);
    if isempty(idx)
      M{3}=nan(8,1);
      break    
    end
    idxM=find(idxM==idx);
    
    M{3} = DB_kp16.data(:,idxM);
else  %% Matlab data load

    for i=1:length(DB_kp16.textdata(1,:))
      idx=strcmp(DB_kp16.textdata{1,i},M{1});
      if idx==1
        M{3}=DB_kp16.data(:,i-1);
        break  
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%