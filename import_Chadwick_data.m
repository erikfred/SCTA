% import_Chadwick_data.m
%
% Builds annual data structures for each lily tiltmeter and bottom pressure
% recorder at Axial, using Bill's direct data server rather than IRIS
% and/or OOI portal. (IRIS data was found to have database-induced gaps)
%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% BUILD  CAPABILITY TO APPEND INCOMING 2019 DATA %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


clear; close all

data_load=0; % 0 if building from scratch, 1 if appending to saved structure
endyr=str2double(datestr(now,'YYYY'));

datdir='../Chadwick_files/data/';
filelabels=['F';'E';'D';'B'];
stations=['AXCC1';'AXEC2';'AXID1';'ASHES'];

%% the code
if data_load==0
    if ~exist([datdir 'matfiles'],'dir')
        mkdir([datdir 'matfiles']) % where we'll be saving later on
    end
    
    for i=1:length(filelabels)
        for j=2014:endyr
            filestring=['MJ03' filelabels(i) num2str(j) 'NANO.Data'];
            yrstr=num2str(j); yrstr=yrstr(3:4);
            strucname=[stations(i,:) '_' yrstr];
            
            % read Bill's text file
            if exist([datdir 'matfiles/' strucname '.mat'],'file')
                disp(' ')
                disp(['matfile already exists for ' stations(i,:) ' in year ' num2str(j)])
                continue
            elseif exist([datdir filestring],'file')
                fid=fopen([datdir filestring],'r');
                
                % this particular file formatted differently
                if i==3 && j==2014
                    D=textscan(fid,'%s %s %f %f %f');
                    fclose(fid);
                    [temp{1:length(D{1})}]=deal(' ');
                    temp=temp';
                    temp2=strcat(D{1},temp,D{2});
                    D(:,1)=[];
                    D{:,1}=temp2;
                else
                    D=textscan(fid,'%s %f %f %f','delimiter',','); % {dates} {psi} {detided m} {C}
                    fclose(fid);
                end
            else
                disp(' ')
                disp(['No data for ' stations(i,:) ' in year ' num2str(j)])
                continue
            end
            
            % format to match other structures I use
            eval([strucname '.time=datenum(cell2mat(D{1}));'])
            eval([strucname '.p_psi=D{2};'])
            eval([strucname '.p_m=D{3};'])
            eval([strucname '.T=D{4};'])
            
            % save structure and clear it for next loop
            eval(['save(''' datdir 'matfiles/' strucname ''',''' strucname ''')'])
            eval(['clear ' strucname])
        end
    end
    
elseif data_load==1
    
end