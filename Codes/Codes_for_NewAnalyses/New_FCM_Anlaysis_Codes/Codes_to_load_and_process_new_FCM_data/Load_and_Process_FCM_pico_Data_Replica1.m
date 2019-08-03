%% Load_and_Process_FCM_pico_Data_Replica1 for the New Pico Measurements (Only for Replica 1)
% Author: Mohammad Tanhaemami
% Organization: Dr. Brian Munsky Research Lab, Department of Chemical and Biological Engineering, Colorado State University, fort Collins, CO, USA
% This script reads *.csv files for the new FCM measurements of Picochlorum sp.
% This script collects the data (ONLY FOR REPLICA 1 of the measurements), then processes and stores them in a structure format and as a *.mat file.
% This script then removes the "zero" values in the data set, because they are considered missed measurements
% The structure name is "DATA" and it is saved as "FCM_Pico_Replica1_Data_NoZero" in the current directory.
clear; clc;
Replica = 1; % For the other 2 replates this code should me modified according to the available *.csv files!!
dir = ['../../../../Data_Files/Data_NewAnalyses/New_FCM_Anlaysis_Codes/FCM_data/Replica',num2str(Replica),'/'];
%Dates at which experiments are performed:
dates = {'033019','040319','040819','041219','041919','042219','040119','040519','041019','041519','042419','041719'};
prefix_un = cell(0); % Pre-alocate prefixes for speed
Alphabet=char('A'+(1:7)-1)'; % Prepare the file names to be read
for i = 1:length(Alphabet)-3
    k = 1;
    for j = 1:3:10
        prefix_un{i,k} = [Alphabet(i),num2str(j,'%02i')]; % Generate beginning of the file names for unlabeled measurements
        k = k + 1;
    end
end
prefix_st = cell(0); % Pre-alocate prefixes for speed
for i = 5:length(Alphabet)
    k = 1;
    for j = 1:3:10
        prefix_st{i-4,k} = [Alphabet(i),num2str(j,'%02i')]; % Generate beginning of the file names for labeled measurements
        k = k + 1;
    end
end
[sorted_dates,rel_dates] = sort(dates);
%Initialize a structure named "DATA".
DATA.Max = -inf*ones(1,14); %Initial maximum.
DATA.Min =  inf*ones(1,14); %Initial minimum.
DATA.Name = {'Pico'};
k = 1;
for Date = sorted_dates
    for unPrefix = prefix_un
        for m = 1:size(unPrefix,1)
            fileName_un = [unPrefix{m,1},' ',Date{1},'_Pico_P_1_un.fcs.csv']; % Generate file names according to data files in the directory (unlabeled)
            try
                X = importdata([dir,fileName_un]); % Import data from the file
                DATA.FeatureNames = X.colheaders; % Get the features names
                DATA.Unstained(k,m).DATA = X.data; % Get the data at each time point
                DATA.Unstained(k,m).FileName = fileName_un; % Record the file name
                DATA.Max = max(DATA.Max,max(X.data)); % Get the maximum of the data
                DATA.Min = min(DATA.Min,min(X.data)); % Get the minimum of the data
                DATA.Time(k) = rel_dates(k); % Get measurement dates
            catch
            end
        end
    end
    for stPrefix = prefix_st
        for m = 1:size(stPrefix,1)
            fileName_st = [stPrefix{m,1},' ',Date{1},'_Pico_P_1_BOD.fcs.csv']; % Generate file names according to data files in the directory (labeled with BODIPY)
            try
                X = importdata([dir,fileName_st]); % Import data from the file
                DATA.FeatureNames = X.colheaders; % Get the features names
                DATA.Stained(k,m).DATA = X.data; % Get the data at each time point
                DATA.Stained(k,m).FileName = fileName_st; % Record the file name
                DATA.Max = max(DATA.Max,max(X.data)); % Get the maximum of the data
                DATA.Min = min(DATA.Min,min(X.data)); % Get the minimum of the data
                DATA.Time(k) = rel_dates(k); % Get measurement dates
            catch
            end
        end
    end
    k=k+1; %Increase counter by one.
end
% Remove zeros (missed measurements):
for i = 1:size(DATA.Stained,1)
    for j = 1:size(DATA.Stained,2)
        DATA.Stained(i,j).DATA(any(DATA.Stained(i,j).DATA==0,2),:) = []; % Search for zero values in each row and remove that row (if any)
    end
    for j = 1:size(DATA.Unstained,2)
        DATA.Unstained(i,j).DATA(any(DATA.Unstained(i,j).DATA==0,2),:) = []; % Search for zero values in each row and remove that row (if any)
    end
end
save(['FCM_Pico_Replica',num2str(Replica),'_Data_NoZero.mat'],'DATA'); %Saves DATA info as *.mat file









