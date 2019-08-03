%% Load_And_Process_NewPico_Data
%Author: Mohammad Tanhaemami
%Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.

% This function collects all of the New Pico data into a single file for each species.  The end result will be a structure where:
% DATA.Name will be the species.
% DATA.FeatureNames will be the list of features.
% DATA.Time(k) will be the time in days since the first time point.
% DATA.Stained(k).DATA will be the data at the kth time point.
clear; clc;
%Dates at which experiments are performed.
dates = {'072716','081516','072816','081616','080116','081716',...
    '080216','081816','080316','081916','080516','082216',...
    '080616','082316','080716','082616','080816','083016',...
    '081016','090216','081116','091116','081216'}; %Deelta_T's are the same between the two experiment batches??
Alphabet=char('A'+(1:7)-1)';
prefix_un = {''};
for i=1:length(Alphabet)-3
    for j=1:12
        prefix_un{i,j} = [Alphabet(i),num2str(j,'%02i')];
    end
end

prefix_st = {''};
for i=1:length(Alphabet)-4
    for j=1:12
        prefix_st{i,j} = [Alphabet(i+4),num2str(j,'%02i')];
    end
end
rel_dates = [0 19 1 20 5 21 6 22 7 23 9 26 10 27 11 30 12 34 14 37 15 46 16];
%Initialize a structure named "DATA".
DATA.Max = -inf*ones(1,14); %Initial maximum.
DATA.Min =  inf*ones(1,14); %Initial minimum.
for Species = {'Pico'} %Initialize species names in DATA.
    k=1; %Initialize counter.
    DATA.Name = Species{1}; %Assign "Nanno" as the name for DATA.
    for Date = dates %For each of the dates experiments are run.
        for Un_Prefix = prefix_un
            for m=1:size(Un_Prefix,1)
                File_unstained = [Un_Prefix{m,1},' ',Date{1},' ',Species{1},' un.fcs.csv']; %Save unstained data in a ".csv" file.
                try
                    X = importdata(File_unstained); %Import saved file for unstained as "X".
                    DATA.FeatureNames = X.colheaders; %Assgin header of each column of X to DATA.FeatureNames.
                    DATA.Unstained(k,m).DATA = X.data; %Assign data stored in X to DATA.Unstained(k).DATA.
                    DATA.Max = max(DATA.Max,max(X.data)); %Compare maximum of DATA and the data stored in X and assign the bigger value to DATA.Max.
                    DATA.Min = min(DATA.Min,min(X.data)); %Compare minimum of DATA and the data stored in X and assign the smaller value to DATA.Min.
                    DATA.Time(k) = rel_dates(k); %DATA.Time elements are the elements of the relative dates.
                catch
                end
            end
        end
        for St_Prefix = prefix_st
            for m=1:size(St_Prefix,1)
                File_stained = [St_Prefix{m,1},' ',Date{1},' ',Species{1},' st.fcs.csv']; %Save stained data in a ".csv" file.
                try
                    X = importdata(File_stained); %Import saved file for stained as "X".
                    DATA.FeatureNames = X.colheaders; %Assgin header of each column of X to DATA.FeatureNames.
                    DATA.Stained(k,m).DATA = X.data; %Assign data stored in X to DATA.Unstained(k).DATA.
                    DATA.Max = max(DATA.Max,max(X.data)); %Compare maximum of DATA and the data stored in X and assign the bigger value to DATA.Max.
                    DATA.Min = min(DATA.Min,min(X.data)); %Compare minimum of DATA and the data stored in X and assign the smaller value to DATA.Min.
                    DATA.Time(k) = rel_dates(k); %Compare maximum of DATA and the data stored in X and assign the bigger value to DATA.Max.
                catch
                end
            end
        end
        k=k+1; %Increase counter by one.
    end
    % Remove zeros (missed measurements):
    for i = 1:size(DATA.Stained,1)
        for j = 1:size(DATA.Stained,2)
            DATA.Stained(i,j).DATA(any(DATA.Stained(i,j).DATA==0,2),:) = [];
        end
        for j = 1:size(DATA.Unstained,2)
            DATA.Unstained(i,j).DATA(any(DATA.Unstained(i,j).DATA==0,2),:) = [];
        end
    end
    save(['New',Species{1},'_Data_NoZero.mat'],'DATA'); %Saves DATA info as *.mat file
end
