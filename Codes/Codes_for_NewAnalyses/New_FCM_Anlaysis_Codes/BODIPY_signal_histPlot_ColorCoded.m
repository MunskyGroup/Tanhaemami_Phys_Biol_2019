%% Same as calling the function makePlotBodipySignal
%Author: Mohammad Tanhaemami
%Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.

% But with color coding with respect to time
% Just run this to see theresulting figure
path(path,'./FCM_Functions'); %Path to Test_Codes.
load('../../../Data_Files/Data_NewAnalyses/FCM_Pico_Replica1_Data_NoZero.mat'); %Load the file Pico_Data.mat.
TargetSet = [3,9]; %Enter the target parmeters' index (here 3 and 9 are targets).
TargetInd = TargetSet(1,1); %Target the user looks for. Alternative: Target_Ind = Targets_Set(1,2);
TimePoint = 12; %Enter number of time points (days after nitorgen starvation) to be selected.
BINS = logspace(log10(DATA.Max(TargetInd))-4,log10(DATA.Max(TargetInd)),500); %Set logarithmic range for x-axis(500 x points).

%%
figure(1); clf;
for k=1:TimePoint %For time points 1 to the final time point.
    try
        N = hist(DATA.Stained(k,1).DATA(:,TargetInd),BINS); %Create histogram of the stained target wrt x axis range(BINS).
        [Nmax(k),ind_max(k)] = max(N);
        semilogx(BINS,N,'Color',[0.4940+0.45-k/40,0.1840+0.45-k/40,0.5560+0.45-k/40]); %Plot data with logarithmic x axis.
        legendInfo{k} = ['TimePoint ',num2str(k)];
        hold on
    catch
    end
end
legend(legendInfo,'Location','bestoutside')
