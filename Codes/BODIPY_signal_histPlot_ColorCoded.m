%% BODIPY_signal_histPlot_ColorCoded
%Author: Mohammad Tanhaemami
%Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.

% Same as calling the function makePlotBodipySignal, but with color coding with respect to time
path(path,'./Functions'); % Path to the functions folder
load('../Data_Files/NewPico_Data_NoZero.mat'); % Load the data file
TargetSet = [3,9]; % Enter the target parmeters' index (here 3 and 9 are targets).
TargetInd = TargetSet(1,1); % Target the user looks for -- here: target 3. Alternative: Target_Ind = Targets_Set(1,2);
TimePoint = 23; % Enter number of time points (days after nitorgen starvation) to be selected.
BINS = logspace(log10(DATA.Max(TargetInd))-4,log10(DATA.Max(TargetInd)),500); % Set logarithmic range for x-axis(500 x points).

%%
figure(1); clf;
for k=1:TimePoint
    try
        N{k} = hist(DATA.Stained(k,1).DATA(:,TargetInd),BINS); % Create histogram of the stained target wrt x axis range(BINS).
        [~,ind_max(k)] = max(N{k}); % Record maximum of hists in each bin for future analyses
        semilogx(BINS,N{k},'Color',[0.4940+0.45-k/40,0.1840+0.45-k/40,0.5560+0.45-k/40]); % Plot data with logarithmic x axis.
        legendInfo{k} = ['TimePoint ',num2str(k)]; % Legend
        hold on
    catch
    end
end
legend(legendInfo,'Location','bestoutside'); % Legend
