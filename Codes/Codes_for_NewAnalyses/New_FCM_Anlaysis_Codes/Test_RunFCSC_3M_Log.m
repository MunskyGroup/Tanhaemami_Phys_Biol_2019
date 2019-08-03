%% Test_RunFCSC_3M_Log
% Author: Mohammad Tanhaemami
% Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.

% This script pulls the results of out fully developed weighted model and tests it on the new data set for FCM measurements on Pico sp.
% This script plot the results as measured and predicted values of the average lipid content following nitrogen starvation.
% THE GOAL here is to show that although our model is not trained based on the new data, it correctly sees the increase in lipid accumulations at each day.

%% INITIALIZATION, LOAD THE DESIRED DATABASE, AND ADD REQUIRED PATHS
clear; close all; clc;
warning('off');
load('../../../Data_Files/Data_NewAnalyses/FCM_Pico_Replica1_Data_NoZero.mat'); %Load the file Pico_Data.mat.
path(path,'./FCM_Functions'); %Path to Test_Codes.
if ~exist('./Results/Object_Results','dir')
    mkdir('./Results/Object_Results');
end
if ~exist('./Results/Figures','dir')
    mkdir('./Results/Figures');
end
if ~exist('./Results/Figures_FINAL','dir')
    mkdir('./Results/Figures_FINAL');
end

%% SECTION I: USER INPUTS & INITIAL VALUES
TargetSet = [3,9]; %Enter the target parmeters' index (here 3 and 9 are targets).
TargetInd = TargetSet(1,1); %Target the user looks for. Alternative: Target_Ind = Targets_Set(1,2);
TimePoint = 12; %Enter number of time points (days after nitorgen starvation) to be selected.
NCells = 3000;  %Enter number of cells to be used from all replication
NCellsL = 1000; %Enter number of cells to be used from each replication at each timepoint (labeled cells).
NCellsUL = 750; %Enter number of cells to be used from each replication at each timepoint (unlabeled cells).
bodipyBins = logspace(log10(DATA.Max(TargetInd))-4,log10(DATA.Max(TargetInd)),500); %Set logarithmic range for x-axis(500 x points).
bodipySignals = MakePlotBodipySignal_FCM(TargetInd,TimePoint,'True',1,bodipyBins,DATA);
% BINS = logspace(log10(DATA.Max(TargetInd))-4,log10(DATA.Max(TargetInd)),500); %Set logarithmic range for x-axis(500 x points).
BINS = linspace(log(DATA.Max(TargetInd))-10,log(DATA.Max(TargetInd))+40,500);
% BINS = linspace(log(DATA.Max(TargetInd))-40,log(DATA.Max(TargetInd))+40,500);
[SignalPeakInds_sort,SortInd] = sort(SignalPeakInd); % Sort the BODIPY signal peaks to get train/valid/test sets
TestSet = 1:TimePoint; % ALL of the data is used for testing
FeatureInds = 1:size(DATA.FeatureNames(~strcmp(DATA.FeatureNames,'Time')),2); %Define Feature_Inds
FeatureInds(TargetSet) = []; %Remove indicies for tragets.
FeatureInds_S = sort([FeatureInds,TargetSet]); %Define Feature_Inds for S test statistic

%% TESTING THE MODELS ON NEW DATA
% clear Alphas
close all;
for TestLRep = 3 % Replica chosen arbitrarily
    for TestULRep = 3 % Replica chosen arbitrarily
        close all;
        clear XTrainL XTrainUL YTrainL XValidL XValidUL YValidL
        [XTestL,YTestL,XTestUL,XTestUL_S] = DataCollector_Test_FCM(NCells,FeatureInds,FeatureInds_S,TargetInd,TestSet,TestLRep,TestULRep,DATA); %Extract data
        [XTestLQ,XTestULQ,XTestULQ_S] = deal({});
        for i=1:size(XTestL,1)
            XTestLQ{i,1} = Get_Quad_Features(XTestL{i,1}); %Getting quadratic Features.
            XTestULQ{i,1} = Get_Quad_Features(XTestUL{i,1}); %Getting quadratic Features.
            XTestULQ_S{i,1} = Get_Quad_Features(XTestUL_S{i,1}); %Getting quadratic Features.
        end
        
        %% TESTING STEP 6:
        %Call the model defined by the genetic algorithm in step 6 and extract the information regarding which features selected:
        close all;
        clear CurrentObject;
        Step = 6;
        load(['../../002SCENARIOS/Results_Scenario8_After2Trials/Object_Results/Step',num2str(Step),'_TRY_2.mat']); % Load the final weighted model
        XTrainLQ = CurrentObject.XTrainL; % Get the training data
        XTrainULQ = CurrentObject.XTrainUL; % Get the training data
        XTrainULQ_S = CurrentObject.XTrainUL_S; % Get the training data
        YTrainL = CurrentObject.YTrainL; % Get the training data
        M = CurrentObject.M; % Get the model
        Q = CurrentObject.Q; % Get the weight quotient
        IforS = CurrentObject.IforS; % Get the result of GA for useful columns of S (rows of Q)
        IvecBest1 = CurrentObject.I(:,1); % Get the results of GA for feature selection for Model 1
        IvecBest2 = CurrentObject.I(:,2); % Get the results of GA for feature selection for Model 2
        IvecBest3 = CurrentObject.I(:,3); % Get the results of GA for feature selection for Model 3
        %Do the testing for the model defined in step 6 by genetic algorithm:
        clear CurrentObject;
        if ~exist(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'file')
            CurrentObject = FCSC_3M_Log;
            CurrentObject.Step = Step;
            CurrentObject.Stage = 'Test';
            CurrentObject.BINS = BINS;
            CurrentObject.XTrainL = XTrainLQ;
            CurrentObject.XTrainUL = XTrainULQ;
            CurrentObject.XTrainUL_S = XTrainULQ_S;
            CurrentObject.YTrainL = YTrainL;
            CurrentObject.XTestL = XTestLQ;
            CurrentObject.YTestL = YTestL;
            CurrentObject.XTestUL = XTestULQ;
            CurrentObject.XTestUL_S = XTestULQ_S;
            CurrentObject.I(:,1) = IvecBest1';
            CurrentObject.I(:,2) = IvecBest2';
            CurrentObject.I(:,3) = IvecBest3';
            CurrentObject.M = M;
            CurrentObject.Q = Q;
            CurrentObject.IforS = IforS;
            CurrentObject.MakeHistPlot = 'False';
            CurrentObject.GenerateFigure = 'False';
            ModelResults = CurrentObject.RegressModel
            save(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'CurrentObject');
        else
            disp(['For the current replication, object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
            load(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat']);
        end
    end
end

%% Average lipid accumulation per day (mean +/- std.dev. vs time with shaded erorr bars)
Object_Test_Scen8 = CurrentObject;
LogYL_Scen8 = cell(TimePoint,1);
LogYPredUL_Scen8 = cell(TimePoint,1);
j = 1;
for i = 1:TimePoint % Testing time points (all of the time points for this new FCM_pico data set are considered as testing)
    LogYL_Scen8{i,1} = Object_Test_Scen8.RegressModel.LogYTestL{j,1};
    LogYPredUL_Scen8{i,1} = Object_Test_Scen8.RegressModel.LogYPredTestUL{j,1};
    j = j + 1;
end
% Compute means and standard deviations for measured and predicted results:
mu = []; stdv = []; muhat = []; stdvhat = [];
for i = 1:TimePoint
    mu = [mu;mean(LogYL_Scen8{i,1})];
    stdv = [stdv;std(LogYL_Scen8{i,1})];
    muhat = [muhat;mean(LogYPredUL_Scen8{i,1})];
    stdvhat = [stdvhat;std(LogYPredUL_Scen8{i,1})];
end
muhatScaled = muhat.*(mu(1)/muhat(1)); % Scale the predicted means to have the same mean at the initial time point
orient(figure(1),'landscape');
Days = DATA.Time; % The days that samples are taken
k = 1:TimePoint;
set(gca,'XTick',k);
set(gca,'XTickLabel',k);
hold on;
AvgLip_Measured = shadedErrorBar(k,mu,stdv,'lineprops',{'b-o','markerfacecolor','b'},'transparent',1); %Measured lipid content
AvgLip_Predicted = shadedErrorBar(k,muhatScaled,stdvhat,'lineprops',{'r-o','markerfacecolor','r'},'transparent',1); %Predicted lipid content
get(gca,'XTickLabel');
set(gca,'FontSize',20);
get(gca,'YTickLabel');
set(gca,'FontSize',20);
xlabel('Days After Nitrogen Starvation','FontSize',26);
ylabel('Average Lipid Content (AUC)','FontSize',26);
legend([AvgLip_Measured.edge AvgLip_Measured.mainLine AvgLip_Predicted.edge AvgLip_Predicted.mainLine],...
    {'+\sigma Measured','-\sigma Measured','\mu Measured',...
    '+\sigma Predicted','-\sigma Predicted','\mu Predicted'},...
    'Location','northwest','FontSize',20);
title('(g) Average Lipid Content for Labeled Replica 3 and Unlabeled Replica 3','FontSize',24);
saveas(figure(1),'Testing_weighted_model_on_FCMpico.fig');
set(gcf, 'PaperPosition', [0 0 15 15]);
set(gcf, 'PaperSize', [15 15]);
print(figure(1),'Testing_weighted_model_on_FCMpico','-dpdf','-fillpage');
