%% RunFCSC_Log
%Author: Mohammad Tanhaemami
%Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University. 

% This script generated a model based on averaged statistics of the data as opposed to the weighted modeling strategy.
% The analysis procedure is carried our in several steps.

%% INITIALIZATION, LOAD THE DESIRED DATABASE, AND ADD REQUIRED PATHS
clear; close all; clc;
warning('off');
load('../../Data_Files/NewPico_Data_NoZero.mat'); %Load the file NewPico_Data_NoZero.mat.
path(path,'../Functions'); %Path to funcitions.
% The resuls will be saved in the following directories:
if ~exist('./Results/Object_Results','dir') % The results saved as matlab objects
    mkdir('./Results/Object_Results');
end
if ~exist('./Results/Figures','dir') % Temporary directory for figures
    mkdir('./Results/Figures');
end
if ~exist('./Results/Figures_FINAL','dir') % The results saved as figures for visualization
    mkdir('./Results/Figures_FINAL');
end

%% SECTION I: USER INPUTS & INITIAL VALUES
TargetSet = [3,9]; %Enter the target parmeters' index (here 3 and 9 are targets).
TargetInd = TargetSet(1,1); %Target the user looks for. Alternative: Target_Ind = Targets_Set(1,2);
TimePoint = 23; %Enter number of time points (days after nitorgen starvation) to be selected.
NCells = 3000;  %Enter number of cells to be used from all replication
NCellsL = 1000; %Enter number of cells to be used from each replication at each timepoint (labeled cells).
NCellsUL = 750; %Enter number of cells to be used from each replication at each timepoint (unlabeled cells).
% BINS = logspace(log10(DATA.Max(TargetInd))-4,log10(DATA.Max(TargetInd)),500); %Set logarithmic range for x-axis(500 x points).
BINS = linspace(log(DATA.Max(TargetInd))-14,log(DATA.Max(TargetInd)),500); %%%%%%%%% Scenario 2 (regular approach)
% BINS = linspace(log(DATA.Max(TargetInd))-10,log(DATA.Max(TargetInd))+40,500); %%%%%%%%% From Scenario 8 (weighted modeling)
TrainSet = [2,10,23]; %[SortInd(1) SortInd(end)]; %Widest range in Bodipy_Signals.
ValidSet = [1,11,22]; %[SortInd(2) SortInd(end-1)]; %Second widest range in Bodipy_Signals.
TestSet = [3:9,12:21]; %sort(SortInd(3:end-2)); %Define testing dataset (not used before).
FeatureInds = 1:size(DATA.FeatureNames(~strcmp(DATA.FeatureNames,'Time')),2); %Define Feature_Inds
FeatureInds(TargetSet) = []; %Remove indicies for tragets.

%% SECTION II: COLLECT DATA FOR ALL COMBINATIONS OF REPLICATIONS IN TRAINING & VALIDATION DATA SETS
% Collect the data for further analysis
[XTrainL,XTrainUL,YTrainL,XValidL,XValidUL,YValidL] = DataCollector_Concat_TrainValid(NCellsL,NCellsUL,FeatureInds,TargetInd,TrainSet,ValidSet,DATA);

%% STEP 1: REGRESSION - PREDICT LIPID CONTENT WITH ALL FEATURES
% At this step, linear regression is used to predict the lipid content using all available features
clear CurrentObject;
Step = 1;
if ~exist(['Results/Object_Results/Step',num2str(Step),'.mat'],'file')
    CurrentObject = FCSC_Log; % Initialize the matlab object
    CurrentObject.Step = Step; % Set the step of the calculations (here stpe 1 for simple linear regression)
    CurrentObject.Stage = 'TrainValid'; % Set the stage for calculations (here is the training and validation stage)
    CurrentObject.BINS = BINS; % Assign bins
    CurrentObject.XTrainL = XTrainL; % Assign labeled training features
    CurrentObject.XTrainUL = XTrainUL; % Assing unlabeled training features
    CurrentObject.YTrainL = YTrainL; % Assign the measured labeled trainsing taragets
    CurrentObject.XTestL = XValidL; % Assign labeled validation features
    CurrentObject.XTestUL = XValidUL; % Assing unlabled validation features
    CurrentObject.YTestL = YValidL; % Assign the measured labeled validation targets
    CurrentObject.I = ones(size(XValidL{1,1},2),1); % Assing the vector of selected features (here is all features)
    CurrentObject.MakeHistPlot = 'True'; % Logical set to make the histogram plots
    CurrentObject.GenerateFigure = 'True'; % Logical set to generate figures for results
    ModelResults = CurrentObject.RegressModel
    save(['Results/Object_Results/Step',num2str(Step),'.mat'],'CurrentObject');
    set(gcf, 'PaperPositionMode', 'auto')
    print(figure(18),['Results/Figures_FINAL/Step',num2str(Step),'_Sep_FINAL'],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
    print(figure(22),['Results/Figures_FINAL/Step',num2str(Step),'_Comb_FINAL'],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
else
    disp(['Object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
    load(['Results/Object_Results/Step',num2str(Step),'.mat']);
end

%% STEP 2: FIND POSSIBLE CORRUPTIONS -- THIS STEP WILL BE SKIPPED FOR NOW
% SelectedTimePoint = 6; %We have 13 time points. We select one of them, arbitrarily.
% Step = 2;
% if ~exist(['Results/Figures_FINAL/Step2_TrainRep',num2str(TrainRep),'_ValidLRep',num2str(ValidLRep),'_ValidULRep',num2str(ValidULRep),'_FINAL.pdf'],'file') %Check whether figure exists
%     [Stained,Unstained,KS] = Corruption_Test(Step,NCells,FeatureInds,SelectedTimePoint,ValidLRep,ValidULRep,DATA);
%     saveas(figure(SelectedTimePoint),['Results/Figures_FINAL/Step2_TrainRep',num2str(TrainRep),'_ValidLRep',num2str(ValidLRep),'_ValidULRep',num2str(ValidULRep),'_FINAL.pdf']); %Saving final figure as a *.pdf
% else
%     disp('For the current replication combination, we already know which features are corrupted. Skipping Step 2.')
% end

%% STEP 3: REGRESSION REDUCED - PREDICT LIPID CONTENT WITH REDUCED FEATURES
% After removing the corrupted features (here features 4 and 10), linear regression is used again at this step to predict the lipid content
CorruptedFeatures = [4,10]; % User input for corrupted features to be removed from analysis
for i=1:length(CorruptedFeatures)
    FeatureInds(FeatureInds==CorruptedFeatures(i)) = [];
end
%Redo data extraction after introducing new feature set (without corrupted features)
clear XTrainL YTrainL XValidL YValidL XValidUL
[XTrainL,XTrainUL,YTrainL,XValidL,XValidUL,YValidL] = DataCollector_Concat_TrainValid(NCellsL,NCellsUL,FeatureInds,TargetInd,TrainSet,ValidSet,DATA);
close all; clear CurrentObject;
Step = 3;
if ~exist(['Results/Object_Results/Step',num2str(Step),'.mat'],'file')
    CurrentObject = FCSC_Log;
    CurrentObject.Step = Step;
    CurrentObject.Stage = 'TrainValid';
    CurrentObject.BINS = BINS;
    CurrentObject.XTrainL = XTrainL;
    CurrentObject.XTrainUL = XTrainUL;
    CurrentObject.YTrainL = YTrainL;
    CurrentObject.XTestL = XValidL;
    CurrentObject.XTestUL = XValidUL;
    CurrentObject.YTestL = YValidL;
    CurrentObject.I = ones(size(XValidL{1,1},2),1);
    CurrentObject.MakeHistPlot = 'True';
    CurrentObject.GenerateFigure = 'True';
    ModelResults = CurrentObject.RegressModel
    save(['Results/Object_Results/Step',num2str(Step),'.mat'],'CurrentObject');
    set(gcf, 'PaperPositionMode', 'auto')
    print(figure(18),['Results/Figures_FINAL/Step',num2str(Step),'_Sep_FINAL'],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
    print(figure(22),['Results/Figures_FINAL/Step',num2str(Step),'_Comb_FINAL'],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
else
    disp(['Object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
    load(['Results/Object_Results/Step',num2str(Step),'.mat']);
end

%% STEP 4: GENETIC ALGORITHM + REGRESSION
% At this step, the genetic algorithm (GA) is used to automatically select the most informative features prior to application of linear regression
% Include all features again before genetic algorithm starts:
close all
clear XTrainL YTrainL XValidL YValidL XValidUL
FeatureInds = 1:size(DATA.FeatureNames(~strcmp(DATA.FeatureNames,'Time')),2); % Redefine FeatureInds back to the original verision to include all features again
FeatureInds(TargetSet) = []; %Remove indicies for tragets.
[XTrainL,XTrainUL,YTrainL,XValidL,XValidUL,YValidL] = DataCollector_Concat_TrainValid(NCellsL,NCellsUL,FeatureInds,TargetInd,TrainSet,ValidSet,DATA);

%Apply genetic algorithm to automatically select the most informative features with respect to errors in prediction of unlabeled cells lipid content:
NFeatures = size(XTrainL{1,1},2);  % Record number of features.
if ~exist(['Results/Object_Results/Step',num2str(Step+1),'.mat'],'file')
    clear CurrentObject;
    load('Results/Object_Results/Step1.mat'); % Load the result for step 1 (linear regression with all features)
    ObjFun = @(I)(CurrentObject.Compute_SepKS(CurrentObject.XTrainL,CurrentObject.YTrainL,CurrentObject.XTestUL,CurrentObject.YTestL,I)); % Define the objective funtion to be used in the GA
    PopInit = 1*(rand(50*NFeatures,NFeatures)<0.1); % Define the initial population for GA iterations
    options = optimoptions('ga','PopulationType','bitstring',...
        'UseParallel',true,'Display','final',...
        'InitialPopulationMatrix',PopInit,...
        'PopulationSize',50*NFeatures,...
        'MutationFcn',{@mutationuniform, 0.5/NFeatures}); % Specify options for the GA
    IvecBest = ga(ObjFun,NFeatures,options); %OR: Ivec_Best = ga(Reg_Error,N_Features,[],[],[],[],[],[],[],options); % Run the GA for the model
else
    clear CurrentObject;
    disp(['Object for Step',num2str(Step+1),' already exists. Loading Step',num2str(Step+1),' object.']);
    load(['Results/Object_Results/Step',num2str(Step+1),'.mat']);
    IvecBest = CurrentObject.I; % Assing the result of the GA (selected features) as an input to the matlab object, to be used in the linear regression
end
% After application of the genetic algorithm, do a regression analysis:
clear CurrentObject;
Step = 4;
if ~exist(['Results/Object_Results/Step',num2str(Step),'.mat'],'file')
    CurrentObject = FCSC_Log;
    CurrentObject.Step = Step;
    CurrentObject.Stage = 'TrainValid';
    CurrentObject.BINS = BINS;
    CurrentObject.XTrainL = XTrainL;
    CurrentObject.XTrainUL = XTrainUL;
    CurrentObject.YTrainL = YTrainL;
    CurrentObject.XTestL = XValidL;
    CurrentObject.XTestUL = XValidUL;
    CurrentObject.YTestL = YValidL;
    CurrentObject.I = IvecBest';
    CurrentObject.MakeHistPlot = 'True';
    CurrentObject.GenerateFigure = 'True';
    ModelResults = CurrentObject.RegressModel
    save(['Results/Object_Results/Step',num2str(Step),'.mat'],'CurrentObject');
    set(gcf, 'PaperPositionMode', 'auto')
    print(figure(18),['Results/Figures_FINAL/Step',num2str(Step),'_Sep_FINAL'],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
    print(figure(22),['Results/Figures_FINAL/Step',num2str(Step),'_Comb_FINAL'],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
else
    disp(['Object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
    load(['Results/Object_Results/Step',num2str(Step),'.mat']);
end

%% STEP 5: QUADRATIC FEATURES + REGRESSION
% The quadratic method: at this step, linear regression is applied on an expanded matrix which includes the features, their quadratic terms, and first order products
close all;
[XTrainLQ,XTrainULQ,XValidLQ,XValidULQ] = deal({});
for i=1:size(XTrainL,1)
    XTrainLQ{i,1} = Get_Quad_Features(XTrainL{i,1}); % Get the quadratic values and add them to the original matrix of features using the function 'Get_Quad_Features' (labeled training)
    XTrainULQ{i,1} = Get_Quad_Features(XTrainUL{i,1}); % Unlabeled training
    XValidLQ{i,1} = Get_Quad_Features(XValidL{i,1}); % Labeled validation
    XValidULQ{i,1} = Get_Quad_Features(XValidUL{i,1}); % Unlabeled validation
end
clear CurrentObject;
Step = 5;
if ~exist(['Results/Object_Results/Step',num2str(Step),'.mat'],'file')
    CurrentObject = FCSC_Log;
    CurrentObject.Step = Step;
    CurrentObject.Stage = 'TrainValid';
    CurrentObject.BINS = BINS;
    CurrentObject.XTrainL = XTrainLQ;
    CurrentObject.XTrainUL = XTrainULQ;
    CurrentObject.YTrainL = YTrainL;
    CurrentObject.XTestL = XValidLQ;
    CurrentObject.XTestUL = XValidULQ;
    CurrentObject.YTestL = YValidL;
    CurrentObject.I = ones(size(XValidLQ{1,1},2),1);
    CurrentObject.MakeHistPlot = 'True';
    CurrentObject.GenerateFigure = 'True';
    ModelResults = CurrentObject.RegressModel
    save(['Results/Object_Results/Step',num2str(Step),'.mat'],'CurrentObject');
    set(gcf, 'PaperPositionMode', 'auto')
    print(figure(18),['Results/Figures_FINAL/Step',num2str(Step),'_Sep_FINAL'],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
    print(figure(22),['Results/Figures_FINAL/Step',num2str(Step),'_Comb_FINAL'],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
else
    disp(['Object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
    load(['Results/Object_Results/Step',num2str(Step),'.mat']);
end

%% STEP 6: GENETIC ALGORITHM + QUADRATIC FEATURES + REGRESSION (WITH L1 REGULARIZATION)
% At this step, the GA is applied on the quadratic matrix of features to select the most informative attributes prior to regression analysis. The rest is the same as step 4.
close all;
NFeatures = size(XTrainLQ{1,1},2);  %Record number of features.
if ~exist(['Results/Object_Results/Step',num2str(Step+1),'.mat'],'file')
    clear CurrentObject;
    load('Results/Object_Results/Step5.mat');
    ObjFun = @(I)(CurrentObject.Compute_SepKS(CurrentObject.XTrainL,CurrentObject.YTrainL,CurrentObject.XTestUL,CurrentObject.YTestL,I)); %NO regularization here!
    PopInit = 1*(rand(50*NFeatures,NFeatures)<0.1);
    options = optimoptions('ga','PopulationType','bitstring',...
        'UseParallel',true,'Display','final',...
        'InitialPopulationMatrix',PopInit,...
        'PopulationSize',50*NFeatures,...
        'MutationFcn',{@mutationuniform, 0.5/NFeatures});
    IvecBest = ga(ObjFun,NFeatures,options);
else
    clear CurrentObject;
    disp(['Object for Step',num2str(Step+1),' already exists. Loading Step',num2str(Step+1),' object.']);
    load(['Results/Object_Results/Step',num2str(Step+1),'.mat']);
    IvecBest = CurrentObject.I;
end
% After application of the genetic algorithm, Do a regression analysis as usual:
clear CurrentObject;
Step = 6;
if ~exist(['Results/Object_Results/Step',num2str(Step),'.mat'],'file')
    CurrentObject = FCSC_Log;
    CurrentObject.Step = Step;
    CurrentObject.Stage = 'TrainValid';
    CurrentObject.BINS = BINS;
    CurrentObject.XTrainL = XTrainLQ;
    CurrentObject.XTrainUL = XTrainULQ;
    CurrentObject.YTrainL = YTrainL;
    CurrentObject.XTestL = XValidLQ;
    CurrentObject.XTestUL = XValidULQ;
    CurrentObject.YTestL = YValidL;
    CurrentObject.I = IvecBest;
    CurrentObject.MakeHistPlot = 'True';
    CurrentObject.GenerateFigure = 'True';
    ModelResults = CurrentObject.RegressModel
    save(['Results/Object_Results/Step',num2str(Step),'.mat'],'CurrentObject');
    set(gcf, 'PaperPositionMode', 'auto')
    print(figure(18),['Results/Figures_FINAL/Step',num2str(Step),'_Sep_FINAL'],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
    print(figure(22),['Results/Figures_FINAL/Step',num2str(Step),'_Comb_FINAL'],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
else
    disp(['Object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
    load(['Results/Object_Results/Step',num2str(Step),'.mat']);
end

%% SECTION IV: TESTING THE MODELS GATHERED FROM CONCATENATED TRAINING & VALIDAION DATA
% In this section, we test our fully developed model on all the remaining data including replications of measurements. No training or validation is performed here.
for TestLRep = 1:size(DATA.Stained,2)
    for TestULRep = 1:size(DATA.Unstained,2)
        close all;
        [XTestL,YTestL,XTestUL] = DataCollector_Test(NCells,FeatureInds,TargetInd,TestSet,TestLRep,TestULRep,DATA); %Extract data

        %% TESTING STEP 1
        clear CurrentObject;
        Step = 1;
        if ~exist(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'file')
            CurrentObject = FCSC_Log;
            CurrentObject.Step = Step;
            CurrentObject.Stage = 'Test';
            CurrentObject.BINS = BINS;
            CurrentObject.XTrainL = XTrainL;
            CurrentObject.XTrainUL = XTrainUL;
            CurrentObject.YTrainL = YTrainL;
            CurrentObject.XTestL = XTestL;
            CurrentObject.XTestUL = XTestUL;
            CurrentObject.YTestL = YTestL;
            CurrentObject.I = ones(size(XTestL{1,1},2),1);
            CurrentObject.MakeHistPlot = 'True';
            CurrentObject.GenerateFigure = 'True';
            ModelResults = CurrentObject.RegressModel
            save(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'CurrentObject');
            set(gcf, 'PaperPositionMode', 'auto')
            print(figure(74),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Labeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(75),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Unlabeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(79),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Comb_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
        else
            disp(['For the current replication, object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
            load(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat']);
        end
        
        %% TESTING STEP 3
        CorruptedFeatures = [4,10]; %User input
        for i=1:length(CorruptedFeatures)
            FeatureInds(FeatureInds==CorruptedFeatures(i)) = [];
        end
        %Redo data extraction after introducing new feature set (without corrupted features)
        clear XTrainL YTrainL XTestL YTestL XTestUL
        [XTrainL,XTrainUL,YTrainL,~,~,~] = DataCollector_Concat_TrainValid(NCellsL,NCellsUL,FeatureInds,TargetInd,TrainSet,ValidSet,DATA);
        [XTestL,YTestL,XTestUL] = DataCollector_Test(NCells,FeatureInds,TargetInd,TestSet,TestLRep,TestULRep,DATA); %Extract datasets.
        close all; clear CurrentObject;
        Step = 3;
        if ~exist(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'file')
            CurrentObject = FCSC_Log;
            CurrentObject.Step = Step;
            CurrentObject.Stage = 'Test';
            CurrentObject.BINS = BINS;
            CurrentObject.XTrainL = XTrainL;
            CurrentObject.XTrainUL = XTrainUL;
            CurrentObject.YTrainL = YTrainL;
            CurrentObject.XTestL = XTestL;
            CurrentObject.XTestUL = XTestUL;
            CurrentObject.YTestL = YTestL;
            CurrentObject.I = ones(size(XTestL{1,1},2),1);
            CurrentObject.MakeHistPlot = 'True';
            CurrentObject.GenerateFigure = 'True';
            ModelResults = CurrentObject.RegressModel
            save(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'CurrentObject');
            set(gcf, 'PaperPositionMode', 'auto')
            print(figure(74),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Labeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(75),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Unlabeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(79),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Comb_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
        else
            disp(['For the current replication, object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
            load(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat']);
        end
        
        %% TESTING STEP 4
        %Include all features again before genetic algorithm starts:
        close all;
        clear XTrainL YTrainL XTestL YTestL XTestUL
        FeatureInds = 1:size(DATA.FeatureNames(~strcmp(DATA.FeatureNames,'Time')),2); %Reform FeatureInds back to the original verision to include all features again.
        FeatureInds(TargetSet) = []; %Remove indicies for tragets.
        [XTrainL,XTrainUL,YTrainL,~,~,~] = DataCollector_Concat_TrainValid(NCellsL,NCellsUL,FeatureInds,TargetInd,TrainSet,ValidSet,DATA);
        [XTestL,YTestL,XTestUL] = DataCollector_Test(NCells,FeatureInds,TargetInd,TestSet,TestLRep,TestULRep,DATA); %Extract datasets.

        %Call the model defined by the genetic algorithm in step 4 and extract the information regarding which features selected:
        Step = 4;
        clear CurrentObject;
        load(['Results/Object_Results/Step',num2str(Step),'.mat']);
        IvecBest = CurrentObject.I;
        %Do the testing for the model defined in step 4 by genetic algorithm:
        clear CurrentObject;
        if ~exist(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'file')
            CurrentObject = FCSC_Log;
            CurrentObject.Step = Step;
            CurrentObject.Stage = 'Test';
            CurrentObject.BINS = BINS;
            CurrentObject.XTrainL = XTrainL;
            CurrentObject.XTrainUL = XTrainUL;
            CurrentObject.YTrainL = YTrainL;
            CurrentObject.XTestL = XTestL;
            CurrentObject.XTestUL = XTestUL;
            CurrentObject.YTestL = YTestL;
            CurrentObject.I = IvecBest;
            CurrentObject.MakeHistPlot = 'True';
            CurrentObject.GenerateFigure = 'True';
            ModelResults = CurrentObject.RegressModel
            save(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'CurrentObject');
            set(gcf, 'PaperPositionMode', 'auto')
            print(figure(74),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Labeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(75),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Unlabeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(79),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Comb_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
        else
            disp(['For the current replication, object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
            load(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat']);
        end
        
        %% TESTING STEP 5
        close all;
        [XTrainLQ,XTrainULQ,XTestLQ,XTestULQ] = deal({});
        for i=1:size(XTrainL,1)
            XTrainLQ{i,1} = Get_Quad_Features(XTrainL{i,1}); %Getting quadratic XTrainL with the function Get_Quad_Features.
            XTrainULQ{i,1} = Get_Quad_Features(XTrainUL{i,1}); %Getting quadratic XTrainUL with the function Get_Quad_Features.
        end
        for i=1:size(XTestL,1)
            XTestLQ{i,1} = Get_Quad_Features(XTestL{i,1}); %Getting quadratic XtestL with the function Get_Quad_Features.
            XTestULQ{i,1} = Get_Quad_Features(XTestUL{i,1}); %Getting quadratic XtestUL with the function Get_Quad_Features.
        end
        clear CurrentObject;
        Step = 5;
        if ~exist(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'file')
            CurrentObject = FCSC_Log;
            CurrentObject.Step = Step;
            CurrentObject.Stage = 'Test';
            CurrentObject.BINS = BINS;
            CurrentObject.XTrainL = XTrainLQ;
            CurrentObject.XTrainUL = XTrainULQ;
            CurrentObject.YTrainL = YTrainL;
            CurrentObject.XTestL = XTestLQ;
            CurrentObject.XTestUL = XTestULQ;
            CurrentObject.YTestL = YTestL;
            CurrentObject.I = ones(size(XTestLQ{1,1},2),1);
            CurrentObject.MakeHistPlot = 'True';
            CurrentObject.GenerateFigure = 'True';
            ModelResults = CurrentObject.RegressModel
            save(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'CurrentObject');
            set(gcf, 'PaperPositionMode', 'auto')
            print(figure(74),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Labeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(75),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Unlabeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(79),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Comb_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
        else
            disp(['For the current replication, object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
            load(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat']);
        end
        
        %% TESTING STEP 6:
        %Call the model defined by the genetic algorithm in step 6 and extract the information regarding which features selected:
        close all;
        Step = 6;
        clear CurrentObject;
        load(['Results/Object_Results/Step',num2str(Step),'.mat']);
        IvecBest = CurrentObject.I;
        %Do the testing for the model defined in step 6 by genetic algorithm:
        clear CurrentObject;
        if ~exist(['Results/Object_Results/Step',num2str(Step),'_ValidLRep',num2str(TestLRep),'_ValidULRep',num2str(TestULRep),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'file')
            CurrentObject = FCSC_Log;
            CurrentObject.Step = Step;
            CurrentObject.Stage = 'Test';
            CurrentObject.BINS = BINS;
            CurrentObject.XTrainL = XTrainLQ;
            CurrentObject.XTrainUL = XTrainULQ;
            CurrentObject.YTrainL = YTrainL;
            CurrentObject.XTestL = XTestLQ;
            CurrentObject.XTestUL = XTestULQ;
            CurrentObject.YTestL = YTestL;
            CurrentObject.I = IvecBest;
            CurrentObject.MakeHistPlot = 'True';
            CurrentObject.GenerateFigure = 'True';
            ModelResults = CurrentObject.RegressModel
            save(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'CurrentObject');
            set(gcf, 'PaperPositionMode', 'auto')
            print(figure(74),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Labeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(75),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Unlabeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(79),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Comb_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
        else
            disp(['For the current replication, object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
            load(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat']);
        end
    end
end
