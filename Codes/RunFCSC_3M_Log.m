%% RunFCSC_3M_Log
%Author: Mohammad Tanhaemami
%Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.

%This script defines a weighted model as M = a1M1 + a2M2 + a3M3 (a1, a2, and a3 are the weights)
%This script does not plot average lipid content by day (for that, find "AvgLipidByDay.m" in the current directory).
%This script does not plot bodipy signals, hence it only accepts manual input for training, validation, and testing data sets (plot bodipy signals, find "AvgLipidByDay.m").
%This script assigns 3 time points for training and validation.

%% INITIALIZATION, LOAD THE DESIRED DATABASE, AND ADD REQUIRED PATHS
clear; close all; clc;
warning('off');
load('../Data_Files/NewPico_Data_NoZero.mat'); %Load the file NewPico_Data_NoZero.mat.
path(path,'./Functions'); %Path to funcitions.
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
BINS = linspace(log(DATA.Max(TargetInd))-10,log(DATA.Max(TargetInd))+40,500);
% BINS = linspace(log(DATA.Max(TargetInd))-40,log(DATA.Max(TargetInd))+40,500);
TrainSet = [2,10,23]; %[SortInd(1) SortInd(end)]; %Widest range in Bodipy_Signals.
ValidSet = [1,11,22]; %[SortInd(2) SortInd(end-1)]; %Second widest range in Bodipy_Signals.
TestSet = [3:9,12:21]; %sort(SortInd(3:end-2)); %Define testing dataset (not used before).
FeatureInds = 1:size(DATA.FeatureNames(~strcmp(DATA.FeatureNames,'Time')),2); %Define Feature_Inds
FeatureInds(TargetSet) = []; %Remove indicies for tragets.
FeatureInds_S = sort([FeatureInds,TargetSet]); %Define Feature_Inds for S test statistic

%% SECTION II: COLLECT DATA FOR ALL COMBINATIONS OF REPLICATIONS IN TRAINING & VALIDATION DATA SETS
for mmm = 1:2
    % Collect the data for further analysis
    clear XTrainL XTrainUL YTrainL XValidL XValidUL YValidL
    [XTrainL,XTrainUL,XTrainUL_S,YTrainL,XValidL,XValidUL,XValidUL_S,YValidL] = DataCollector_Concat_TrainValid(NCellsL,NCellsUL,FeatureInds,FeatureInds_S,TargetInd,TrainSet,ValidSet,DATA);
    
    %% STEP 1: REGRESSION - PREDICT LIPID CONTENT WITH ALL FEATURES
    % At this step, linear regression is used to predict the lipid content using all available features
    close all; clear Alphas CurrentObject;
    Step = 1;
    if ~exist(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat'],'file')
        FeaturesType = 'Linear'; % 'Linear' or 'Quadratic' options for the matrix of features
        I = ones(size(XValidL{1,1},2),size(XValidL,1)); % Initial value for the vector that indicates which features are to be used by the GA
        if mmm == 1
            PrevBestColS_Step1 = [];
        else
            PrevBestColS_Step1 = PrevBestColS_Step1(mmm-1,:);
        end
        Alphas = Alpha_Finder(DATA,I,FeatureInds,FeatureInds_S,FeaturesType,PrevBestColS_Step1); % Compute the regression coefficient, optimized weight quotient, and the best columns of the test statistics matrix (S) (or rows of the weight quotient (Q))
        PrevBestColS_Step1(mmm,:) = Alphas.IvecBest_is; % Set the initial value of the best columns of S (rows of Q)
        CurrentObject = FCSC_3M_Log; % Initialize the matlab object
        CurrentObject.Step = Step; % Set the step of the calculations (here stpe 1 for simple linear regression)
        CurrentObject.Stage = 'TrainValid'; % Set the stage for calculations (here is the training and validation stage)
        CurrentObject.BINS = BINS; % Assign bins
        CurrentObject.XTrainL = XTrainL; % Assign labeled training features
        CurrentObject.XTrainUL = XTrainUL; % Assing unlabeled training features
        CurrentObject.XTrainUL_S = XTrainUL_S; % Assign unlabeled training features including the main channel for the secondary regression
        CurrentObject.YTrainL = YTrainL; % Assign the measured labeled trainsing taragets
        CurrentObject.XTestL = XValidL; % Assign labeled validation features
        CurrentObject.XTestUL = XValidUL; % Assing unlabled validation features
        CurrentObject.XTestUL_S = XValidUL_S; % Assing unlabeled validation feature including the main channel
        CurrentObject.YTestL = YValidL; % Assign the measured labeled validation targets
        CurrentObject.I = I; % Assing the vector of selected features
        CurrentObject.M = Alphas.M; % Assign regression coeffcients (the 3 models from training from the function 'Alpha_Finder')
        CurrentObject.Q = Alphas.Q; % Assign the weight qoutients from the function 'Alpha_Finder'
        CurrentObject.IforS = Alphas.IvecBest_is; % Assign the vector of selected columns of S (rows of Q) from the function 'Alpha_Finder'
        CurrentObject.MakeHistPlot = 'True'; % Logical set to make the histogram plots
        CurrentObject.GenerateFigure = 'True'; % Logical set to generate figures for results
        ModelResults = CurrentObject.RegressModel
        save(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat'],'CurrentObject');
        set(gcf, 'PaperPositionMode', 'auto');
        print(figure(25),['Results/Figures_FINAL/Step',num2str(Step),'_Sep_FINAL','_TRY_',num2str(mmm)],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
    else
        disp(['Object for Step',num2str(Step),'_TRY_',num2str(mmm),' already exists. Loading Step',num2str(Step),'_TRY_',num2str(mmm),' object.']);
        load(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat']);
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
    close all; clear Alphas CurrentObject;
    CorruptedFeatures = [4,10]; % User input for corrupted features to be removed from analysis
    for i = 1:length(CorruptedFeatures)
        FeatureInds(FeatureInds == CorruptedFeatures(i)) = [];
    end
    FeatureInds_S = sort([FeatureInds,TargetSet]); %Redefine Feature_Inds for S test statistic after corrupted features removed
    %Redo data extraction after introducing new feature set (without corrupted features)
    clear XTrainL XTrainUL YTrainL XValidL XValidUL YValidL
    [XTrainL,XTrainUL,XTrainUL_S,YTrainL,XValidL,XValidUL,XValidUL_S,YValidL] = DataCollector_Concat_TrainValid(NCellsL,NCellsUL,FeatureInds,FeatureInds_S,TargetInd,TrainSet,ValidSet,DATA);
    Step = 3;
    if ~exist(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat'],'file')
        I = ones(size(XValidL{1,1},2),size(XValidL,1));
        if mmm == 1
            PrevBestColS_Step3 = [];
        else
            PrevBestColS_Step3 = PrevBestColS_Step3(mmm-1,:);
        end
        Alphas = Alpha_Finder(DATA,I,FeatureInds,FeatureInds_S,FeaturesType,PrevBestColS_Step3);
        PrevBestColS_Step3(mmm,:) = Alphas.IvecBest_is;
        CurrentObject = FCSC_3M_Log;
        CurrentObject.Step = Step;
        CurrentObject.Stage = 'TrainValid';
        CurrentObject.BINS = BINS;
        CurrentObject.XTrainL = XTrainL;
        CurrentObject.XTrainUL = XTrainUL;
        CurrentObject.XTrainUL_S = XTrainUL_S;
        CurrentObject.YTrainL = YTrainL;
        CurrentObject.XTestL = XValidL;
        CurrentObject.XTestUL = XValidUL;
        CurrentObject.XTestUL_S = XValidUL_S;
        CurrentObject.YTestL = YValidL;
        CurrentObject.I = I;
        CurrentObject.M = Alphas.M;
        CurrentObject.Q = Alphas.Q;
        CurrentObject.IforS = Alphas.IvecBest_is;
        CurrentObject.MakeHistPlot = 'True';
        CurrentObject.GenerateFigure = 'True';
        ModelResults = CurrentObject.RegressModel
        save(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat'],'CurrentObject');
        set(gcf, 'PaperPositionMode', 'auto');
        print(figure(25),['Results/Figures_FINAL/Step',num2str(Step),'_Sep_FINAL','_TRY_',num2str(mmm)],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
    else
        disp(['Object for Step',num2str(Step),'_TRY_',num2str(mmm),' already exists. Loading Step',num2str(Step),'_TRY_',num2str(mmm),' object.']);
        load(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat']);
    end
    
    %% STEP 4: GENETIC ALGORITHM + REGRESSION
    % At this step, the genetic algorithm (GA) is used to automatically select the most informative features prior to application of linear regression
    % Include all features again before genetic algorithm starts:
    close all;
    clear XTrainL XTrainUL YTrainL XValidL XValidUL YValidL
    FeatureInds = 1:size(DATA.FeatureNames(~strcmp(DATA.FeatureNames,'Time')),2); % Redefine FeatureInds back to the original verision to include all features again
    FeatureInds(TargetSet) = []; % Remove indices for tragets.
    FeatureInds_S = sort([FeatureInds,TargetSet]); % Define Feature_Inds for S test statistic
    [XTrainL,XTrainUL,XTrainUL_S,YTrainL,XValidL,XValidUL,XValidUL_S,YValidL] = DataCollector_Concat_TrainValid(NCellsL,NCellsUL,FeatureInds,FeatureInds_S,TargetInd,TrainSet,ValidSet,DATA); % Collect the data features (all features) and targets
    
    % Apply genetic algorithm to automatically select the most informative features with respect to errors in prediction of unlabeled cells lipid content:
    Step = 4;
    NFeatures = size(XTrainL{1,1},2); % Record number of features.
    if ~exist(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat'],'file')
        clear CurrentObject;
        load(['Results/Object_Results/Step1_TRY_',num2str(mmm),'.mat']); % Load the result for step 1 (linear regression with all features)
        ObjFun1 = @(I)(CurrentObject.Compute_KS(CurrentObject.XTrainL{1,1},CurrentObject.YTrainL{1,1},CurrentObject.XTestUL{1,1},CurrentObject.YTestL{1,1},I)); % Define the objective funtion for model 1, to be used in the GA
        ObjFun2 = @(I)(CurrentObject.Compute_KS(CurrentObject.XTrainL{2,1},CurrentObject.YTrainL{2,1},CurrentObject.XTestUL{2,1},CurrentObject.YTestL{2,1},I)); % Define the objective funtion for model 2, to be used in the GA
        ObjFun3 = @(I)(CurrentObject.Compute_KS(CurrentObject.XTrainL{3,1},CurrentObject.YTrainL{3,1},CurrentObject.XTestUL{3,1},CurrentObject.YTestL{3,1},I)); % Define the objective funtion for model 3, to be used in the GA
        PopInit = 1*(rand(50*NFeatures,NFeatures)<0.1); % Define the initial population for GA iterations
        options = optimoptions('ga','PopulationType','bitstring',...
            'UseParallel',true,'Display','final',...
            'InitialPopulationMatrix',PopInit,...
            'PopulationSize',50*NFeatures,...
            'MutationFcn',{@mutationuniform, 0.5/NFeatures}); % Specify options for the GA
        IvecBest1 = ga(ObjFun1,NFeatures,options); % OR: Ivec_Best = ga(Reg_Error,N_Features,[],[],[],[],[],[],[],options); % Run the GA for the 1st  model
        IvecBest2 = ga(ObjFun2,NFeatures,options); % OR: Ivec_Best = ga(Reg_Error,N_Features,[],[],[],[],[],[],[],options); % Run the GA for the 2nd  model
        IvecBest3 = ga(ObjFun3,NFeatures,options); % OR: Ivec_Best = ga(Reg_Error,N_Features,[],[],[],[],[],[],[],options); % Run the GA for the 3rd  model
    else
        clear CurrentObject;
        disp(['Object for Step',num2str(Step),'_TRY_',num2str(mmm),' already exists. Loading Step',num2str(Step),'_TRY_',num2str(mmm),' object.']);
        load(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat']);
        IvecBest1 = CurrentObject.I(:,1); % Assing the result of the GA (selected features) as an input to the matlab object, to be used in the linear regression (Model 1)
        IvecBest2 = CurrentObject.I(:,2); % Assing the result of the GA (selected features) as an input to the matlab object, to be used in the linear regression (Model 2)
        IvecBest3 = CurrentObject.I(:,3); % Assing the result of the GA (selected features) as an input to the matlab object, to be used in the linear regression (Model 3)
    end
    % After application of the genetic algorithm, do a regression analysis:
    clear I Alphas CurrentObject;
    if ~exist(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat'],'file')
        I(:,1) = IvecBest1';
        I(:,2) = IvecBest2';
        I(:,3) = IvecBest3';
        if mmm == 1
            PrevBestColS_Step4 = [];
        else
            PrevBestColS_Step4 = PrevBestColS_Step4(mmm-1,:);
        end
        Alphas = Alpha_Finder(DATA,I,FeatureInds,FeatureInds_S,FeaturesType,PrevBestColS_Step4);
        PrevBestColS_Step4(mmm,:) = Alphas.IvecBest_is;
        CurrentObject = FCSC_3M_Log;
        CurrentObject.Step = Step;
        CurrentObject.Stage = 'TrainValid';
        CurrentObject.BINS = BINS;
        CurrentObject.XTrainL = XTrainL;
        CurrentObject.XTrainUL = XTrainUL;
        CurrentObject.XTrainUL_S = XTrainUL_S;
        CurrentObject.YTrainL = YTrainL;
        CurrentObject.XTestL = XValidL;
        CurrentObject.XTestUL = XValidUL;
        CurrentObject.XTestUL_S = XValidUL_S;
        CurrentObject.YTestL = YValidL;
        CurrentObject.I(:,1) = IvecBest1';
        CurrentObject.I(:,2) = IvecBest2';
        CurrentObject.I(:,3) = IvecBest3';
        CurrentObject.M = Alphas.M;
        CurrentObject.Q = Alphas.Q;
        CurrentObject.IforS = Alphas.IvecBest_is;
        CurrentObject.MakeHistPlot = 'True';
        CurrentObject.GenerateFigure = 'True';
        ModelResults = CurrentObject.RegressModel
        save(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat'],'CurrentObject');
        set(gcf, 'PaperPositionMode', 'auto');
        print(figure(25),['Results/Figures_FINAL/Step',num2str(Step),'_Sep_FINAL','_TRY_',num2str(mmm)],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
    else
        disp(['Object for Step',num2str(Step),'_TRY_',num2str(mmm),' already exists. Loading Step',num2str(Step),'_TRY_',num2str(mmm),' object.']);
        load(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat']);
    end
    
    %% STEP 5: QUADRATIC FEATURES + REGRESSION
    % The quadratic method: at this step, linear regression is applied on an expanded matrix which includes the features, their quadratic terms, and first order products
    close all;
    [XTrainLQ,XTrainULQ,XTrainULQ_S,XValidLQ,XValidULQ,XValidULQ_S] = deal({});
    for i=1:size(XTrainL,1)
        XTrainLQ{i,1} = Get_Quad_Features(XTrainL{i,1}); % Get the quadratic values and add them to the original matrix of features using the function 'Get_Quad_Features' (labeled training)
        XTrainULQ{i,1} = Get_Quad_Features(XTrainUL{i,1}); % Unlabeled training
        XTrainULQ_S{i,1} = Get_Quad_Features(XTrainUL_S{i,1}); % Unlabeled training including main channel information for secondary regression analysis
        XValidLQ{i,1} = Get_Quad_Features(XValidL{i,1}); % Labeled validation
        XValidULQ{i,1} = Get_Quad_Features(XValidUL{i,1}); % Unlabeled validation
        XValidULQ_S{i,1} = Get_Quad_Features(XValidUL_S{i,1}); % Unlabeled validation including main channel information for secondary regression analysis
    end
    clear I Alphas CurrentObject;
    Step = 5;
    if ~exist(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat'],'file')
        FeaturesType = 'Quadratic'; % 'Linear' or 'Quadratic' features type
        I = ones(size(XValidLQ{1,1},2),size(XValidLQ,1));
        if mmm == 1
            PrevBestColS_Step5 = [];
        else
            PrevBestColS_Step5 = PrevBestColS_Step5(mmm-1,:);
        end
        Alphas = Alpha_Finder(DATA,I,FeatureInds,FeatureInds_S,FeaturesType,PrevBestColS_Step5);
        PrevBestColS_Step5(mmm,:) = Alphas.IvecBest_is;
        CurrentObject = FCSC_3M_Log;
        CurrentObject.Step = Step;
        CurrentObject.Stage = 'TrainValid';
        CurrentObject.BINS = BINS;
        CurrentObject.XTrainL = XTrainLQ;
        CurrentObject.XTrainUL = XTrainULQ;
        CurrentObject.XTrainUL_S = XTrainULQ_S;
        CurrentObject.YTrainL = YTrainL;
        CurrentObject.XTestL = XValidLQ;
        CurrentObject.XTestUL = XValidULQ;
        CurrentObject.XTestUL_S = XValidULQ_S;
        CurrentObject.YTestL = YValidL;
        CurrentObject.I = I;
        CurrentObject.M = Alphas.M;
        CurrentObject.Q = Alphas.Q;
        CurrentObject.IforS = Alphas.IvecBest_is;
        CurrentObject.MakeHistPlot = 'True';
        CurrentObject.GenerateFigure = 'True';
        ModelResults = CurrentObject.RegressModel
        save(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat'],'CurrentObject');
        set(gcf, 'PaperPositionMode', 'auto');
        print(figure(25),['Results/Figures_FINAL/Step',num2str(Step),'_Sep_FINAL','_TRY_',num2str(mmm)],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
    else
        disp(['Object for Step',num2str(Step),'_TRY_',num2str(mmm),' already exists. Loading Step',num2str(Step),'_TRY_',num2str(mmm),' object.']);
        load(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat']);
    end
    
    %% STEP 6: GENETIC ALGORITHM + QUADRATIC FEATURES + REGRESSION (WITH L1 REGULARIZATION)
    % At this step, the GA is applied on the quadratic matrix of features to select the most informative attributes prior to regression analysis. The rest is the same as step 4.
    close all;
    Step = 6;
    NFeatures = size(XTrainLQ{1,1},2);  %Record number of features.
    if ~exist(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat'],'file')
        clear CurrentObject;
        load(['Results/Object_Results/Step5_TRY_',num2str(mmm),'.mat']);
        ObjFun1 = @(I)(CurrentObject.Compute_KS(CurrentObject.XTrainL{1,1},CurrentObject.YTrainL{1,1},CurrentObject.XTestUL{1,1},CurrentObject.YTestL{1,1},I));
        ObjFun2 = @(I)(CurrentObject.Compute_KS(CurrentObject.XTrainL{2,1},CurrentObject.YTrainL{2,1},CurrentObject.XTestUL{2,1},CurrentObject.YTestL{2,1},I));
        ObjFun3 = @(I)(CurrentObject.Compute_KS(CurrentObject.XTrainL{3,1},CurrentObject.YTrainL{3,1},CurrentObject.XTestUL{3,1},CurrentObject.YTestL{3,1},I));
        PopInit = 1*(rand(50*NFeatures,NFeatures)<0.1);
        options = optimoptions('ga','PopulationType','bitstring',...
            'UseParallel',true,'Display','final',...
            'InitialPopulationMatrix',PopInit,...
            'PopulationSize',50*NFeatures,...
            'MutationFcn',{@mutationuniform, 0.5/NFeatures});
        IvecBest1 = ga(ObjFun1,NFeatures,options);
        IvecBest2 = ga(ObjFun2,NFeatures,options);
        IvecBest3 = ga(ObjFun3,NFeatures,options);
    else
        clear CurrentObject;
        disp(['Object for Step',num2str(Step),'_TRY_',num2str(mmm),' already exists. Loading Step',num2str(Step),'_TRY_',num2str(mmm),' object.']);
        load(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat']);
        IvecBest = CurrentObject.I;
    end
    % After application of the genetic algorithm, Do a regression analysis as usual:
    clear I Alphas CurrentObject;
    if ~exist(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat'],'file')
        I(:,1) = IvecBest1';
        I(:,2) = IvecBest2';
        I(:,3) = IvecBest3';
        if mmm == 1
            PrevBestColS_Step6 = [];
        else
            PrevBestColS_Step6 = PrevBestColS_Step6(mmm-1,:);
        end
        Alphas = Alpha_Finder(DATA,I,FeatureInds,FeatureInds_S,FeaturesType,PrevBestColS_Step6);
        PrevBestColS_Step6(mmm,:) = Alphas.IvecBest_is;
        CurrentObject = FCSC_3M_Log;
        CurrentObject.Step = Step;
        CurrentObject.Stage = 'TrainValid';
        CurrentObject.BINS = BINS;
        CurrentObject.XTrainL = XTrainLQ;
        CurrentObject.XTrainUL = XTrainULQ;
        CurrentObject.XTrainUL_S = XTrainULQ_S;
        CurrentObject.YTrainL = YTrainL;
        CurrentObject.XTestL = XValidLQ;
        CurrentObject.XTestUL = XValidULQ;
        CurrentObject.XTestUL_S = XValidULQ_S;
        CurrentObject.YTestL = YValidL;
        CurrentObject.I(:,1) = IvecBest1';
        CurrentObject.I(:,2) = IvecBest2';
        CurrentObject.I(:,3) = IvecBest3';
        CurrentObject.M = Alphas.M;
        CurrentObject.Q = Alphas.Q;
        CurrentObject.IforS = Alphas.IvecBest_is;
        CurrentObject.MakeHistPlot = 'True';
        CurrentObject.GenerateFigure = 'True';
        ModelResults = CurrentObject.RegressModel
        save(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat'],'CurrentObject');
        set(gcf, 'PaperPositionMode', 'auto');
        print(figure(25),['Results/Figures_FINAL/Step',num2str(Step),'_Sep_FINAL','_TRY_',num2str(mmm)],'-dpdf','-bestfit'); %Saving final figure as a *.pdf
    else
        disp(['Object for Step',num2str(Step),'_TRY_',num2str(mmm),' already exists. Loading Step',num2str(Step),'_TRY_',num2str(mmm),' object.']);
        load(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat']);
    end
end
%% SECTION IV: TESTING THE MODELS GATHERED FROM CONCATENATED TRAINING & VALIDAION DATA
% In this section, we test our fully developed model on all the remaining data including replications of measurements. No training or validation is performed here.
clear Alphas
for TestLRep = 1:size(DATA.Stained,2)
    for TestULRep = 1:size(DATA.Unstained,2)
        close all;
        [XTestL,YTestL,XTestUL,XTestUL_S] = DataCollector_Test(NCells,FeatureInds,FeatureInds_S,TargetInd,TestSet,TestLRep,TestULRep,DATA); %Extract datasets.
        
        %% TESTING STEP 1
        clear CurrentObject;
        Step = 1;
        load(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat']);
        M = CurrentObject.M;
        Q = CurrentObject.Q;
        IforS = CurrentObject.IforS;
        clear CurrentObject;
        if ~exist(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'file')
            CurrentObject = FCSC_3M_Log;
            CurrentObject.Step = Step;
            CurrentObject.Stage = 'Test';
            CurrentObject.BINS = BINS;
            CurrentObject.XTrainL = XTrainL;
            CurrentObject.XTrainUL = XTrainUL;
            CurrentObject.XTrainUL_S = XTrainUL_S;
            CurrentObject.YTrainL = YTrainL;
            CurrentObject.XTestL = XTestL;
            CurrentObject.YTestL = YTestL;
            CurrentObject.XTestUL = XTestUL;
            CurrentObject.XTestUL_S = XTestUL_S;
            CurrentObject.I = ones(size(XTestL{1,1},2),size(XTestL,1));
            CurrentObject.M = M;
            CurrentObject.Q = Q;
            CurrentObject.IforS = IforS;
            CurrentObject.MakeHistPlot = 'True';
            CurrentObject.GenerateFigure = 'True';
            ModelResults = CurrentObject.RegressModel
            save(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'CurrentObject');
            set(gcf, 'PaperPositionMode', 'auto');
            print(figure(81),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Labeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(82),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Unlabeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
        else
            disp(['For the current replication, object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
            load(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat']);
        end
        
        %% TESTING STEP 3
        CorruptedFeatures = [4,10]; % User input
        for i=1:length(CorruptedFeatures)
            FeatureInds(FeatureInds==CorruptedFeatures(i)) = [];
        end
        FeatureInds_S = sort([FeatureInds,TargetSet]); %Redefine Feature_Inds for S test statistic after corrupted features removed
        % Redo data extraction after introducing new feature set (without corrupted features)
        clear XTrainL XTrainUL YTrainL XTestL YTestL XTestUL
        [XTrainL,XTrainUL,XTrainUL_S,YTrainL,~,~,~] = DataCollector_Concat_TrainValid(NCellsL,NCellsUL,FeatureInds,FeatureInds_S,TargetInd,TrainSet,ValidSet,DATA);
        [XTestL,YTestL,XTestUL,XTestUL_S] = DataCollector_Test(NCells,FeatureInds,FeatureInds_S,TargetInd,TestSet,TestLRep,TestULRep,DATA); %Extract datasets.
        close all; clear CurrentObject;
        Step = 3;
        load(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat']);
        M = CurrentObject.M;
        Q = CurrentObject.Q;
        IforS = CurrentObject.IforS;
        clear CurrentObject;
        if ~exist(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'file')
            CurrentObject = FCSC_3M_Log;
            CurrentObject.Step = Step;
            CurrentObject.Stage = 'Test';
            CurrentObject.BINS = BINS;
            CurrentObject.XTrainL = XTrainL;
            CurrentObject.XTrainUL = XTrainUL;
            CurrentObject.XTrainUL_S = XTrainUL_S;
            CurrentObject.YTrainL = YTrainL;
            CurrentObject.XTestL = XTestL;
            CurrentObject.YTestL = YTestL;
            CurrentObject.XTestUL = XTestUL;
            CurrentObject.XTestUL_S = XTestUL_S;
            CurrentObject.I = ones(size(XTestL{1,1},2),size(XTestL,1));
            CurrentObject.M = M;
            CurrentObject.Q = Q;
            CurrentObject.IforS = IforS;
            CurrentObject.MakeHistPlot = 'True';
            CurrentObject.GenerateFigure = 'True';
            ModelResults = CurrentObject.RegressModel
            save(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'CurrentObject');
            set(gcf, 'PaperPositionMode', 'auto');
            print(figure(81),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Labeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(82),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Unlabeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
        else
            disp(['For the current replication, object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
            load(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat']);
        end
        
        %% TESTING STEP 4
        %Include all features again before genetic algorithm starts:
        close all;
        clear XTrainL XTrainUL YTrainL XTestL YTestL XTestUL
        FeatureInds = 1:size(DATA.FeatureNames(~strcmp(DATA.FeatureNames,'Time')),2); %Reform FeatureInds back to the original verision to include all features again.
        FeatureInds(TargetSet) = []; %Remove indicies for tragets.
        FeatureInds_S = sort([FeatureInds,TargetSet]); %Define Feature_Inds for S test statistic
        [XTrainL,XTrainUL,XTrainUL_S,YTrainL,~,~,~] = DataCollector_Concat_TrainValid(NCellsL,NCellsUL,FeatureInds,FeatureInds_S,TargetInd,TrainSet,ValidSet,DATA);
        [XTestL,YTestL,XTestUL,XTestUL_S] = DataCollector_Test(NCells,FeatureInds,FeatureInds_S,TargetInd,TestSet,TestLRep,TestULRep,DATA); %Extract datasets.
        
        %Call the model defined by the genetic algorithm in step 4 and extract necessary information regarding selected features, M, and Q:
        clear CurrentObject;
        Step = 4;
        load(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat']);
        M = CurrentObject.M;
        Q = CurrentObject.Q;
        IforS = CurrentObject.IforS;
        IvecBest1 = CurrentObject.I(:,1);
        IvecBest2 = CurrentObject.I(:,2);
        IvecBest3 = CurrentObject.I(:,3);
        %Do the testing for the model defined in step 4 by genetic algorithm:
        clear CurrentObject;
        if ~exist(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'file')
            CurrentObject = FCSC_3M_Log;
            CurrentObject.Step = Step;
            CurrentObject.Stage = 'Test';
            CurrentObject.BINS = BINS;
            CurrentObject.XTrainL = XTrainL;
            CurrentObject.XTrainUL = XTrainUL;
            CurrentObject.XTrainUL_S = XTrainUL_S;
            CurrentObject.YTrainL = YTrainL;
            CurrentObject.XTestL = XTestL;
            CurrentObject.YTestL = YTestL;
            CurrentObject.XTestUL = XTestUL;
            CurrentObject.XTestUL_S = XTestUL_S;
            CurrentObject.I(:,1) = IvecBest1';
            CurrentObject.I(:,2) = IvecBest2';
            CurrentObject.I(:,3) = IvecBest3';
            CurrentObject.M = M;
            CurrentObject.Q = Q;
            CurrentObject.IforS = IforS;
            CurrentObject.MakeHistPlot = 'True';
            CurrentObject.GenerateFigure = 'True';
            ModelResults = CurrentObject.RegressModel
            save(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'CurrentObject');
            set(gcf, 'PaperPositionMode', 'auto');
            print(figure(81),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Labeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(82),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Unlabeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
        else
            disp(['For the current replication, object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
            load(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat']);
        end
        
        %% TESTING STEP 5
        close all;
        [XTrainLQ,XTrainULQ,XTrainULQ_S,XTestLQ,XTestULQ,XTestULQ_S] = deal({});
        for i=1:size(XTrainL,1)
            XTrainLQ{i,1} = Get_Quad_Features(XTrainL{i,1}); %Getting quadratic X1 with the function Get_Quad_Features.
            XTrainULQ{i,1} = Get_Quad_Features(XTrainUL{i,1});
            XTrainULQ_S{i,1} = Get_Quad_Features(XTrainUL_S{i,1});
        end
        for i=1:size(XTestL,1)
            XTestLQ{i,1} = Get_Quad_Features(XTestL{i,1}); %Getting quadratic XtestL with the function Get_Quad_Features.
            XTestULQ{i,1} = Get_Quad_Features(XTestUL{i,1}); %Getting quadratic XtestUL with the function Get_Quad_Features.
            XTestULQ_S{i,1} = Get_Quad_Features(XTestUL_S{i,1}); %Getting quadratic XtestUL with the function Get_Quad_Features.
        end
        clear CurrentObject;
        Step = 5;
        load(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat']);
        M = CurrentObject.M;
        Q = CurrentObject.Q;
        IforS = CurrentObject.IforS;
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
            CurrentObject.I = ones(size(XTestLQ{1,1},2),size(XTestLQ,1));
            CurrentObject.M = M;
            CurrentObject.Q = Q;
            CurrentObject.IforS = IforS;
            CurrentObject.MakeHistPlot = 'True';
            CurrentObject.GenerateFigure = 'True';
            ModelResults = CurrentObject.RegressModel
            save(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'CurrentObject');
            set(gcf, 'PaperPositionMode', 'auto');
            print(figure(81),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Labeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(82),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Unlabeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
        else
            disp(['For the current replication, object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
            load(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat']);
        end
        
        %% TESTING STEP 6:
        %Call the model defined by the genetic algorithm in step 6 and extract the information regarding which features selected:
        close all;
        clear CurrentObject;
        Step = 6;
        load(['Results/Object_Results/Step',num2str(Step),'_TRY_',num2str(mmm),'.mat']);
        M = CurrentObject.M;
        Q = CurrentObject.Q;
        IforS = CurrentObject.IforS;
        IvecBest1 = CurrentObject.I(:,1);
        IvecBest2 = CurrentObject.I(:,2);
        IvecBest3 = CurrentObject.I(:,3);
        %Do the testing for the model defined in step 6 by genetic algorithm:
        clear CurrentObject;
        if ~exist(['Results/Object_Results/Step',num2str(Step),'_ValidLRep',num2str(TestLRep),'_ValidULRep',num2str(TestULRep),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'file')
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
            CurrentObject.MakeHistPlot = 'True';
            CurrentObject.GenerateFigure = 'True';
            ModelResults = CurrentObject.RegressModel
            save(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat'],'CurrentObject');
            set(gcf, 'PaperPositionMode', 'auto');
            print(figure(81),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Labeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
            print(figure(82),['Results/Figures_FINAL/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'_Sep_Unlabeled_FINAL'],'-dpdf','-fillpage'); %Saving final figure as a *.pdf
        else
            disp(['For the current replication, object for Step',num2str(Step),' already exists. Loading Step',num2str(Step),' object.']);
            load(['Results/Object_Results/Step',num2str(Step),'_TestLRep',num2str(TestLRep),'_TestULRep',num2str(TestULRep),'.mat']);
        end
    end
end
