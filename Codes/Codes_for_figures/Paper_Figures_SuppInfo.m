%% Paper_Figures_SuppInfo
% Author: Mohammad Tanhaemami
% Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.
% Focus is on scenario 2 & scenario 8 (after 2 trials)
% This script creates figures for the supplementary information:
%   - Figure_S2_LinReg_trainvalid
%   - Figure_S3_LinReg_Quad
%   - Figure_S4_Blasi_GBML
%   - Figure_S5_MLPNN
%   - Figure_S6_lab_vs_unlab_features
%   - Figure_S7_Steps_3_4_6
%   - Figure_S8_Linreg_reducedFeats
%   - Figure_S9_GA_linear
%   - Figure_S10_GA_quad
%   - Figure_S11_WeightedModel_TrainValid
%   - Figure_S12_WeightedModel_Test
clear; close all; clc;
warning('off');

%% Call Data
% Call Results_Scenario2 (the latest version):
i = 1;
Object_Scen2 = cell(0);
for s = [1,5,3,4,6] % Calling results in the order they are mentioned in the paper (steps 1,5,3,4,6 in the analysis)
    load(['../../Data_Files/Data_CodesForFigures_Regular_Object_Results/Object_Results/Step',num2str(s),'.mat']);
    CurrentObject.MakeHistPlot = 'False';
    CurrentObject.GenerateFigure = 'False';
    Object_Scen2{i} = CurrentObject;
    i = i + 1;
end

% Call Results_Scenario8_After2Trials (step 6 only)
load('../../Data_Files/Data_CodesForFigures_Weighted_Object_Results/Object_Results/Step6_TRY_2.mat');
CurrentObject.MakeHistPlot = 'False';
CurrentObject.GenerateFigure = 'False';
Object_Scen8_Step6 = CurrentObject;

%% Figures S2, S3, S8, S9, and S10 (old strategy, i.e., the regular approach with the averaged)
% S2: Linear regression
% S3: Linear regression on quadratic features - no feature selection
% S8: Linear regression on reduced features
% S9: Linear regression with GA on linear features
% S10: Linear regression with GA on quadratic features

titles_S238910 = {'Labeled Training Day 1','Labeled Training Day 14','Labeled Training Day 46',...
    'Labeled Training Day 1','Labeled Training Day 14','Labeled Training Day 46',...
    'Labeled Validation Day 0','Labeled Validation Day 15','Labeled Validation Day 37',...
    'Unlabeled Validation Day 0','Unlabeled Validation Day 15','Unlabeled Validation Day 37'};
suptitles_S238910 = {'Fig. S2. Linear regression','Fig. S3. Linear regression on quadratic features',...
    'Fig. S8. Linear regression on reduced features',' Fig. S9. Linear regression with the genetic algorithm on linear features',...
    'Fig. S10. Linear regression with the genetic algorithm on quadratic features'};
FigS = [2,3,8,9,10]; % Supplementary figure orders according to the paper
DaysTrain = [1,14,46];
DaysValid = [0,15,37];
for i = 1:5 % Representing 5 figures for steps 2,3,8,9,10
    orient(figure(i),'landscape');
    for j = 1:3 % Representing the 3 time points, either the training ([2,10,23]) or the validation ([1,11,22])
        % Scatter plots of labeled training data:
        h_S238910(j) = subplot(4,3,j);
        Scen2_YtrL{j,i} = exp(Object_Scen2{1,i}.RegressModel.LogYTrainL{j,1}); % Measured labeled training targets
        Scen2_YprtrL{j,i} = exp(Object_Scen2{1,i}.RegressModel.LogYPredTrainL{j,1}); % Predicted labeled training targets
        scatter(Scen2_YtrL{j,i}(1:1000),Scen2_YprtrL{j,i}(1:1000),'MarkerEdgeColor','none','MarkerFaceColor','flat','MarkerFaceAlpha',0.2); % Replica 1 scatter plot
        Pearson_Corr_Rep1(j) = corr(Scen2_YtrL{j,i}(1:1000),Scen2_YprtrL{j,i}(1:1000)); % Correlation coefficient
        hold on
        scatter(Scen2_YtrL{j,i}(1001:2000),Scen2_YprtrL{j,i}(1001:2000),'MarkerEdgeColor','none','MarkerFaceColor','flat','MarkerFaceAlpha',0.2); % Replica 2 scatter plot
        Pearson_Corr_Rep2(j) = corr(Scen2_YtrL{j,i}(1001:2000),Scen2_YprtrL{j,i}(1001:2000)); % Correlation coefficient
        scatter(Scen2_YtrL{j,i}(2001:3000),Scen2_YprtrL{j,i}(2001:3000),'MarkerEdgeColor','none','MarkerFaceColor','flat','MarkerFaceAlpha',0.2); % Replica 3 scatter plot
        Pearson_Corr_Rep3(j) = corr(Scen2_YtrL{j,i}(2001:3000),Scen2_YprtrL{j,i}(2001:3000)); % Correlation coefficient
        xlim([0,2e6]);
        ylim([0,2e6]);
        Z = get(gca,'XLim');
        plot(Z,Z,'k--');
        text(Z(1,2)-3e5,Z(1,2)-5e5,'m = 1','Rotation',40,'FontSize',18);
        text(1.2e6,7e5,['r_{1} = ',num2str(Pearson_Corr_Rep1(j))],'FontSize',17);
        text(1.2e6,5e5,['r_{2} = ',num2str(Pearson_Corr_Rep2(j))],'FontSize',17);
        text(1.2e6,3e5,['r_{3} = ',num2str(Pearson_Corr_Rep3(j))],'FontSize',17);
        get(gca,'XTickLabel');
        set(gca,'FontSize',17);
        get(gca,'YTickLabel');
        set(gca,'FontSize',17);
        xlabel('Measured','FontSize',18);
        ylabel('Predicted','FontSize',18);
        title(titles_S238910(j),'FontSize',18);
        % Histograms of labeled training data:
        h_S238910(j+3) = subplot(4,3,j+3);
        Scen2_Hist_YtrL{i,j} = hist(Scen2_YtrL{j,i},exp(Object_Scen2{1,i}.BINS))/sum(hist(Scen2_YtrL{j,i},exp(Object_Scen2{1,i}.BINS))); % Histogram of measured labeled training targets
        Scen2_HistYprtrL{i,j} = hist(Scen2_YprtrL{j,i},exp(Object_Scen2{1,i}.BINS))/sum(hist(Scen2_YprtrL{j,i},exp(Object_Scen2{1,i}.BINS))); % Histogram of predicted labeled training targets
        Hist_Train = plot(exp(Object_Scen2{1,i}.BINS)',[Scen2_Hist_YtrL{i,j};Scen2_HistYprtrL{i,j}]'); % plot the above histograms
        hold on
        DKStr{i,j} = Object_Scen2{1,i}.RegressModel.ErrTrainL{j,1}; % Get the KS
        text(1.5e2,0.025,sprintf(['KS_{',num2str(DaysTrain(j)),'} = %1.4f'],eval(num2str(DKStr{i,j}))),'FontSize',16); % Display the KS
        xlim([1e2,1e7]);
        set(gca,'XScale','log');
        ylim([0,0.03]);
        get(gca,'XTickLabel');
        set(gca,'FontSize',17);
        get(gca,'YTickLabel');
        set(gca,'FontSize',17);
        xlabel('Lipid Content (AUC)','FontSize',18);
        ylabel('Probability','FontSize',18);
        title(titles_S238910(j+3),'FontSize',18);
        % Histograms of labeled validation data:
        h_S238910(j+6) = subplot(4,3,j+6);
        Scen2_YvaL{j,i} = exp(Object_Scen2{1,i}.RegressModel.LogYTestL{j,1}); % Measured labeled validation targets
        Scen2_YprvaL{j,i} = exp(Object_Scen2{1,i}.RegressModel.LogYPredTestL{j,1}); % Predicted labeled validation targets
        Scen2_Hist_YvaL{i,j} = hist(Scen2_YvaL{j,i},exp(Object_Scen2{1,i}.BINS))/sum(hist(Scen2_YvaL{j,i},exp(Object_Scen2{1,i}.BINS))); % Histograms for the measured labeled validation targets
        Scen2_HistYprvaL{i,j} = hist(Scen2_YprvaL{j,i},exp(Object_Scen2{1,i}.BINS))/sum(hist(Scen2_YprvaL{j,i},exp(Object_Scen2{1,i}.BINS))); % Histograms for the predicted labeled validation targets
        Hist_ValidL = plot(exp(Object_Scen2{1,i}.BINS)',[Scen2_Hist_YvaL{i,j};Scen2_HistYprvaL{i,j}]'); % Plot the histograms
        hold on
        DKSvaL{i,j} = Object_Scen2{1,i}.RegressModel.ErrTestL{j,1}; % Get the KS
        text(1.5e2,0.025,sprintf(['KS_{',num2str(DaysValid(j)),'} = %1.4f'],eval(num2str(DKSvaL{i,j}))),'FontSize',16); % Display the KS
        xlim([1e2,1e7]);
        set(gca,'XScale','log');
        ylim([0,0.03]);
        get(gca,'XTickLabel');
        set(gca,'FontSize',17);
        get(gca,'YTickLabel');
        set(gca,'FontSize',17);
        xlabel('Lipid Content (AUC)','FontSize',18);
        ylabel('Probability','FontSize',18);
        title(titles_S238910(j+6),'FontSize',18);
        % Histograms of unlabeled validation data:
        Scen2_YprvaUL{j,i} = exp(Object_Scen2{1,i}.RegressModel.LogYPredTestUL{j,1}); % Predicted
        h_S238910(j+9) = subplot(4,3,j+9);
        Scen2_HistYprvaUL{i,j} = hist(Scen2_YprvaUL{j,i},exp(Object_Scen2{1,i}.BINS))/sum(hist(Scen2_YprvaUL{j,i},exp(Object_Scen2{1,i}.BINS)));
        Hist_ValidUL = plot(exp(Object_Scen2{1,i}.BINS)',[Scen2_Hist_YvaL{i,j};Scen2_HistYprvaUL{i,j}]');
        hold on
        DKSvaUL{i,j} = Object_Scen2{1,i}.RegressModel.ErrTestUL{j,1};
        text(1.5e2,0.025,sprintf(['KS_{',num2str(DaysValid(j)),'} = %1.4f'],eval(num2str(DKSvaUL{i,j}))),'FontSize',16);
        xlim([1e2,1e7]);
        set(gca,'XScale','log');
        ylim([0,0.03]);
        get(gca,'XTickLabel');
        set(gca,'FontSize',17);
        get(gca,'YTickLabel');
        set(gca,'FontSize',17);
        xlabel('Lipid Content (AUC)','FontSize',18);
        ylabel('Probability','FontSize',18);
        title(titles_S238910(j+9),'FontSize',18);
    end
    h_S238910_L = legend(h_S238910(1),'Labeled Replica 1','Labeled Replica 2','Labeled Replica 3');
    newPosition = [0.102 0.83 0.2 0.15];
    newUnits = 'normalized';
    set(h_S238910_L,'Position', newPosition,'Units', newUnits);
    suptitle(suptitles_S238910(i));
    saveas(figure(i),['Paper_Figures/suppFigs/Figure_S',num2str(FigS(i)),'.fig']);
    set(gcf, 'PaperPosition', [0 0 19 20]);
    set(gcf, 'PaperSize', [19 20]);
    print(figure(i),['Paper_Figures/suppFigs/Figure_S',num2str(FigS(i))],'-dpdf','-fillpage');
end

%% Figure S4 (GBML by Blasi et al.)
figure(6);
load('../../Data_Files/Data_Thomas_Blasi_Cell_Cycle_Analysis/BLASSI_Results.mat');
H = logspace(log10(min(train_ground_truth)+1),log10(max(train_ground_truth)),100); % Define bins for histograms

N1 = hist(train_ground_truth,H)/sum(hist(train_ground_truth,H)); % Histogram of the training ground truth
N2 = hist(estimated_lipid_content_S_TRAIN,H)/sum(hist(estimated_lipid_content_S_TRAIN,H)); % Histogram of the predicted labeled training data
[~,~,KStrainL] = kstest2(train_ground_truth,estimated_lipid_content_S_TRAIN); % The KS distance
h_S4(1) = subplot(2,2,1);
plot(H,N1,H,N2)
text(1.2e3,0.1,['KS = ',num2str(KStrainL)],'FontSize',10);
text(1.2e3,0.08,['r = ',num2str(corr_S_train(1,2))],'FontSize',10);
xlim([1e3,1e7]);
set(gca,'XScale','log');
legend('Measured','Predicted','Location','southwest');
xlabel('Lipid Content (AUC)');
ylabel('Probability');
title('Training Blassi - Labeled');

N3 = hist(test_ground_truth,H)/sum(hist(test_ground_truth,H)); % Histogram of the testing ground truth
N4 = hist(estimated_lipid_content_S_TEST,H)/sum(hist(estimated_lipid_content_S_TEST,H)); % histogram of the predicted labeled testing data
[~,~,KStestL] = kstest2(test_ground_truth,estimated_lipid_content_S_TEST); % The KS distance
h_S4(2) = subplot(2,2,2);
plot(H,N3,H,N4)
text(1.2e3,0.05,['KS = ',num2str(KStestL)],'FontSize',10);
text(1.2e3,0.04,['r = ',num2str(corr_S_test(1,2))],'FontSize',10);
xlim([1e3,1e7]);
set(gca,'XScale','log');
xlabel('Lipid Content (AUC)');
ylabel('Probability');
title('Testing Blassi - Labeled');


N1 = hist(train_ground_truth,H)/sum(hist(train_ground_truth,H)); % Histogram of the training ground truth
NU2 = hist(estimated_lipid_content_U_TRAIN,H)/sum(hist(estimated_lipid_content_U_TRAIN,H)); % Histogram of the predicted unlabeled training data
[~,~,KStrainU] = kstest2(train_ground_truth,estimated_lipid_content_U_TRAIN);
h_S4(3) = subplot(2,2,3);
plot(H,N1,H,NU2)
text(1.2e3,0.22,['KS = ',num2str(KStrainU)],'FontSize',10);
text(1.2e3,0.18,['r = ',num2str(corr_U_train(1,2))],'FontSize',10);
xlim([1e3,1e7]);
set(gca,'XScale','log');
xlabel('Lipid Content (AUC)');
ylabel('Probability');
title('Training Blassi - Unlabeled');

N3 = hist(test_ground_truth,H)/sum(hist(test_ground_truth,H)); % Histogram of the testing ground truth
NU4 = hist(estimated_lipid_content_U_TEST,H)/sum(hist(estimated_lipid_content_U_TEST,H)); % histogram of the predicted unlabeled testing data
[~,~,KStestU] = kstest2(test_ground_truth,estimated_lipid_content_U_TEST);
h_S4(4) = subplot(2,2,4);
plot(H,N3,H,NU4)
text(1.2e3,0.18,['KS = ',num2str(KStestU)],'FontSize',10);
text(1.2e3,0.15,['r = ',num2str(corr_U_test(1,2))],'FontSize',10);
xlim([1e3,1e7]);
set(gca,'XScale','log');
xlabel('Lipid Content (AUC)');
ylabel('Probability');
title('Testing Blassi - Unlabeled');

suptitle('Fig. S4. Gradiet boosting machine learning by Blasi et al.');

saveas(figure(6),'Paper_Figures/suppFigs/Figure_S4.fig');
set(gcf, 'PaperPositionMode', 'auto');
print(figure(6),'Paper_Figures/suppFigs/Figure_S4','-dpdf','-fillpage');

%% Figure S5 (MLPNN)
figure(7);
% Get the results of the MLPNN:
File_Name = '../../Data_Files/Data_MLPNN/AllPico.csv';
FileContent = importdata(File_Name);
NNmodelDATA = FileContent.data;
NN_YPredTrain = NNmodelDATA(:,1); % Predicted training data
NN_YTrain = NNmodelDATA(:,2); % Trainging ground truth
NN_YPredValidL = NNmodelDATA(:,3); % Predicted validation data (labeled)
NN_YValidL = NNmodelDATA(:,4); % Validation ground truth
NN_YPredValidUL = NNmodelDATA(:,5); %Predicted validation data (unlabeled)
NN_BINS = logspace(log10(max(NN_YValidL))-4,log10(max(NN_YValidL)),500); %Set logarithmic range for x-axis(500 x points).

% Training histograms:
h_S5(1) = subplot(2,2,1);
NN_histYTrain = hist(NN_YTrain,NN_BINS)/sum(hist(NN_YTrain,NN_BINS)); % Measured histograms
NN_histYPredTrain = hist(NN_YPredTrain,NN_BINS)/sum(hist(NN_YPredTrain,NN_BINS)); % Predicted histograms
plot_histNNtrain = plot(NN_BINS',[NN_histYTrain;NN_histYPredTrain]'); % Plot histograms
[~,~,KS_NNtrain] = kstest2(NN_histYTrain,NN_histYPredTrain); % KS distance
set(gca,'XScale','log');
xlim([1e3,1e7]);
xlabel('Lipid Content (AUC)');
ylim([0,0.017]);
ylabel('Probability');
text(1.2e3,0.014,sprintf('KS = %1.4f',eval(num2str(KS_NNtrain))),'FontSize',9); % Display KS on plot
title('Training the MLPNN');

% Validation histograms:
NN_HistYValidL = hist(NN_YValidL,NN_BINS)/sum(hist(NN_YValidL,NN_BINS)); % Measured labeled histograms
NN_HistYPredValidL = hist(NN_YPredValidL,NN_BINS)/sum(hist(NN_YPredValidL,NN_BINS)); % Predicted labeled histograms
NN_HistYPredValidUL = hist(NN_YPredValidUL,NN_BINS)/sum(hist(NN_YPredValidUL,NN_BINS)); % Predicted unlabeled histograms
h_S5(2) = subplot(2,2,3);
plotHist_NN_ValidL = plot(NN_BINS',[NN_HistYValidL;NN_HistYPredValidL]'); % Plot histograms
[~,~,KS_NN_ValidL] = kstest2(NN_HistYValidL,NN_HistYPredValidL); % KS distance
set(gca,'XScale','log');
xlim([1e3,1e7]);
xlabel('Lipid Content (AUC)');
ylim([0,0.017]);
ylabel('Probability');
text(1.2e3,0.014,sprintf('KS = %1.4f',eval(num2str(KS_NN_ValidL))),'FontSize',9); % Display KS on plot
title('Validation of the MLPNN - Labeled');
h_S5(3) = subplot(2,2,4);
plotHist_NN_ValidUL = plot(NN_BINS',[NN_HistYValidL;NN_HistYPredValidUL]');
[~,~,KS_NN_ValidUL] = kstest2(NN_HistYValidL,NN_HistYPredValidUL);
set(gca,'XScale','log');
xlim([1e3,1e7]);
xlabel('Lipid Content (AUC)');
ylim([0,0.017]);
ylabel('Probability');
text(1.2e3,0.014,sprintf('KS = %1.4f',eval(num2str(KS_NN_ValidUL))),'FontSize',9);
title('Validation of the MLPNN - Unlabeled');

suptitle('Fig. S5. Multilayer perceptron neural network');
saveas(figure(7),'Paper_Figures/suppFigs/Figure_S5.fig');
set(gcf, 'PaperPositionMode', 'auto');
print(figure(7),'Paper_Figures/suppFigs/Figure_S5','-dpdf','-fillpage');

%% Figure S6: FIND POSSIBLE CORRUPTIONS -- Step 2 in the main script which was commented out there
load('../../Data_Files/NewPico_Data_NoZero.mat'); %Load the file Pico_Data.mat.
TargetSet = [3,9]; %Enter the target parmeters' index (here 3 and 9 are targets).
TimePoint = 23; %Enter number of time points (days after nitorgen starvation) to be selected.
NCells = 3000;  %Enter number of cells to be used from all replication
FeatureInds = 1:size(DATA.FeatureNames(~strcmp(DATA.FeatureNames,'Time')),2); %Define Feature_Inds
FeatureInds(TargetSet) = []; %Remove indicies for tragets.
SelectedTimePoint = 6; %We have 23 time points. We select one of them, arbitrarily.
ValidLRep = 3;
ValidULRep = 3;
[Stained,Unstained,KS] = Paper_Corruption_Test(NCells,FeatureInds,SelectedTimePoint,ValidLRep,ValidULRep,DATA);
saveas(figure(8),'Paper_Figures/suppFigs/Figure_S6.fig');
saveas(figure(8),'Paper_Figures/suppFigs/Figure_S6.pdf');

%% Figure S11: FINAL method: training and validation the weighted model
orient(figure(9),'landscape');
titles_S11 = {'Labeled Training Day 1','Labeled Training Day 14','Labeled Training Day 46',...
    'Labeled Training Day 1','Labeled Training Day 14','Labeled Training Day 46',...
    'Labeled Validation Day 0','Labeled Validation Day 15','Labeled Validation Day 37',...
    'Unlabeled Validation Day 0','Unlabeled Validation Day 15','Unlabeled Validation Day 37'};
DaysTrain = [1,14,46];
DaysValid = [0,15,37];
for j = 1:3
   h_S11(j) = subplot(4,3,j);
   % Scatter plots of labeled training data:
   Scen8_YtrL{j,1} = exp(Object_Scen8_Step6.RegressModel.LogYTrainL{j,1}); % Measured
   Scen8_YprtrL{j,1} = exp(Object_Scen8_Step6.RegressModel.LogYPredTrainL{j,1}); % Labeled predicted
   scatter(Scen8_YtrL{j,1}(1:1000),Scen8_YprtrL{j,1}(1:1000),'MarkerEdgeColor','none','MarkerFaceColor','flat','MarkerFaceAlpha',0.2); % Replica 1
   Pearson_Corr_Rep1(j) = corr(Scen8_YtrL{j,1}(1:1000),Scen8_YprtrL{j,1}(1:1000));
   hold on
   scatter(Scen8_YtrL{j,1}(1001:2000),Scen8_YprtrL{j,1}(1001:2000),'MarkerEdgeColor','none','MarkerFaceColor','flat','MarkerFaceAlpha',0.2); % Replica 2
   Pearson_Corr_Rep2(j) = corr(Scen8_YtrL{j,1}(1001:2000),Scen8_YprtrL{j,1}(1001:2000));
   scatter(Scen8_YtrL{j,1}(2001:3000),Scen8_YprtrL{j,1}(2001:3000),'MarkerEdgeColor','none','MarkerFaceColor','flat','MarkerFaceAlpha',0.2); % Replica 3
   Pearson_Corr_Rep3(j) = corr(Scen8_YtrL{j,1}(2001:3000),Scen8_YprtrL{j,1}(2001:3000));
   xlim([0,2.5e6]);
   ylim([0,2.5e6]);
   Z = get(gca,'XLim');
   plot(Z,Z,'k--');
   text(Z(1,2)-4e5,Z(1,2)-6e5,'m = 1','Rotation',41,'FontSize',18);
   text(1.5e6,7e5,['r_{1} = ',num2str(Pearson_Corr_Rep1(j))],'FontSize',17);
   text(1.5e6,5e5,['r_{2} = ',num2str(Pearson_Corr_Rep2(j))],'FontSize',17);
   text(1.5e6,3e5,['r_{3} = ',num2str(Pearson_Corr_Rep3(j))],'FontSize',17);
   get(gca,'XTickLabel');
   set(gca,'FontSize',17);
   get(gca,'YTickLabel');
   set(gca,'FontSize',17);
   xlabel('Measured','FontSize',18);
   ylabel('Predicted','FontSize',18);
   title(titles_S11(j),'FontSize',18);
   % Histograms of labeled training data:
   h_S11(j+3) = subplot(4,3,j+3);
   Scen8_Step6_Hist_YtrL{1,j} = hist(Scen8_YtrL{j,1},exp(Object_Scen8_Step6.BINS))/sum(hist(Scen8_YtrL{j,1},exp(Object_Scen8_Step6.BINS))); % Measured
   Scen8_Step8_HistYprtrL{1,j} = hist(Scen8_YprtrL{j,1},exp(Object_Scen8_Step6.BINS))/sum(hist(Scen8_YprtrL{j,1},exp(Object_Scen8_Step6.BINS))); % Labeled predicted
   Hist_Scen8_Step6_Train = plot(exp(Object_Scen8_Step6.BINS)',[Scen8_Step6_Hist_YtrL{1,j};Scen8_Step8_HistYprtrL{1,j}]');
   hold on
   DKStr_Scen8_Step6(j) = Object_Scen8_Step6.RegressModel.ErrTrainL(j); % KS distance
   text(1.5e2,0.1,sprintf(['KS_{',num2str(DaysTrain(j)),'} = %1.4f'],eval(num2str(DKStr_Scen8_Step6(j)))),'FontSize',16); % Display KS on plot
   xlim([1e2,1e7]);
   set(gca,'XScale','log');
   ylim([0,0.12]);
   get(gca,'XTickLabel');
   set(gca,'FontSize',17);
   get(gca,'YTickLabel');
   set(gca,'FontSize',17);
   xlabel('Lipid Content (AUC)','FontSize',18);
   ylabel('Probability','FontSize',18);
   title(titles_S11(j+3),'FontSize',18);
   % Histograms of labeled validation data:
   h_S11(j+6) = subplot(4,3,j+6);
   Scen8_YvaL{j,1} = exp(Object_Scen8_Step6.RegressModel.LogYTestL{j,1}); % Measured
   Scen8_YprvaL{j,1} = exp(Object_Scen8_Step6.RegressModel.LogYPredTestL{j,1}); % Labeled predicted
   Scen8_Step6_Hist_YvaL{1,j} = hist(Scen8_YvaL{j,1},exp(Object_Scen8_Step6.BINS))/sum(hist(Scen8_YvaL{j,1},exp(Object_Scen8_Step6.BINS))); % Measured histogram
   Scen8_Step6_HistYprvaL{1,j} = hist(Scen8_YprvaL{j,1},exp(Object_Scen8_Step6.BINS))/sum(hist(Scen8_YprvaL{j,1},exp(Object_Scen8_Step6.BINS))); % Labeled predicted histogram
   Hist_Scen8_Step6_ValidL = plot(exp(Object_Scen8_Step6.BINS)',[Scen8_Step6_Hist_YvaL{1,j};Scen8_Step6_HistYprvaL{1,j}]');
   hold on
   DKSvaL_Scen8_Step6(j) = Object_Scen8_Step6.RegressModel.ErrTestL(j); % KS distance
   text(1.5e2,0.1,sprintf(['KS_{',num2str(DaysValid(j)),'} = %1.4f'],eval(num2str(DKSvaL_Scen8_Step6(j)))),'FontSize',16); % Display KS on plot
   xlim([1e2,1e7]);
   set(gca,'XScale','log');
   ylim([0,0.12]);
   get(gca,'XTickLabel');
   set(gca,'FontSize',17);
   get(gca,'YTickLabel');
   set(gca,'FontSize',17);
   xlabel('Lipid Content (AUC)','FontSize',18);
   ylabel('Probability','FontSize',18);
   title(titles_S11(j+6),'FontSize',18);
   % Histograms of unlabeled validation data:
   h_S11(j+9) = subplot(4,3,j+9);
   Scen8_YprvaUL{j,1} = exp(Object_Scen8_Step6.RegressModel.LogYPredTestUL{j,1}); % Unlabeled predicted
   Scen8_Step6_HistYprvaUL{1,j} = hist(Scen8_YprvaUL{j,1},exp(Object_Scen8_Step6.BINS))/sum(hist(Scen8_YprvaUL{j,1},exp(Object_Scen8_Step6.BINS))); % Unlabeled predicted histogram
   Hist_Scen8_Step6_ValidUL = plot(exp(Object_Scen8_Step6.BINS)',[Scen8_Step6_Hist_YvaL{1,j};Scen8_Step6_HistYprvaUL{1,j}]');
   hold on
   DKSvaUL_Scen8_Step6(j) = Object_Scen8_Step6.RegressModel.ErrTestUL(j); % KS distance
   text(1.5e2,0.1,sprintf(['KS_{',num2str(DaysValid(j)),'} = %1.4f'],eval(num2str(DKSvaUL_Scen8_Step6(j)))),'FontSize',16); % Display KS on plot
   xlim([1e2,1e7]);
   set(gca,'XScale','log');
   ylim([0,0.12]);
   get(gca,'XTickLabel');
   set(gca,'FontSize',17);
   get(gca,'YTickLabel');
   set(gca,'FontSize',17);
   xlabel('Lipid Content (AUC)','FontSize',18);
   ylabel('Probability','FontSize',18);
   title(titles_S11(j+9),'FontSize',18);
end
h_S8_L = legend(h_S11(1),'Labeled Replica 1','Labeled Replica 2','Labeled Replica 3');
newPosition = [0.102 0.83 0.2 0.15];
newUnits = 'normalized';
set(h_S8_L,'Position', newPosition,'Units', newUnits);
suptitle('Fig. S11. Training and validation of the proposed strategy with weighted model');
saveas(figure(9),'Paper_Figures/suppFigs/Figure_S11.fig');
set(gcf, 'PaperPosition', [0 0 19 20]);
set(gcf, 'PaperSize', [19 20]);
print(figure(9),'Paper_Figures/suppFigs/Figure_S11','-dpdf','-fillpage');

%% Figure S12: Testing the final method on all testing days
orient(figure(10),'landscape');
Object_TestReps_Scen8 = cell(0);
Scen8_HistYtesL = cell(0);
Scen8_HistYprtesUL = cell(0);
TestDays = [5,6,7,9,10,11,12,16,19,20,21,22,23,26,27,30,34]; % Testing days
for t = 1:17 % 17 testing time points at each replication combination
    for i = 1:3 % 3 labeled replicates
        for j = 1:4 % 4 unlabeled replicates
            load(['../../Data_Files/Data_CodesForFigures_Weighted_Object_Results/Object_Results/Step6_TestLRep',num2str(i),'_TestULRep',num2str(j),'.mat']);
            CurrentObject.MakeHistPlot = 'False'; % Prevent object from plotting histograms
            CurrentObject.GenerateFigure = 'False'; % Prevent object from generatig figures
            Object_TestReps_Scen8{i,j,t} = CurrentObject;
            Scen8_YtsL{i,j,t} = exp(Object_TestReps_Scen8{i,j,t}.RegressModel.LogYTestL{t,1}); % Measured testing
            Scen8_YprtsUL{i,j,t} = exp(Object_TestReps_Scen8{i,j,t}.RegressModel.LogYPredTestUL{t,1}); % Predicted unlabeled testing
            Scen8_HistYtesL{i,j,t} = hist(Scen8_YtsL{i,j,t},exp(Object_TestReps_Scen8{i,j,t}.BINS))/sum(hist(Scen8_YtsL{i,j,t},exp(Object_TestReps_Scen8{i,j,t}.BINS))); % Measured targets (labeled) histograms
            Scen8_HistYprtesUL{i,j,t} = hist(Scen8_YprtsUL{i,j,t},exp(Object_TestReps_Scen8{i,j,t}.BINS))/sum(hist(Scen8_YprtsUL{i,j,t},exp(Object_TestReps_Scen8{i,j,t}.BINS))); % Predicted (unlabeled) histograms
            HistPlotsYtesL_Scen8{i,j} = plot3(4*t*ones(1,500),exp(Object_TestReps_Scen8{i,j,t}.BINS)',[Scen8_HistYtesL{i,j,t}]','b','LineWidth',1); % Measured 3D plot
            hold on
            HistPlotsYprtesUL_Scen8{i,j} = plot3(4*t*ones(1,500),exp(Object_TestReps_Scen8{i,j,t}.BINS)',[Scen8_HistYprtesUL{i,j,t}]','r','LineWidth',1); % Predicted 3D plot
            KStestScen8(i,j,t) = Object_TestReps_Scen8{i,j,t}.RegressModel.ErrTestUL(t);
        end
    end
    KS_avg(t) = mean(mean(KStestScen8(:,:,t))); % Average KS distance for replications
    text(4*t,1.5e6,0.083,sprintf(['KS_{',num2str(TestDays(t)),'} = %1.4f'],eval(num2str(KS_avg(t)))),'FontSize',24); % Display KS on plots
end
xticks(4:4:68);
xticklabels(1:17);
get(gca,'XTickLabel');
set(gca,'FontSize',20);
get(gca,'YTickLabel');
set(gca,'FontSize',20);
xlabel('Time Point','FontSize',24,'Rotation',-62);
ylim([1e2,1e7]);
set(gca,'YScale','log');
ylabel('Lipid Content (AUC)','FontSize',24,'Rotation',32);
zlabel('Probability','FontSize',24);
title('Fig. S12. Testing the proposed weighted model for all 17 time points','FontSize',24);
grid on
view(60,75);

saveas(figure(10),'Paper_Figures/suppFigs/Figure_S12.fig');
set(gcf, 'PaperPosition', [0 0 15 20]);
set(gcf, 'PaperSize', [15 20]);
print(figure(10),'Paper_Figures/suppFigs/Figure_S12','-dpdf','-fillpage');

%% Figure S7: 3D representation for histograms of labeled training data vs. unlabeled validation data for steps 3, 4, 6 for the regular approach
% Step 3: Linear regression with reduced features
% Step 4: Linear regression with the genetic algorithm for automated feature selection on linear features
% Step 6: Linear regression with the genetic algorithm for automated feature selection on quadratic features
% 3 timepoints
% 2 colors (measured & predicted)
clear; close all; clc;
i = 1;
Object_Scen2 = cell(0);
for s = [1,3,4,6]
    load(['../../Data_Files/Data_CodesForFigures_Regular_Object_Results/Object_Results/Step',num2str(s),'.mat']);
    CurrentObject.MakeHistPlot = 'False';
    CurrentObject.GenerateFigure = 'False';
    Object_Scen2{i} = CurrentObject;
    i = i + 1;
end
figure(11);
titles = {'(A) Labeled Training with Reduced Features','(B) Unlabeled Validation with Reduced Features','(C) Labeled Training with the GA on Linear Features',...
    '(D) Unlabeled Validation with the GA on Linear Features','(E) Labeled Training with the GA on Quadratic Features','(F) Unlabeled Validation with the GA on Quadratic Features'};
DaysTrain = [1,14,46];
for i = 2:4 % Steps 3, 4, and 6
    for j = 1:3 % Time points
        h_S7(2*i-3) = subplot(3,2,2*i-3);
        Scen2_YtrL{j,i} = exp(Object_Scen2{1,i}.RegressModel.LogYTrainL{j,1}); % Measured training
        Scen2_YprtrL{j,i} = exp(Object_Scen2{1,i}.RegressModel.LogYPredTrainL{j,1}); % Labeled predicted training
        Scen2_346_HistYtrL{i-1,j} = hist(Scen2_YtrL{j,i},exp(Object_Scen2{1,i}.BINS))/sum(hist(Scen2_YtrL{j,i},exp(Object_Scen2{1,i}.BINS)));
        Scen2_346_HistYprtrL{i-1,j} = hist(Scen2_YprtrL{j,i},exp(Object_Scen2{1,i}.BINS))/sum(hist(Scen2_YprtrL{j,i},exp(Object_Scen2{1,i}.BINS)));
        HistPlots_346_YtrL = plot3(j*ones(1,500),exp(Object_Scen2{1,i}.BINS)',[Scen2_346_HistYtrL{i-1,j}]','b','LineWidth',2);
        hold on
        HistPlots_346_YprtrL = plot3(j*ones(1,500),exp(Object_Scen2{1,i}.BINS)',[Scen2_346_HistYprtrL{i-1,j}]','r','LineWidth',2);
        DKStr_Scen2_346{i-1,j} = Object_Scen2{1,i}.RegressModel.ErrTrainL{j,1};
        text(j,1.41e6,0.025,sprintf(['KS_{',num2str(DaysTrain(j)),'} = %1.4f'],eval(num2str(DKStr_Scen2_346{i-1,j}))),'FontSize',16);
    end
    xticklabels([1,2,3]);
    get(gca,'XTickLabel');
    set(gca,'FontSize',16);
    get(gca,'YTickLabel');
    set(gca,'FontSize',16);
    xlabel('Time Point','FontSize',16,'Rotation',-51);
    ylim([1e2,1e7]);
    set(gca,'YScale','log');
    ylabel('Lipid Content (AUC)','FontSize',16,'Rotation',8);
    zlabel('Probability','FontSize',16);
    title(titles{2*i-3},'FontSize',16);
    grid on
    view(70,30);
end
legend(h_S7(1),{'Measured','Label-free'},'Location','northwest','FontSize',15);
DaysValid = [0,15,37];
for i = 2:4
    for j = 1:3
        h_S7(2*i-2) = subplot(3,2,2*i-2);
        Scen2_YvaL{j,i} = exp(Object_Scen2{1,i}.RegressModel.LogYTestL{j,1});
        Scen2_YprvaUL{j,i} = exp(Object_Scen2{1,i}.RegressModel.LogYPredTestUL{j,1});
        Scen2_346_HistYvaL{i-1,j} = hist(Scen2_YvaL{j,i},exp(Object_Scen2{1,i}.BINS))/sum(hist(Scen2_YvaL{j,i},exp(Object_Scen2{1,i}.BINS)));
        Scen2_346_HistYprvaUL{i-1,j} = hist(Scen2_YprvaUL{j,i},exp(Object_Scen2{1,i}.BINS))/sum(hist(Scen2_YprvaUL{j,i},exp(Object_Scen2{1,i}.BINS)));
        HistPlots_346_YvaL = plot3(j*ones(1,500),exp(Object_Scen2{1,i}.BINS)',[Scen2_346_HistYvaL{i-1,j}]','b','LineWidth',2);
        hold on
        HistPlots_346_YprvaUL = plot3(j*ones(1,500),exp(Object_Scen2{1,i}.BINS)',[Scen2_346_HistYprvaUL{i-1,j}]','r','LineWidth',2);
        DKSva_Scen2_346{i-1,j} = Object_Scen2{1,i}.RegressModel.ErrTestUL{j,1};
        text(j,1.41e6,0.025,sprintf(['KS_{',num2str(DaysValid(j)),'} = %1.4f'],eval(num2str(DKSva_Scen2_346{i-1,j}))),'FontSize',16);
    end
    xticklabels([1,2,3]);
    get(gca,'XTickLabel');
    set(gca,'FontSize',16);
    get(gca,'YTickLabel');
    set(gca,'FontSize',16);
    xlabel('Time Point','FontSize',16,'Rotation',-51);
    ylim([1e2,1e7]);
    set(gca,'YScale','log');
    ylabel('Lipid Content (AUC)','FontSize',16,'Rotation',8);
    zlabel('Probability','FontSize',16);
    title(titles{2*i-2},'FontSize',16);
    grid on
    view(70,30);
end
suptitle('Further Analyses with Reduced Features, and the Genetic Algorithm (Linear and Quadratic)');
saveas(figure(11),'Paper_Figures/suppFigs/Figure_S7_Steps_3_4_6.fig');
set(gcf, 'PaperPosition', [0 0 19 22]);
set(gcf, 'PaperSize', [19 22]);
print(figure(11),'Paper_Figures/suppFigs/Figure_S7_Steps_3_4_6','-dpdf','-fillpage');
