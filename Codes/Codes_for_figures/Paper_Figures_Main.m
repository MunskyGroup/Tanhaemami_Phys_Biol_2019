%% Paper_Figures_Main
% Author: Mohammad Tanhaemami
% Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.
% This scripts calls the final results and plots them in figures that for the paper.
% This script calls selected features at each step (1,3,4,6) and stores them in an Excel spreadsheet.
% This script also counts the number of each selected feature / feature combination. 
% Generates figures 3, 4, 5, 6, and 7 in the paper
clear; close all; clc;
path(path,'../Functions');

%% Call Data
% Call Results_regular_approach:
i = 1;
Object_Scen2 = cell(0);
for s = [1,3,4,6]
    load(['../../Data_Files/Data_CodesForFigures_Regular_Object_Results/Object_Results/Step',num2str(s),'.mat']);
    CurrentObject.MakeHistPlot = 'False';
    CurrentObject.GenerateFigure = 'False';
    Object_Scen2{i} = CurrentObject;
    i = i + 1;
end

% Call Results_weighted_model:
i = 1;
Object_Scen8 = cell(0);
for s = [1,3,4,6]
    load(['../../Data_Files/Data_CodesForFigures_Weighted_Object_Results/Object_Results/Step',num2str(s),'_TRY_2.mat']);
    CurrentObject.MakeHistPlot = 'False';
    CurrentObject.GenerateFigure = 'False';
    Object_Scen8{i} = CurrentObject;
    i = i + 1;
end

%% Figure 3: Linear regression scatter plot & 3D representation for histograms (Panels 1&2 have moved to SI in the paper)
% This figure contains the following panels (3 time points in each panel):
% Panel 1: Scatter plot of labeled training data (3 colors for 3 time points) 
% Panel 2: Histogram of labeled training data (2 colors for measured and predicted) 
% Panel 3: Histogram of labeled validation data (2 colors for measured and predicted)
% Panel 4: Histogram of unlabeled validation data (2 colors for measured and predicted)

orient(figure(3),'landscape');
titles = {'(a) Scatter Plot of Labeled Training Data','(b) Labeled Training','(c) Labeled Validation','(d) Unlabeled Validation'};
% Panel 1:
h1(1) = subplot(2,2,1);
for j = 1:3 % Get the results for the regular approach (averaged model)
    Scen2_YtrL{j,1} = exp(Object_Scen2{1,1}.RegressModel.LogYTrainL{j,1}); % Measured labeled training targets
    Scen2_YprtrL{j,1} = exp(Object_Scen2{1,1}.RegressModel.LogYPredTrainL{j,1}); % Predicted labeled training targets
    scatter(Scen2_YtrL{j,1},Scen2_YprtrL{j,1},'MarkerEdgeColor','none','MarkerFaceColor','flat','MarkerFaceAlpha',0.2); % Scatter plot
    hold on
    Pearson_Corr(j) = corr(Scen2_YtrL{j,1},Scen2_YprtrL{j,1}); % Correlation coefficient
end
xlim([0,2.5e6]);
ylim([0,2.5e6]);
Z = get(gca,'XLim');
plot(Z,Z,'k--','LineWidth',1);
text(1.65e6,1.45e6,'m = 1','Rotation',47,'FontSize',16);
legend('off');
text(1e5,2.3e6,['r_{1} = ',num2str(Pearson_Corr(1))],'FontSize',17);
text(1e5,2.1e6,['r_{14} = ',num2str(Pearson_Corr(2))],'FontSize',17);
text(1e5,1.9e6,['r_{46} = ',num2str(Pearson_Corr(3))],'FontSize',17);
get(gca,'XTickLabel');
set(gca,'FontSize',17);
get(gca,'YTickLabel');
set(gca,'FontSize',17);
xlabel('Measured','FontSize',17);
ylabel('Predicted','FontSize',17);
legend(h1(1),{'Day 1', 'Day 14', 'Day 46'},'Location','southeast','FontSize',16);
title(titles(1),'FontSize',17);
% Panel 2:
h1(2) = subplot(2,2,2);
DaysTrain = [1,14,46];
for j = 1:3
    Scen2_LinReg_HistYtrL{1,j} = hist(Scen2_YtrL{j,1},exp(Object_Scen2{1,1}.BINS))/sum(hist(Scen2_YtrL{j,1},exp(Object_Scen2{1,1}.BINS))); % Histogram of the measured labeled training targets
    Scen2_LinReg_HistYprtrL{1,j} = hist(Scen2_YprtrL{j,1},exp(Object_Scen2{1,1}.BINS))/sum(hist(Scen2_YprtrL{j,1},exp(Object_Scen2{1,1}.BINS))); % Histogram of the predicted labeled training targets
    HistPlots_LinReg_YtrL = plot3(j*ones(1,500),exp(Object_Scen2{1,1}.BINS)',[Scen2_LinReg_HistYtrL{1,j}]','b','LineWidth',2); % 3D plotting the histograms
    hold on
    HistPlots_LinReg_YprtrL = plot3(j*ones(1,500),exp(Object_Scen2{1,1}.BINS)',[Scen2_LinReg_HistYprtrL{1,j}]','r','LineWidth',2);
    DKStr_Scen2_LinReg{1,j} = Object_Scen2{1,1}.RegressModel.ErrTrainL{j,1}; % Get the KS distances
    text(j,1.41e6,0.025,sprintf(['KS_{',num2str(DaysTrain(j)),'} = %1.4f'],eval(num2str(DKStr_Scen2_LinReg{1,j}))),'FontSize',17); % Display the KS on plots
end
legend(h1(2),{'Measured','Estimated from non-label channels'},'Location','northwest','FontSize',16);
get(gca,'XTickLabel');
set(gca,'FontSize',17);
get(gca,'YTickLabel');
set(gca,'FontSize',17);
xlabel('Time Point','FontSize',17,'Rotation',-49);
ylim([1e2,1e7]);
set(gca,'YScale','log');
ylabel('Lipid Content (AUC)','FontSize',17,'Rotation',7);
zlabel('Probability','FontSize',17);
title(titles(2),'FontSize',17);
grid on
view(70,30);
% Panel 3:
h1(3) = subplot(2,2,3);
DaysValid = [0,15,37];
for j = 1:3 % Get the validation results
    Scen2_YvaL{j,1} = exp(Object_Scen2{1,1}.RegressModel.LogYTestL{j,1}); % Measured labeled validaiton targets
    Scen2_YprvalL{j,1} = exp(Object_Scen2{1,1}.RegressModel.LogYPredTestL{j,1}); % Predicted labeled validation targets
    Scen2_LinReg_HistYvaL{1,j} = hist(Scen2_YvaL{j,1},exp(Object_Scen2{1,1}.BINS))/sum(hist(Scen2_YvaL{j,1},exp(Object_Scen2{1,1}.BINS))); % Histogram of the measured labeled validation targets
    Scen2_LinReg_HistYprvaL{1,j} = hist(Scen2_YprvalL{j,1},exp(Object_Scen2{1,1}.BINS))/sum(hist(Scen2_YprvalL{j,1},exp(Object_Scen2{1,1}.BINS))); % Histogram of the predicted labeled validation targets
    HistPlots_LinReg_YvaL = plot3(j*ones(1,500),exp(Object_Scen2{1,1}.BINS)',[Scen2_LinReg_HistYvaL{1,j}]','b','LineWidth',2); % 3D plot histograms
    hold on
    HistPlots_LinReg_YprvaL = plot3(j*ones(1,500),exp(Object_Scen2{1,1}.BINS)',[Scen2_LinReg_HistYprvaL{1,j}]','r','LineWidth',2); % 3D plot histograms
    DKSvaL_Scen2_LinReg{1,j} = Object_Scen2{1,1}.RegressModel.ErrTestL{j,1}; % Get the KS distance result
    text(j,1.41e6,0.025,sprintf(['KS_{',num2str(DaysValid(j)),'} = %1.4f'],eval(num2str(DKSvaL_Scen2_LinReg{1,j}))),'FontSize',17); % Dosplay the KS on plots
end
get(gca,'XTickLabel');
set(gca,'FontSize',17);
get(gca,'YTickLabel');
set(gca,'FontSize',17);
xlabel('Time Point','FontSize',17,'Rotation',-49);
ylim([1e2,1e7]);
set(gca,'YScale','log');
ylabel('Lipid Content (AUC)','FontSize',17,'Rotation',7);
zlabel('Probability','FontSize',17);
title(titles(3),'FontSize',17);
grid on
view(70,30);
% Panel 4:
h1(4) = subplot(2,2,4);
for j = 1:3 % Get the unlabeled validation results
    Scen2_YprvaUL{j,1} = exp(Object_Scen2{1,1}.RegressModel.LogYPredTestUL{j,1}); % Predicted unlabeled validation targets
    Scen2_LinReg_HistYprvaUL{1,j} = hist(Scen2_YprvaUL{j,1},exp(Object_Scen2{1,1}.BINS))/sum(hist(Scen2_YprvaUL{j,1},exp(Object_Scen2{1,1}.BINS))); % Histogram of the predicted unlabeled validation targets
    HistPlots_LinReg_YvaL = plot3(j*ones(1,500),exp(Object_Scen2{1,1}.BINS)',[Scen2_LinReg_HistYvaL{1,j}]','b','LineWidth',2); % 3D plot the histograms of measured laebeld validation targets
    hold on
    HistPlots_LinReg_YprvaUL = plot3(j*ones(1,500),exp(Object_Scen2{1,1}.BINS)',[Scen2_LinReg_HistYprvaUL{1,j}]','r','LineWidth',2);% 3D plot the histogram of the predicted unlabeled validation targets
    DKSvaUL_Scen2_LinReg{1,j} = Object_Scen2{1,1}.RegressModel.ErrTestUL{j,1}; % Get the KS distance
    text(j,1.41e6,0.025,sprintf(['KS_{',num2str(DaysValid(j)),'} = %1.4f'],eval(num2str(DKSvaUL_Scen2_LinReg{1,j}))),'FontSize',17); % Display KS on plot
end
get(gca,'XTickLabel');
set(gca,'FontSize',17);
get(gca,'YTickLabel');
set(gca,'FontSize',17);
xlabel('Time Point','FontSize',17,'Rotation',-49);
ylim([1e2,1e7]);
set(gca,'YScale','log');
ylabel('Lipid Content (AUC)','FontSize',17,'Rotation',7);
zlabel('Probability','FontSize',17);
title(titles(4),'FontSize',17);
grid on
view(70,30);
saveas(figure(3),'Paper_Figures/mainFigs/Figure_3_Linear_Regression.fig');
set(gcf, 'PaperPosition', [0 0 18 15]);
set(gcf, 'PaperSize', [18 15]);
print(figure(3),'Paper_Figures/mainFigs/Figure_3_Linear_Regression','-dpdf','-fillpage');

%% Figure 4: Results of the FINAL method (weighted model)
% This figure contains the following panels from scnenario 8 (after 2 trials):
% Panel 1: Labeled training results (3 time points - 2 colors: measured & predicted) in 3D
% Panel 2: Unlabeled validation results (3 time points - 2 colors: measured & predicted) in 3D
% Panel 3: Unlabeled testing results for 4 time points (4 subplots) separately - TestLRep 3 & TestULRep 3 (2 colors: measured & predicted) 

% Call Results_Scenario8_After2Trials:
load('../../Data_Files/Data_CodesForFigures_Weighted_Object_Results/Object_Results/Step6_TRY_2.mat'); % load the weighted model for the case where we apply GA on quadratic features
CurrentObject.MakeHistPlot = 'False';
CurrentObject.GenerateFigure = 'False';

% Panel 1:
orient(figure(4),'landscape');
h4(1) = subplot(3,4,[1,2,5,6]);
DaysTrain = [1,14,46];
for j = 1:3 % Get the resuls for the optimized & final weighted modeling strategy:
    Scen8_YtrL{j,1} = exp(CurrentObject.RegressModel.LogYTrainL{j,1}); % Measured labeled training targetes
    Scen8_YprtrL{j,1} = exp(CurrentObject.RegressModel.LogYPredTrainL{j,1}); % Predicted labeled training targets
    Scen8_HistYtrL{1,j} = hist(Scen8_YtrL{j,1},exp(CurrentObject.BINS))/sum(hist(Scen8_YtrL{j,1},exp(CurrentObject.BINS))); % Histogram of measured labeled training targets
    Scen8_HistYprtrL{1,j} = hist(Scen8_YprtrL{j,1},exp(CurrentObject.BINS))/sum(hist(Scen8_YprtrL{j,1},exp(CurrentObject.BINS))); % Histogram of predicted labeled training targets
    HistPlotsYtrL_Scen8 = plot3(j*ones(1,500),exp(CurrentObject.BINS)',[Scen8_HistYtrL{1,j}]','b','LineWidth',2); % 3D plot of histogram for measured labeled training targets
    hold on
    HistPlotsYprtrL_Scen8 = plot3(j*ones(1,500),exp(CurrentObject.BINS)',[Scen8_HistYprtrL{1,j}]','r','LineWidth',2); % 3D plot of histogram for predicted labeled training targets
    DKStr_Scen8(j) = CurrentObject.RegressModel.ErrTrainL(j); % Get the KS distance
    text(j,1.41e6,0.083,sprintf(['KS_{',num2str(DaysTrain(j)),'} = %1.4f'],eval(num2str(DKStr_Scen8(j)))),'FontSize',24); % Display the KS on plots
end
xticklabels([1,2,3]);
get(gca,'XTickLabel');
set(gca,'FontSize',20);
get(gca,'YTickLabel');
set(gca,'FontSize',20);
xlabel('Time Point','FontSize',24,'Rotation',-51);
ylim([1e2,1e7]);
set(gca,'YScale','log');
ylabel('Lipid Content (AUC)','FontSize',24,'Rotation',7);
zlabel('Probability','FontSize',24);
title('A','FontSize',24);
grid on
view(70,30);
h4L = legend(h4(1),[HistPlotsYtrL_Scen8,HistPlotsYprtrL_Scen8],{'Measured','Estimated from non-label channels'},'Location','northwest','FontSize',22);
% newPosition = [0.09 0.895 0.2 0.15];
% newUnits = 'normalized';
% set(h4L,'Position', newPosition,'Units', newUnits);

% Panel 2:
h4(2) = subplot(3,4,[3,4,7,8]);
DaysValid = [0,15,37];
for j = 1:3 % Get the validation result for our optimzed & final weighted modeling strategy:
    Scen8_YvaL{j,1} = exp(CurrentObject.RegressModel.LogYTestL{j,1}); % Measured labeled validation
    Scen8_YprvaUL{j,1} = exp(CurrentObject.RegressModel.LogYPredTestUL{j,1}); % Predicted unlabeled validation
    Scen8_HistYvaL{1,j} = hist(Scen8_YvaL{j,1},exp(CurrentObject.BINS))/sum(hist(Scen8_YvaL{j,1},exp(CurrentObject.BINS))); % Histogram of Measured labeled validation
    Scen8_HistYprvaUL{1,j} = hist(Scen8_YprvaUL{j,1},exp(CurrentObject.BINS))/sum(hist(Scen8_YprvaUL{j,1},exp(CurrentObject.BINS))); % Histogram of Predicted unlabeled validation
    HistPlotsYvaL_Scen8 = plot3(j*ones(1,500),exp(CurrentObject.BINS)',[Scen8_HistYvaL{1,j}]','b','LineWidth',2); %3D plot of histograms of Measured labeled validation
    hold on
    HistPlotsYvaUL_Scen8 = plot3(j*ones(1,500),exp(CurrentObject.BINS)',[Scen8_HistYprvaUL{1,j}]','r','LineWidth',2); % 3D plot of histograms of predicted unlabeled validation
    DKSva_Scen8(j) = CurrentObject.RegressModel.ErrTestUL(j); % get the KS results
    text(j,1.41e6,0.083,sprintf(['KS_{',num2str(DaysValid(j)),'} = %1.4f'],eval(num2str(DKSva_Scen8(j)))),'FontSize',24); % Display the KS results
end
xticklabels([1,2,3]);
get(gca,'XTickLabel');
set(gca,'FontSize',20);
get(gca,'YTickLabel');
set(gca,'FontSize',20);
xlabel('Time Point','FontSize',24,'Rotation',-51);
ylim([1e2,1e7]);
set(gca,'YScale','log');
ylabel('Lipid Content (AUC)','FontSize',24,'Rotation',7);
zlabel('Probability','FontSize',24);
title('(b) Unlabeled Validation','FontSize',24);
grid on
view(70,30);

% Get the KS between days 0 & 37 for measured targets (for use in Figure 5 - panel 3):
[~,~,KS037meas] = kstest2(CurrentObject.RegressModel.LogYTestL{1,1},CurrentObject.RegressModel.LogYTestL{3,1});
% Get the KS between days 0 & 37 for predicted unlabeled targets (for use in Figure 5 - panel 3):
[~,~,KS037pred] = kstest2(CurrentObject.RegressModel.LogYPredTestUL{1,1},CurrentObject.RegressModel.LogYPredTestUL{3,1});

% Panel 3:
TilePosition = 9:12; % Subplot tile position in the figure
% TestTimePoints = [5,12,14,21]; % Selected actual time points
DaysTest = [7,16,20,34]; % Selected testing days
j = 1;
for i = [3,8,10,17] % Arbitrarily selected relative time points in the testing data to show on plots
    h4(j+2) = subplot(3,4,TilePosition(j));
    for l = 1:3 % Labeled reps
        for u = 1:4 % Unlabeled reps
            load(['../../Data_Files/Data_CodesForFigures_Weighted_Object_Results/Object_Results/Step6_TestLRep',num2str(l),'_TestULRep',num2str(u),'.mat']);
            CurrentObject.MakeHistPlot = 'False';
            CurrentObject.GenerateFigure = 'False';
            Object_Test_Scen8{l,u,j} = CurrentObject;
            Scen8_YtesL{j,1} = exp(Object_Test_Scen8{l,u,j}.RegressModel.LogYTestL{i,1});
            Scen8_YprtesUL{j,1} = exp(Object_Test_Scen8{l,u,j}.RegressModel.LogYPredTestUL{i,1});
            Scen8_HistYtesL{l,u,j} = hist(Scen8_YtesL{j,1},exp(Object_Test_Scen8{l,u,j}.BINS))/sum(hist(Scen8_YtesL{j,1},exp(Object_Test_Scen8{l,u,j}.BINS)));
            Scen8_HistYprtesUL{l,u,j} = hist(Scen8_YprtesUL{j,1},exp(Object_Test_Scen8{l,u,j}.BINS))/sum(hist(Scen8_YprtesUL{j,1},exp(Object_Test_Scen8{l,u,j}.BINS)));
            HistPlotsYtesL_Scen8{l,u} = plot(exp(Object_Test_Scen8{l,u,j}.BINS)',[Scen8_HistYtesL{l,u,j}]','b','LineWidth',2);
            hold on
            HistPlotsYprtesUL_Scen8{l,u} = plot(exp(Object_Test_Scen8{l,u,j}.BINS)',[Scen8_HistYprtesUL{l,u,j}]','r','LineWidth',2);
            KSts_Scen8(l,u,j) = Object_Test_Scen8{l,u,j}.RegressModel.ErrTestUL(i);
        end
    end
    
    KSavgTesScen8(j) = mean(mean(KSts_Scen8(:,:,j)));
    text(1.2e2,0.12,sprintf('<KS>_{pred} = %1.4f',eval(num2str(KSavgTesScen8(j)))),'FontSize',22);
    
    [~,~,KS_ideal(1,j)] = kstest2(Object_Test_Scen8{1,1,j}.RegressModel.LogYTestL{i,1},Object_Test_Scen8{2,1,j}.RegressModel.LogYTestL{i,1}); % The ideal KS for reps 1 & 2
    [~,~,KS_ideal(2,j)] = kstest2(Object_Test_Scen8{1,1,j}.RegressModel.LogYTestL{i,1},Object_Test_Scen8{3,1,j}.RegressModel.LogYTestL{i,1}); % The ideal KS for reps 1 & 3
    [~,~,KS_ideal(3,j)] = kstest2(Object_Test_Scen8{2,1,j}.RegressModel.LogYTestL{i,1},Object_Test_Scen8{3,1,j}.RegressModel.LogYTestL{i,1}); % The ideal KS for reps 2 & 3
    KSavgIdeal(j) = mean(KS_ideal(:,j)); % Average of ideal KS's
    text(1.2e2,0.1,sprintf('<KS>_{data} = %1.4f',eval(num2str(KSavgIdeal(j)))),'FontSize',22);
    
    xlim([1e2,1e7]);
    set(gca,'XScale','log');
    ylim([0,0.13]);
    get(gca,'XTickLabel');
    set(gca,'FontSize',20);
    get(gca,'YTickLabel');
    set(gca,'FontSize',20);
    xlabel('Lipid Content (AUC)','FontSize',22);
    ylabel('Probability','FontSize',22);
    title(['(',char('b'+(j)),') Unlabeled Testing Day ',num2str(DaysTest(j))],'FontSize',24);
    j = j + 1;
end
saveas(figure(4),'Paper_Figures/mainFigs/Figure_4_FINAL_weightedModel_Valid_Test.fig');
set(gcf, 'PaperPosition', [0 0 27 20]);
set(gcf, 'PaperSize', [27 20]);
print(figure(4),'Paper_Figures/mainFigs/Figure_4_FINAL_weightedModel_Valid_Test','-dpdf','-fillpage');

%% Figure 5: Average Lipid Accumulation per Day, Weights Based on Label-Free Information, KS Comparisons b/w 0 & 37, and KS Comparisons wrt Time
% This figure contains 4 panels:
% Panel 1: Average lipid accumulation per day (mean +/- std.dev. vs time with shaded erorr bars)
% Panel 2: The change in the weights alhpha 1, alpha 2, and alpha 3 with respect to testing time points
% Panel 3: The KS distance comparison between measured values of the lipid contents for days 0 vs 37, and the KS distance comparison between predicted values of the lipid contents for days 0 vs 37 (VALIDATION TIME POINTS)
% Panel 4: Comparison of the KS distances with respect to time

% Panel 1:
% Pull data from Results_Scenario8_After2Trials for training and validation results:
clear CurrentObject;
load('../../Data_Files/Data_CodesForFigures_Weighted_Object_Results/Object_Results/Step6_TRY_2.mat');
CurrentObject.MakeHistPlot = 'False';
CurrentObject.GenerateFigure = 'False';
Object_Scen8 = CurrentObject;
% Pull data from Results_Scenario8_After2Trials for testLrep 3 and testULrep 3 (arbitrarily):
clear CurrentObject;
load('../../Data_Files/Data_CodesForFigures_Weighted_Object_Results/Object_Results/Step6_TestLRep3_TestULRep3.mat');
CurrentObject.MakeHistPlot = 'False';
CurrentObject.GenerateFigure = 'False';
Object_Test_Scen8 = CurrentObject;
YL_Scen8 = cell(23,1);
YPredUL_Scen8 = cell(23,1);
j = 1;
for i = [2,10,23] % Training time points target values (measured & predicted)
    YL_Scen8{i,1} = exp(Object_Scen8.RegressModel.LogYTrainL{j,1}); % Measured labeled
    YPredUL_Scen8{i,1} = exp(Object_Scen8.RegressModel.LogYPredTrainUL{j,1}); % Predicted unlabeld
    j = j + 1;
end
j = 1;
for i = [1,11,22] % Validaiton time points target values (measured & predicted)
    YL_Scen8{i,1} = exp(Object_Scen8.RegressModel.LogYTestL{j,1});
    YPredUL_Scen8{i,1} = exp(Object_Scen8.RegressModel.LogYPredTestUL{j,1});
    j = j + 1;
end
j = 1;
for i = [3:9,12:21] % Testing time points target values (measured & predicted)
    YL_Scen8{i,1} = exp(Object_Test_Scen8.RegressModel.LogYTestL{j,1});
    YPredUL_Scen8{i,1} = exp(Object_Test_Scen8.RegressModel.LogYPredTestUL{j,1});
    j = j + 1;
end
% Compute means and standard deviations for the measured (mu & stdv) and predicted (muhat & stdvhat) targets:
mu = []; stdv = []; muhat = []; stdvhat = [];
for i = 1:23
    mu = [mu;mean(YL_Scen8{i,1})];
    stdv = [stdv;std(YL_Scen8{i,1})];
    muhat = [muhat;mean(YPredUL_Scen8{i,1})];
    stdvhat = [stdvhat;std(YPredUL_Scen8{i,1})];
end
% Plot the average lipid per day for measured and unlabeled predicted result:
orient(figure(5),'landscape');
% h5(1) = subplot(2,5,1:3);
h5(1) = subplot(2,4,1:2);
Days = [0,1,5,6,7,9,10,11,12,14,15,16,19,20,21,22,23,26,27,30,34,37,46]; % The days that samples are taken
TrainSet = [2,10,23]; %Days 1, 14, and 46
ValidSet = [1,11,22]; %Days 0, 15, and 37
AvgLip_Measured = shadedErrorBar(Days,mu,stdv,'lineprops',{'b-o','markerfacecolor','b','MarkerSize',10},'transparent',1); %Measured lipid content
hold on;
AvgLip_Predicted = shadedErrorBar(Days,muhat,stdvhat,'lineprops',{'r-o','markerfacecolor','r','MarkerSize',10},'transparent',1); %Predicted lipid content
TrainPlot = plot(Days(TrainSet),mu(TrainSet),'s','MarkerFaceColor',[0.22,0.22,0.22],'MarkerEdgeColor','none','MarkerSize',20);
ValidPlot = plot(Days(ValidSet),muhat(ValidSet),'^','MarkerFaceColor',[0.4660,0.6740,0.1880],'MarkerEdgeColor','none','MarkerSize',20);
xticks(0:5:50);
xlim([-1,50]);
ylim([5e3,1.5e6]);
set(gca,'YScale','log');
xlabel('Days After Nitrogen Starvation','FontSize',26);
ylabel('Average Lipid Content (AUC)','FontSize',26);
get(gca,'XTickLabel');
set(gca,'FontSize',24);
get(gca,'YTickLabel');
set(gca,'FontSize',24);
legend([AvgLip_Measured.mainLine AvgLip_Predicted.mainLine TrainPlot ValidPlot],...
    {'\mu Measured','\mu Label-free','Training','Validation',},'Location','northwest','FontSize',22);
title('A','FontSize',26);

% Panel 2:
% h5(2) = subplot(2,5,4:5);
h5(2) = subplot(2,4,3:4);
TestSet = [3:9,12:21]; % Testing time points
[alpha_1,alpha_2,alpha_3] = deal(cell(0));
for testLrep = 1:3
    for testUrep = 1:4
        clear CurrentObject;
        load(['../../Data_Files/Data_CodesForFigures_Weighted_Object_Results/Object_Results/Step6_TestLRep',num2str(testLrep),'_TestULRep',num2str(testUrep),'.mat']);
        CurrentObject.MakeHistPlot = 'False';
        CurrentObject.GenerateFigure = 'False';
        RegMod = CurrentObject.RegressModel;
        alpha_1{testLrep,testUrep} = plot(Days(TestSet),RegMod.AlphaTsUL(:,1),'Color',[0,0.4470,0.7410]); % Blue color for alpha_1
        hold on;
        alpha_2{testLrep,testUrep} = plot(Days(TestSet),RegMod.AlphaTsUL(:,2),'Color',[0.6350,0.0780,0.1840]); % Red color for alpha_2
        alpha_3{testLrep,testUrep} = plot(Days(TestSet),RegMod.AlphaTsUL(:,3),'Color',[0.4660,0.6740,0.1880]); % Green color for alpha_3
    end
end
xlim([0,50]);
legend({'\alpha_1','\alpha_2','\alpha_3'},'Location','southeast','FontSize',22);
title(legend,'Model weights');
% set(gca,'XTick',Days(TestSet));
xticks(0:5:50);
get(gca,'XTickLabel');
set(gca,'FontSize',24);
get(gca,'YTickLabel');
set(gca,'FontSize',24);
xlabel('Days After Nitrogen Starvation','FontSize',26);
ylabel('Model Weights','FontSize',26);
title('B','FontSize',26);

% Panel 3:
% h5(3) = subplot(2,5,6:7);
h5(3) = subplot(2,4,5:6);
KScompar_histPlotYvaL_Day0 = plot(exp(CurrentObject.BINS)',Scen8_HistYvaL{1,1}','b','LineWidth',2); % Plot the measured histograms for day 0
hold on;
KScompar_histPlotYvaL_Day37 = plot(exp(CurrentObject.BINS)',Scen8_HistYvaL{1,3}','b','LineWidth',2); % Plot the measured histograms for day 37
KScompar_histPlotYvaprUL_Day0 = plot(exp(CurrentObject.BINS)',Scen8_HistYprvaUL{1,1}','r','LineWidth',2); % Plot the predicted histograms for day 0
KScompar_histPlotYvaprUL_Day37 = plot(exp(CurrentObject.BINS)',Scen8_HistYprvaUL{1,3}','r','LineWidth',2); % Plot the predicted histograms for day 37
text(9,0.03,'Day 0','FontSize',24);
text(14.5,0.03,'Day 37','FontSize',24);
% text(5.1,0.095,sprintf('KS_{m_{0, 37}} = %1.4f',eval(num2str(KS037meas))),'FontSize',24);
% text(5.1,0.085,sprintf('KS_{p_{0, 37}} = %1.4f',eval(num2str(KS037pred))),'FontSize',24);
ylim([0,0.1]);
xlim([1.5e3,1e7]);
set(gca,'XScale','log');
get(gca,'XTickLabel');
set(gca,'FontSize',24);
get(gca,'YTickLabel');
set(gca,'FontSize',24);
xlabel('Lipid Content (AUC)','FontSize',26);
ylabel('Probability','FontSize',26);
title('C','FontSize',26);

% Panel 4:
clear CurrentObject;
% Load the weighted model
load('../../Data_Files/Data_CodesForFigures_Weighted_Object_Results/Object_Results/Step6_TRY_2.mat');
CurrentObject.MakeHistPlot = 'False';
CurrentObject.GenerateFigure = 'False';
Model = CurrentObject.RegressModel;
% Load testing results
testResults = cell(0);
for testLrep = 1:3
    for testUrep = 1:4
        load(['../../Data_Files/Data_CodesForFigures_Weighted_Object_Results/Object_Results/Step6_TestLRep',num2str(testLrep),'_TestULRep',num2str(testUrep),'.mat']);
        CurrentObject.MakeHistPlot = 'False';
        CurrentObject.GenerateFigure = 'False';
        testResults{testLrep,testUrep} = CurrentObject.RegressModel;
    end
end
trainTimes = [2,10,23];
validTimes = [1,11,22];
testTimes = [3:9,12:21];

% KS comparisons between day 0 and every other day following nitrogen starvation:
[KSm,KSp] = deal(zeros(23,1));
% KS between measured dates:
i = 1;
for k = trainTimes
    [~,~,KSm(k)] = kstest2(Model.LogYTestL{1,1},Model.LogYTrainL{i,1});
    i = i + 1;
end
i = 1;
for k = validTimes
    [~,~,KSm(k)] = kstest2(Model.LogYTestL{1,1},Model.LogYTestL{i,1});
    i = i + 1;
end
i = 1;
for k = testTimes
    [~,~,KSm(k)] = kstest2(Model.LogYTestL{1,1},testResults{1,1}.LogYTestL{i,1});
    i = i + 1;
end

% KS between predicted dates:
i = 1;
for k = trainTimes
    [~,~,KSp(k)] = kstest2(Model.LogYPredTestUL{1,1},Model.LogYPredTrainUL{i,1});
    i = i + 1;
end
i = 1;
for k = validTimes
    [~,~,KSp(k)] = kstest2(Model.LogYPredTestUL{1,1},Model.LogYPredTestUL{i,1});
    i = i + 1;
end
i = 1;
for k = testTimes
    [~,~,KSp(k)] = kstest2(Model.LogYPredTestUL{1,1},testResults{1,1}.LogYPredTestUL{i,1});
    i = i + 1;
end

% KS comparisons between the last day (day 46) and each day following nitrogen starvation (reverse):
[KSm_rev,KSp_rev] = deal(zeros(23,1));
% KS between measured dates:
i = 1;
for k = trainTimes
    [~,~,KSm_rev(k)] = kstest2(Model.LogYTrainL{3,1},Model.LogYTrainL{i,1});
    i = i + 1;
end
i = 1;
for k = validTimes
    [~,~,KSm_rev(k)] = kstest2(Model.LogYTrainL{3,1},Model.LogYTestL{i,1});
    i = i + 1;
end
i = 1;
for k = testTimes
    [~,~,KSm_rev(k)] = kstest2(Model.LogYTrainL{3,1},testResults{1,1}.LogYTestL{i,1});
    i = i + 1;
end

% KS between predicted dates:
i = 1;
for k = trainTimes
    [~,~,KSp_rev(k)] = kstest2(Model.LogYPredTrainUL{3,1},Model.LogYPredTrainUL{i,1});
    i = i + 1;
end
i = 1;
for k = validTimes
    [~,~,KSp_rev(k)] = kstest2(Model.LogYPredTrainUL{3,1},Model.LogYPredTestUL{i,1});
    i = i + 1;
end
i = 1;
for k = testTimes
    [~,~,KSp_rev(k)] = kstest2(Model.LogYPredTrainUL{3,1},testResults{1,1}.LogYPredTestUL{i,1});
    i = i + 1;
end

% Plot the KS and KS_rev:
% h5(4) = subplot(2,5,8:10);
h5(4) = subplot(2,4,7:8);
days = [0,1,5,6,7,9,10,11,12,14,15,16,19,20,21,22,23,26,27,30,34,37,46];
plot(days,KSm,'LineWidth',3);
hold on;
plot(days,KSp,'LineWidth',3);
plot(days,KSm_rev,'LineWidth',3);
plot(days,KSp_rev,'LineWidth',3);
xticks(0:5:50);
get(gca,'XTickLabel');
set(gca,'FontSize',24);
get(gca,'YTickLabel');
set(gca,'FontSize',24);
xlabel('Days After Nitrogen Starvation','FontSize',26);
ylabel('KS Distance','FontSize',26);
legend({'Measured (day 0)','Label free (day 0)',...
    'Measured (day 46)','Label free (day 46)'},'Location','east','FontSize',22);
title('D','FontSize',26);
saveas(figure(5),'Paper_Figures/mainFigs/Figure_5_AvgLip_Alphas_KScomparisons.fig');
set(gcf, 'PaperPosition', [0 0 25 25]);
set(gcf, 'PaperSize', [25 25]);
print(figure(5),'Paper_Figures/mainFigs/Figure_5_AvgLip_Alphas_KScomparisons','-dpdf','-fillpage');

%% Figure 6: Sorting based on our Approach
% This figure has 2 panels:
% Panel 1: Sorting of the labeled data for the original data (random samples are taken on a 50%/50% basis)
% Panel 2: Sorting of the unlabeled data for the new FCM pico data (random samples are taken on a 50%/50% basis)

% Panels 1 & 2:
% Here we take random samples from labeled and unlabeled cells in the days corresponding to the lowest and the highest signal intensities (days 0 & 37 -- validation time points).
% We then take 2500 cells out of each time point.
% The random samples are taken on a 50%/50% basis.
% We test our developed method on the above data set by using previously computed M1, M2, M3, and Q (no new learning from these data).
% We then generate historgrams for the above predictions of the target signal.
clear CurrentObject;
% Load the data set and extract 5000 cells from low and high signals:
load('../../Data_Files/NewPico_Data_NoZero.mat');
TargetSet = [3,9]; %Enter the target parmeters' index (here 3 and 9 are targets).
TargetInd = TargetSet(1,1); %Target the user looks for. Alternative: Target_Ind = Targets_Set(1,2);
TimePoint = 23; %Enter number of time points (days after nitorgen starvation) to be selected.
NCellsUL = 1250; %Enter number of cells to be used from each replication at each timepoint (unlabeled cells).
NCellsL = 1667; %Enter number of cells to be used from each replication at each timepoint (labeled cells).
TrainSet = [2,10,23]; %[SortInd(1) SortInd(end)]; %Widest range in Bodipy_Signals.
ValidSet = [1,11,22]; %[SortInd(2) SortInd(end-1)]; %Second widest range in Bodipy_Signals.
FeatureInds = 1:size(DATA.FeatureNames(~strcmp(DATA.FeatureNames,'Time')),2); %Define Feature_Inds
FeatureInds(TargetSet) = []; %Remove indicies for tragets.
FeatureInds_S = sort([FeatureInds,TargetSet]); %Define Feature_Inds for S test statistic
[~,~,~,~,XValidL,XValidUL,XValidUL_S,~] = DataCollector_Concat_TrainValid(NCellsL,NCellsUL,FeatureInds,FeatureInds_S,TargetInd,TrainSet,ValidSet,DATA); % Collect the 5000 cells from unlabeled data
% Add quad features:
[XValidLQ,XValidULQ,XValidULQ_S] = deal(cell(0));
for i=1:size(XValidUL,1)
    XValidLQ{i,1} = Get_Quad_Features(XValidL{i,1}); % Get quadratic feqatures with the function Get_Quad_Features
    XValidULQ{i,1} = Get_Quad_Features(XValidUL{i,1});
    XValidULQ_S{i,1} = Get_Quad_Features(XValidUL_S{i,1});
end
% Convert to log:
[LogXValidLQ,LogXValidULQ,LogXValidULQ_S] = deal(cell(0));
for i = 1:size(XValidUL,1)
    LogXValidLQ{i,1} = log(XValidLQ{i,1});
    LogXValidULQ{i,1} = log(XValidULQ{i,1});
    LogXValidULQ_S{i,1} = log(XValidULQ_S{i,1});
end
% Load the final model and extract the required information:
load('../../Data_Files/Data_CodesForFigures_Weighted_Object_Results/Object_Results/Step6_TRY_2.mat'); % Load the final weighted model to get M1, M2, M3, and Q
CurrentObject.MakeHistPlot = 'False'; % Prevent the current object from make histogram plots when run
CurrentObject.GenerateFigure = 'False'; % Prevent the current object from generating figures
M = CurrentObject.M; % Get the models M1, M2, and M3
Q = CurrentObject.Q; % Get the weight quotient
IvecBest1 = CurrentObject.I(:,1); % Seletced features for M1
IvecBest2 = CurrentObject.I(:,2); % Seletced features for M2
IvecBest3 = CurrentObject.I(:,3); % Seletced features for M3
IforS = CurrentObject.IforS; % Selected columns of S (rows of Q)
BINS = linspace(6,57,1000); % Specify bins

% 50% from low lipids and 50% from high lipids:
rng(1); % Fix the random number generator
% The 50% from LOW signal intensities:
R1 = ceil(5000*rand(2500,1)); % Generate 2500 random indices from the range 1 to 3000 (the low 50%)
Xl1 = LogXValidLQ{1,1}(R1,:); % Get the 50% data by taking random samples from the LOW signals (time point 0) -- labeled
Xu1 = LogXValidULQ{1,1}(R1,:); % Get the 50% data by taking random samples from the LOW signals (time point 0) -- unlabeled
Xu1_s = LogXValidULQ_S{1,1}(R1,:); % Get the data for statistics matrix by the SAME random indicies at the LOW signal (time point 0) -- unlabeled
% The 50% from HIGH signal intensities:
R2 = ceil(5000*rand(2500,1)); % Generate 2500 random indices from the range 1 to 3000 (the high 50%)
Xl2 = LogXValidLQ{3,1}(R2,:); % Get the 50% data by taking random samples from the LOW signals (time point 0) -- labeled
Xu2 = LogXValidULQ{3,1}(R2,:); % Get the 50% data by taking random samples from the LOW signals (time point 0) -- unlabeled
Xu2_s = LogXValidULQ_S{3,1}(R2,:); % Get the data for statistics matrix by the SAME random indicies at the LOW signal (time point 0) -- unlabeled
% Stack the two populations on top of each other:
Xl = [Xl1;Xl2];
Xu = [Xu1;Xu2];
Xu_s = [Xu1_s;Xu2_s];
% Compute the weights based on unlabeled data:
S = [mean(Xu_s),sqrt(var(Xu_s))]; % Get the test statistic matrix
alpha = S(:,IforS==1)*Q(IforS==1,:); % Compute the weights as alpha = S\Q
Mw = alpha(1,1)*M{1,1} + alpha(1,2)*M{2,1} + alpha(1,3)*M{3,1}; % Compute the weighted model
% Predict the 90% and 10% separately:
Ypredl1 = exp(Xl1*Mw); % Predict the labeled targets
Ypredl2 = exp(Xl2*Mw); % Predict the labeled targets
YpredlAll = exp(Xl*Mw); % Predict the total MIXED population of labeled targets
Ypredu1 = exp(Xu1*Mw); % Predict the unlabeled targets
Ypredu2 = exp(Xu2*Mw); % Predict the unlabeled targets
YpreduAll = exp(Xu*Mw); % Predict the total MIXED population of unlabeled targets

% Plot the predicted 50% and 50% in a figure
% orient(figure(6),'landscape');
figure(6);
% Panel 1: sorting labeled cells
h6(1) = subplot(2,1,1);
hist_50_low_l = hist(Ypredl1,exp(BINS))/sum(hist(YpredlAll,exp(BINS))); % Histogram for the 50% from low signal intensities
sortHistPlot_50_low_l = plot(exp(BINS),hist_50_low_l,'Color',[0.4660,0.6740,0.1880],'LineWidth',1); %Plot histograms for low singal
hold on;
hist_50_high_l = hist(Ypredl2,exp(BINS))/sum(hist(YpredlAll,exp(BINS))); % Histogram for the 50% from the high signal intensities
sortHistPlot_50_high_l = plot(exp(BINS),hist_50_high_l,'Color',[0.4940,0.1840,0.5560],'LineWidth',1); %Plot histograms for high singal
% Plot the mixed population:
sortHistAll_l = hist(YpredlAll,exp(BINS))/sum(hist(YpredlAll,exp(BINS)));
sortHistAllPlot_l = plot(exp(BINS),sortHistAll_l,'k','LineWidth',0.5); %Plot histograms for the mixed population
xlim([3e3 3.5e6]);
set(gca,'XScale','log');
ylim([0,0.025]);
% Find approximate location where the two samples can be separated:
intersect_indexL = find(abs(hist_50_low_l-hist_50_high_l)<0.0005 & hist_50_low_l>0.006 & hist_50_high_l>0.006); % Find on which index the two distributions intersect
x_IntersectL = exp(BINS(intersect_indexL)); % Find the intersect location in bins (lipid content at intersection: x axis)
YY = get(gca,'YLim');
plot(x_IntersectL*ones(1,2),YY,'k--');
text(x_IntersectL,-0.0006,'x','FontWeight','bold','FontSize',27);
% Find the percent of correclty sorted cells:
correctLow_l = sum(hist_50_low_l(1:intersect_indexL))/sum(hist_50_low_l)*100; % Low lipid cells
correctHigh_l = sum(hist_50_high_l(intersect_indexL:end))/sum(hist_50_high_l)*100; % High lipid cells
% text(5.2,0.022,['Correctly sorted low-lipid cells = ',num2str(correctLow_l),'%']);
% text(5.2,0.017,['Correctly sorted high-lipid cells = ',num2str(correctHigh_l),'%']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BELOW METHOD WOULD BE CORRECT IF WE ARE ALLOWED TO ASSUME THAT THE  PRIOR PROBABILITIES ARE EQUAL:
% [~,~,KS_l] = kstest2(Ypredl1,Ypredl2); % Compute the KS distance between low and high
% pErr_l = (1 - KS_l)/2; % Probability of error in separating the two distributions (according to https://www.nature.com/articles/s41598-017-16166-y)
% correctSort = (1 - 2*pErr_l) * 100; % Percent of correctly sorted cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('Lipid Content (AUC)','FontSize',39);
ylabel('Probability','FontSize',39);
get(gca,'XTickLabel');
set(gca,'FontSize',36);
get(gca,'YTickLabel');
set(gca,'FontSize',36);
legend({'Low lipids (50%) -- Day 0','High lipids (50%) -- Day 37','Mixed population'},'Location','northeast','FontSize',32);
title('A');

% Panel 2: sorting unlabeled cells
h6(2) = subplot(2,1,2);
hist_50_low_u = hist(Ypredu1,exp(BINS))/sum(hist(YpreduAll,exp(BINS))); % Histogram for the 50% from low signal intensities
sortHistPlot_50_low_u = plot(exp(BINS),hist_50_low_u,'Color',[0.4660,0.6740,0.1880],'LineWidth',1); %Plot histograms for low signals
hold on;
hist_50_high_u = hist(Ypredu2,exp(BINS))/sum(hist(YpreduAll,exp(BINS))); % Histogram for the 50% from the high signal intensities
sortHistPlot_50_high_u = plot(exp(BINS),hist_50_high_u,'Color',[0.4940,0.1840,0.5560],'LineWidth',1); %Plot histograms for high signals
% Plot the mixed population:
sortHistAll_u = hist(YpreduAll,exp(BINS))/sum(hist(YpreduAll,exp(BINS)));
sortHistAllPlot_u = plot(exp(BINS),sortHistAll_u,'k','LineWidth',0.5); %Plot histograms for the mixed population
xlim([3e3 3.5e6]);
set(gca,'XScale','log');
ylim([0,0.025]);
% Find approximate location where the two samples can be separated:
intersect_indexU = find(abs(hist_50_low_u-hist_50_high_u)<0.0025 & hist_50_low_u>0.005 & hist_50_high_u>0.005); % Find on which index the two distributions intersect
x_IntersectU = exp(BINS(intersect_indexU)); % Find the intersect location in bins (lipid content at intersection: x axis)
YY = get(gca,'YLim');
plot(x_IntersectU*ones(1,2),YY,'k--');
text(x_IntersectU,-0.0006,'x','FontWeight','bold','FontSize',27);
% Find the percent of correclty sorted cells:
correctLow_u = sum(hist_50_low_u(1:intersect_indexU))/sum(hist_50_low_u)*100; % Low lipid cells
correctHigh_u = sum(hist_50_high_u(intersect_indexU:end))/sum(hist_50_high_u)*100; % High lipid cells
% text(5.2,0.022,['Correctly sorted low-lipid cells = ',num2str(correctLow_u),'%']);
% text(5.2,0.017,['Correctly sorted high-lipid cells = ',num2str(correctHigh_u),'%']);
xlabel('Lipid Content (AUC)','FontSize',39);
ylabel('Probability','FontSize',39);
get(gca,'XTickLabel');
set(gca,'FontSize',36);
get(gca,'YTickLabel');
set(gca,'FontSize',36);
title('B');
saveas(figure(6),'Paper_Figures/mainFigs/Figure_6_Sorting.fig');
set(gcf, 'PaperPosition', [0 0 15 18]);
set(gcf, 'PaperSize', [15 18]);
print(figure(6),'Paper_Figures/mainFigs/Figure_6_Sorting','-dpdf','-fillpage');

%% Figure 7: Generality of our Approach
% Panel 1: Original data: unlabeled features and labeled lipid content (median values)
% Panel 2: New data: unlabeled features and labeled lipid content (median values)
% panel 3: Average lipid content per day for the new FCM pico data (measured vs predicted)

% Collect the requried data for panels 1 & 2:
selectedFeats = {'SSC-A','FL3-A','FL4-A','FL1-H'}; % Selected features by the GA in our method for the primary and secondary regression analyses
DATAorig = load('../NewPico_Data_NoZero.mat'); % Original data set
DATAnew = load('../NewAnalyses/FCM_Pico_Replica1_Data_NoZero.mat'); % New data set

for i=1:length(selectedFeats)
    IND_FET(i) = find(ismember(DATAorig.DATA.FeatureNames,selectedFeats{i}));
end
for i=1:23
    data_origl = DATAorig.DATA.Stained(i,3).DATA(:,IND_FET); % Get selected features and targets data (original data set; labeled cells; first day)
    data_origu = DATAorig.DATA.Unstained(i,3).DATA(:,IND_FET); % Get selected features and targets data (original data set; unlabeled cells; first day)
    mu_origl(i,:) = median(data_origl);
    mu_origu(i,:) = median(data_origu);
    try
        data_newl = DATAnew.DATA.Stained(i,3).DATA(:,IND_FET); % Get selected features and targets data (original data set; labeled cells; first day)
        data_newu = DATAnew.DATA.Unstained(i,3).DATA(:,IND_FET); % Get selected features and targets data (original data set; unlabeled cells; first day)
        mu_newl(i,:) = median(data_newl);
        mu_newu(i,:) = median(data_newu);
    catch
    end
end


% Panel 1:
orient(figure(7),'landscape');
h7(1) = subplot(2,2,1);
plot(DATAorig.DATA.Time,mu_origu(:,1:end-1),'linewidth',3);
hold on
plot(DATAorig.DATA.Time,mu_origl(:,end),'--','linewidth',3);
set(gca,'YScale','log');
ylim([3E3,1E6]);
xlim([0,50]);
get(gca,'XTickLabel');
set(gca,'FontSize',36);
get(gca,'YTickLabel');
set(gca,'FontSize',36);
xlabel('Time (days)','FontSize',39);
ylabel('Median Feature (AU)','FontSize',39);
title('A');

% Panel 2:
days_FCMpico = [0,1,3,5,8,10,12,15,17,19,22,24]; % Actual days after nitrogen starvation, assuming that the first day of measurement is the day 0 after nitrogen starvation:
                                                 % Correspnding to the dates {'033019','040119','040319','040519','040819','041019','041219','041519','041719','041919','042219','042419'}
h7(2) = subplot(2,2,2);
plot(days_FCMpico,mu_newu(:,1:end-1),'linewidth',3);
hold on
plot(days_FCMpico,mu_newl(:,end),'--','linewidth',3);
set(gca,'YScale','log');
ylim([3E3,1E6]);
xlim([0,25]);
get(gca,'XTickLabel');
set(gca,'FontSize',36);
get(gca,'YTickLabel');
set(gca,'FontSize',36);
xlabel('Time (days)','FontSize',39);
legend({selectedFeats{1:end-1},'Lipid Content'},'Location','NorthWest','Fontsize',32);
title('B');

% Panel 3: Average lipid content (predicted vs measured) for the new data set
clear CurrentObject;
TimePoint_FCMpico = 12; % Number of measurement time points (12 days of measurements)
load('../../Data_Files/Data_NewAnalyses/Step6_TestLRep3_TestULRep3.mat'); % Load the results object for testing our method on new FCM pico data
CurrentObject.MakeHistPlot = 'False';
CurrentObject.GenerateFigure = 'False';
Object_Test_Scen8_FCMpico = CurrentObject;
LogYL_Scen8_FCMpico = cell(TimePoint_FCMpico,1);
LogYPredUL_Scen8_FCMpico = cell(TimePoint_FCMpico,1);
j = 1;
for i = 1:TimePoint_FCMpico % Testing time points (all of the time points for this new FCM_pico data set are considered as testing)
    LogYL_Scen8_FCMpico{i,1} = Object_Test_Scen8_FCMpico.RegressModel.LogYTestL{j,1};
    LogYPredUL_Scen8_FCMpico{i,1} = Object_Test_Scen8_FCMpico.RegressModel.LogYPredTestUL{j,1};
    LogYL_Scen8_FCMpico{i,1} = exp(LogYL_Scen8_FCMpico{i,1} );
    LogYPredUL_Scen8_FCMpico{i,1} = exp(LogYPredUL_Scen8_FCMpico{i,1});
    j = j + 1;
end
% Compute means and standard deviations for measured and predicted results:
mu_FCMpico = []; stdv_FCMpico = []; muhat_FCMpico = []; stdvhat_FCMpico = [];
for i = 1:TimePoint_FCMpico
    mu_FCMpico = [mu_FCMpico;mean(LogYL_Scen8_FCMpico{i,1})];
    %     stdv_FCMpico = [stdv_FCMpico;std(LogYL_Scen8_FCMpico{i,1})];
    muhat_FCMpico = [muhat_FCMpico;mean(LogYPredUL_Scen8_FCMpico{i,1})];
    %     stdvhat_FCMpico = [stdvhat_FCMpico;std(LogYPredUL_Scen8_FCMpico{i,1})];
end
% Scale the predicted means to have the same mean at the initial time point:
% muhatScaled = muhat.*(mu(1)/muhat(1)); % Original data (predicted)
% muScaled_FCMpico = mu_FCMpico.*(mu(1)/mu_FCMpico(1));% New FCM pico data (measured)
muhatScaled_FCMpico = muhat_FCMpico.*(mu_FCMpico(1)/muhat_FCMpico(1)); % New FCM pico data (predicted)

h7(3) = subplot(2,2,[3,4]);
days_FCMpico = [0,1,3,5,8,10,12,15,17,19,22,24]; % Actual days after nitrogen starvation, assuming that the first day of measurement is the day 0 after nitrogen starvation:
                                                 % Correspnding to the dates {'033019','040119','040319','040519','040819','041019','041219','041519','041719','041919','042219','042419'}
% k = 1:TimePoint_FCMpico;
% set(gca,'XTick',k);
% set(gca,'XTickLabel',k);
% avgLipids_original_new = plot(Days,[mu,muhatScaled]','-*'); % Plot the mean values of the measured and predicted lipid content for the original data
% hold on;
avgLipids_new_l = plot(days_FCMpico,mu_FCMpico','-*','LineWidth',3,'MarkerSize',20); % Plot the mean values of the measured lipid content for the new FCM pico data
hold on
avgLipids_new_u = plot(days_FCMpico,muhatScaled_FCMpico','-o','LineWidth',3,'MarkerSize',20); % Plot the mean values of the predicted lipid content for the new FCM pico data
xticks(0:5:50);
xlim([-1,26]);
set(gca,'YScale','log');
get(gca,'XTickLabel');
set(gca,'FontSize',36);
get(gca,'YTickLabel');
set(gca,'FontSize',36);
xlabel('Days After Nitrogen Starvation','FontSize',39);
ylabel('Average Lipid Content (AUC)','FontSize',39);
title('C');
legend({'Measured','Label free'},'Location','southeast','FontSize',32);

saveas(figure(7),'Paper_Figures/mainFigs/Figure_7_Generality.fig');
set(gcf, 'PaperPosition', [0 0 16 22]);
set(gcf, 'PaperSize', [16 22]);
print(figure(7),'Paper_Figures/mainFigs/Figure_7_Generality','-dpdf','-fillpage');
