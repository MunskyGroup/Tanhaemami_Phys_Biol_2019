%% Paper_Extract_features_models_columns_NEW
%Author: Mohammad Tanhaemami
%Organization: Brian Munsky Group, Department of Chemical and Biological Engineering, Colorado State University.
% This script extract the features selected by the GA applied on linear and quadratic features (scen 2 & scen 8)
% This script then uses the extracted features to get the model that was used in each scenario
% This script extracts the columns selected by the GA applied on the secondary statistics (scen 8)
% This script then uses the extracted columns and the test statistics matrix (S) or the weight quotient (Q) to get the reduced S (or Q)

clear; clc;

%% Call the Objects for the regular approach and the weighted-model strategy (GA linear and GA quadratic)
% Call Results_regular_approach:
i = 1;
scen2 = cell(0);
for s = [4,6]
    scen2{i}.step = s;
    load(['../../Data_Files/Data_CodesForFigures_Regular_Object_Results/Object_Results/Step',num2str(s),'.mat']);
    CurrentObject.MakeHistPlot = 'False';
    CurrentObject.GenerateFigure = 'False';
    scen2{i}.features = CurrentObject.I;
    scen2{i}.maps = CurrentObject.RegressModel.MAll;
    if size(scen2{i}.features,2) > 1
        scen2{i}.features = scen2{i}.features(:);
    end
    scen2{i}.model = zeros(size(scen2{i}.features));
    scen2{i}.model(scen2{i}.features==1) = scen2{i}.model(scen2{i}.features==1) + scen2{i}.maps;
    i = i + 1;
end
% Call Results_weighted_model:
i = 1;
scen8 = cell(0);
for s = [4,6]
    scen8{i}.step = s;
    load(['../../Data_Files/Data_CodesForFigures_Weighted_Object_Results/Object_Results/Step',num2str(s),'_TRY_2.mat']);
    CurrentObject.MakeHistPlot = 'False';
    CurrentObject.GenerateFigure = 'False';
    scen8{i}.fetures = CurrentObject.I;
    scen8{i}.model = CurrentObject.M;
    scen8{i}.ifors = CurrentObject.IforS;
    scen8{i}.quotient = zeros(size(CurrentObject.Q));
    scen8{i}.quotient(CurrentObject.IforS==1,:) = scen8{i}.quotient(CurrentObject.IforS==1,:) + CurrentObject.Q(CurrentObject.IforS==1,:);
    i = i + 1;
end

%% Save the Above Information
% scen2: sheet 1
xlswrite('Paper_FeatModCol',scen2{1}.model,'Scen2','D2'); % Write scen2 GA model starting at A2
xlswrite('Paper_FeatModCol',scen2{2}.model,'Scen2','F2'); % Write scen2 GAQ model starting at E2
% scen8: sheet 2
xlswrite('Paper_FeatModCol',scen8{1}.model{1},'Scen8','D2'); % Write scen8 GA model 1 starting at A2
xlswrite('Paper_FeatModCol',scen8{1}.model{2},'Scen8','E2'); % Write scen8 GA model 2 starting at B2
xlswrite('Paper_FeatModCol',scen8{1}.model{3},'Scen8','F2'); % Write scen8 GA model 3 starting at C2
xlswrite('Paper_FeatModCol',scen8{2}.model{1},'Scen8','H2'); % Write scen8 GAQ model 1 starting at I2
xlswrite('Paper_FeatModCol',scen8{2}.model{2},'Scen8','I2'); % Write scen8 GAQ model 2 starting at J2
xlswrite('Paper_FeatModCol',scen8{2}.model{3},'Scen8','J2'); % Write scen8 GAQ model 3 starting at K2
xlswrite('Paper_FeatModCol',scen8{1}.quotient,'Scen8','L2'); % Write scen8 GA quoitient starting at E2
xlswrite('Paper_FeatModCol',scen8{2}.quotient,'Scen8','Q2'); % Write scen8 GAQ quoitient starting at M2


