%% Load NewPico Data and Initial Values
% Below codes are from Blasi et al., but modified to be applied on our data to predict the lipid content
clear; close all; clc;
load('../../Data_Files/NewPico_Data.mat') %Load the file Pico_Data.mat.Target_Ind = 3; %Choose the target parmeter index(3 and 10 are targets).
N_Cells = 3000;  % Number of cells to use from each data set.
Feature_Inds = [1 2 4 5 6 7 8 10 11 12 13]; %Features to use(3 and 9 are targets).
Target_Ind = 3; %Choose the target parmeter index(3 and 10 are targets).
warning('off')

%% Collect Ground Truth for the Lipid Content
train_data = []; % Labeled training data
train_data_U = []; % Unlabeled training data
train_ground_truth = []; % Training ground truth
for k=[2 10 23] % Train time points
    train_data = [train_data;DATA.Stained(k,2).DATA(1:N_Cells,Feature_Inds)]; % (we chose replica 2 arbitrarily)
    train_ground_truth = [train_ground_truth;DATA.Stained(k,2).DATA(1:N_Cells,Target_Ind)]; % Measured labeled targets (we chose replica 2 arbitrarily).
    train_data_U = [train_data_U;DATA.Unstained(k,2).DATA(1:N_Cells,Feature_Inds)]; % (we chose replica 2 arbitrarily)
end

test_data = []; % Testing data
test_data_U = []; % Unlabeled testing data
test_ground_truth = []; % Testing grouns truth
for k=[3:9,12:21] % Test time points
    test_data = [test_data;DATA.Stained(k,2).DATA(1:N_Cells,Feature_Inds)]; % (we chose replica 2 arbitrarily)
    test_ground_truth = [test_ground_truth;DATA.Stained(k,2).DATA(1:N_Cells,Target_Ind)]; % Measured labeled targets (we chose replica 2 arbitrarily).
    test_data_U = [test_data_U;DATA.Unstained(k,2).DATA(1:N_Cells,Feature_Inds)]; % (we chose replica 2 arbitrarily)
end

%% Machine Leanrning Using LSBoost Regression Ensemble
i=1;
% internal cross-validation
cp_int = cvpartition(train_ground_truth,'k',5);
mstop_int = 0;

for j=1:5
    
    train_train_ind = cp_int.training(j);
    train_valid_ind = cp_int.test(j);
    train_train_data = train_data(train_train_ind,:); % Labaled features for training (internal)
    train_valid_data = train_data(train_valid_ind,:); % Labaled features for testing (internal)
    train_train_ground_truth = train_ground_truth(train_train_ind); % Ground truth for training (internal)
    train_valid_ground_truth = train_ground_truth(train_valid_ind); % Ground truth for testing (internal)
    
    LSTree_int = fitensemble(train_train_data,train_train_ground_truth,'LSBoost',max_trees,'Tree','LearnRate',0.1);
    
    reg_error_int{i}{j} = loss(LSTree_int,train_valid_data,train_valid_ground_truth,'mode','cumulative');
    mstop_int = mstop_int + find(reg_error_int{i}{j}==min(reg_error_int{i}{j}),1,'first');
    
    LSTree_int = [];
end

mstop(i)=round(mstop_int./5);

LSTree = fitensemble(train_data,train_ground_truth,'LSBoost',mstop(i),'Tree','LearnRate',0.1);

reg_error{i} = loss(LSTree,train_data,train_ground_truth,'mode','cumulative');

estimated_lipid_content_S_TRAIN = predict(LSTree,train_data);
estimated_lipid_content_S_TEST = predict(LSTree,test_data);
estimated_lipid_content_U_TRAIN = predict(LSTree,train_data_U);
estimated_lipid_content_U_TEST = predict(LSTree,test_data_U);
corr_S_train = corrcoef(estimated_lipid_content_S_TRAIN,train_ground_truth)
corr_S_test = corrcoef(estimated_lipid_content_S_TEST,test_ground_truth)
corr_U_train = corrcoef(estimated_lipid_content_U_TRAIN,train_ground_truth)
corr_U_test = corrcoef(estimated_lipid_content_U_TEST,test_ground_truth)

save('BLASSI_Results','est*','train_ground_truth','test_ground_truth','corr*'); % Save in the current directory
save('../../Data_Files/Data_Thomas_Blasi_Cell_Cycle_Analysis/BLASSI_Results','est*','train_ground_truth','test_ground_truth','corr*'); % Save in the 'Data_Files' folder
