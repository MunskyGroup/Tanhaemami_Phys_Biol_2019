function [XTrainL,XTrainUL,XTrainUL_S,YTrainL,XValidUL,XValidUL_S,YValidL] = Alpha_CollectData(NCellsL,NCellsUL,FeatureInds,FeatureInds_S,FeatsType,TargetInd,TrainSet,ValidSet,DATA)
%Author: Mohammad Tanhaemami
%Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.

% Collects Data for training and validation timepoints
% Combines replicates at each timepoint, separately. For example, Timepoint 1: [Rep1;Rep2;Rep3] and so on.
% This function can take care of quadratic features too.

i = 1;
[XTrainL,XTrainUL,XTrainUL_S,YTrainL,YValidL,XValidUL,XValidUL_S] = deal(cell(size(TrainSet,2),1));
%TRAINING:
for k = TrainSet  % Run through the different sets of data for training(each set corresponds to a time point).
    for TrainRep = 1:size(DATA.Stained,2)
        XTrainL{i} = [XTrainL{i};DATA.Stained(k,TrainRep).DATA(1:NCellsL,FeatureInds)];   % Labeled features of training data
        YTrainL{i} = [YTrainL{i};DATA.Stained(k,TrainRep).DATA(1:NCellsL,TargetInd)];     % Labeled targets of training data
    end
    i = i + 1;
end
i=1;
for k = TrainSet  % Run through the different sets of data for training(each set corresponds to a time point).
    for TrainRep = 1:size(DATA.Unstained,2)
        XTrainUL{i} = [XTrainUL{i};DATA.Unstained(k,TrainRep).DATA(1:NCellsUL,FeatureInds)];   % Unlabeled features of training data
        XTrainUL_S{i} = [XTrainUL_S{i};DATA.Unstained(k,TrainRep).DATA(1:NCellsUL,FeatureInds_S)];   % Unlabeled features of training data for S
    end
    i = i + 1;
end
%VALIDATION:
i = 1;
for k = ValidSet  % Run through the different sets of data for validation.
    for TestLRep = 1:size(DATA.Stained,2)
        YValidL{i} = [YValidL{i};DATA.Stained(k,TestLRep).DATA(1:NCellsL,TargetInd)];     % Labeled targets of validation data
    end
    i = i +1;
end
i = 1;
for k = ValidSet  % Run through the different sets of data for validation.
    for TestULRep = 1:size(DATA.Unstained,2)
        XValidUL{i} = [XValidUL{i};DATA.Unstained(k,TestULRep).DATA(1:NCellsUL,FeatureInds)];   % Unlabeled features of validation data
        XValidUL_S{i} = [XValidUL_S{i};DATA.Unstained(k,TestULRep).DATA(1:NCellsUL,FeatureInds_S)];   % Unlabeled features of validation data for S
    end
    i = i +1;
end
if strcmp(FeatsType,'Quadratic') %If the input features include quadratic terms, calculate and add the quadratic values to the original matrix of features.
    [XTrainLQ,XTrainULQ,XTrainULQ_S,XValidULQ,XValidULQ_S] = deal({});
    for i = 1:size(XTrainL,1)
        XTrainLQ{i,1} = Get_Quad_Features(XTrainL{i,1}); %Getting quadratic X1 with the function Get_Quad_Features.
        XTrainULQ{i,1} = Get_Quad_Features(XTrainUL{i,1});
        XTrainULQ_S{i,1} = Get_Quad_Features(XTrainUL_S{i,1});
        XValidULQ{i,1} = Get_Quad_Features(XValidUL{i,1}); %Getting quadratic XtestUL with the function Get_Quad_Features.
        XValidULQ_S{i,1} = Get_Quad_Features(XValidUL_S{i,1}); %Getting quadratic XtestUL with the function Get_Quad_Features for S.
    end
    XTrainL = XTrainLQ;
    XTrainUL = XTrainULQ;
    XTrainUL_S = XTrainULQ_S;
    XValidUL = XValidULQ;
    XValidUL_S = XValidULQ_S;
end
end