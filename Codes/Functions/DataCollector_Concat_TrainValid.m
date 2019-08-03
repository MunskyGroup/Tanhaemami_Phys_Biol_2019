function [XTrainL,XTrainUL,XTrainUL_S,YTrainL,XValidL,XValidUL,XValidUL_S,YValidL] = DataCollector_Concat_TrainValid(NCellsL,NCellsUL,FeatureInds,FeatureInds_S,TargetInd,TrainSet,TestSet,DATA)
%Author: Mohammad Tanhaemami
%Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.
% This function collects the data for training and validation purposes.

i = 1;
[XTrainL,XTrainUL,XTrainUL_S,YTrainL,XValidL,XValidUL,XValidUL_S,YValidL] = deal(cell(size(TrainSet,2),1));
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
for k = TestSet  % Run through the different sets of data for validation.
    for TestLRep = 1:size(DATA.Stained,2)
        XValidL{i} = [XValidL{i};DATA.Stained(k,TestLRep).DATA(NCellsL+1:NCellsL*2,FeatureInds)];   % Labeled features of validation data
        YValidL{i} = [YValidL{i};DATA.Stained(k,TestLRep).DATA(NCellsL+1:NCellsL*2,TargetInd)];     % Labeled targets of validation data
    end
    i = i +1;
end
i = 1;
for k = TestSet  % Run through the different sets of data for validation.
    for TestULRep = 1:size(DATA.Unstained,2)
        XValidUL{i} = [XValidUL{i};DATA.Unstained(k,TestULRep).DATA(NCellsUL+1:NCellsUL*2,FeatureInds)];   % Unlabeled features of validation data
        XValidUL_S{i} = [XValidUL_S{i};DATA.Unstained(k,TestULRep).DATA(NCellsUL+1:NCellsUL*2,FeatureInds_S)];   % Unlabeled features of validation data for S
    end
    i = i +1;
end
end