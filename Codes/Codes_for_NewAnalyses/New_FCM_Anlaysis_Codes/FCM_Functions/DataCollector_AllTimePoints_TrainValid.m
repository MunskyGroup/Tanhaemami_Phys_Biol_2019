function [XTrainL,YTrainL,XValidL,YValidL,XValidUL] = DataCollector_AllTimePoints_TrainValid(NCellsL,NCellsUL,FeatureInds,TargetInd,TrainSet,TestSet,DATA)
i = 1;
[XTrainL,YTrainL,XValidL,YValidL,XValidUL] = deal(cell(size(TrainSet,2),1));
%TRAINING: (cells 1:3000)
for k = TrainSet  % Run through the different sets of data for training(each set corresponds to a time point).
    for TrainRep = 1:size(DATA.Stained,2)
        XTrainL{i} = [XTrainL{i};DATA.Stained(k,TrainRep).DATA(1:NCellsL,FeatureInds)];   % Labeled features
        YTrainL{i} = [YTrainL{i};DATA.Stained(k,TrainRep).DATA(1:NCellsL,TargetInd)];     % Labeled targets
    end
    i = i + 1;
end
%VALIDATION: (cells 3001:6000)
i = 1;
for k = TestSet  % Run through the different sets of data for validation.
    for TestLRep = 1:size(DATA.Stained,2)
        XValidL{i} = [XValidL{i};DATA.Stained(k,TestLRep).DATA(NCellsL+1:NCellsL*2,FeatureInds)];   % Labeled features
        YValidL{i} = [YValidL{i};DATA.Stained(k,TestLRep).DATA(NCellsL+1:NCellsL*2,TargetInd)];     % Labeled targets
    end
    i = i +1;
end
i = 1;
for k = TestSet  % Run through the different sets of data for validation.
    for TestULRep = 1:size(DATA.Unstained,2)
        XValidUL{i} = [XValidUL{i};DATA.Unstained(k,TestULRep).DATA(NCellsUL+1:NCellsUL*2,FeatureInds)];   % Unlabeled features
    end
    i = i +1;
end
end