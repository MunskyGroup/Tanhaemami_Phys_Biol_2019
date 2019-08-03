function [XtrainL,YtrainL,XtestL,YtestL,XtestUL] = Data_Collector(NCells,FeatureInds,TargetInd,TrainSet,TestSet,TrainRep,TestLRep,TestULRep,DATA)
XtrainL = [];   %Labeled training data
YtrainL = [];   %Measured targets
XtestL = [];    %Lebeled testing data
YtestL = [];    %Labeled target
XtestUL = [];   %Unlabeled testing data

for k = TrainSet  % Run through the different sets of data for training(each set corresponds to a time point).
    XtrainL = [XtrainL;DATA.Stained(k,TrainRep).DATA(1:NCells,FeatureInds)];   % Labeled features
    YtrainL = [YtrainL;DATA.Stained(k,TrainRep).DATA(1:NCells,TargetInd)];     % Labeled targets
end

for k = TestSet  % Run through the different sets of data for testing.
    XtestL = [XtestL;DATA.Stained(k,TestLRep).DATA(NCells+1:NCells*2,FeatureInds)];   % Labeled features
    YtestL = [YtestL;DATA.Stained(k,TestLRep).DATA(NCells+1:NCells*2,TargetInd)];     % Labeled targets
    XtestUL = [XtestUL;DATA.Unstained(k,TestULRep).DATA(NCells+1:NCells*2,FeatureInds)];   % Unlabeled features
end
end