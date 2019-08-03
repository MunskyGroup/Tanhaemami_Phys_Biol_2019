function [XTestL,YTestL,XTestUL,XTestUL_S] = DataCollector_Test(NCells,FeatureInds,FeatureInds_S,TargetInd,TestSet,TestLRep,TestULRep,DATA)
%Author: Mohammad Tanhaemami
%Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.
% This function collects the data for testing purposes
i=1;
[XTestL,YTestL,XTestUL,XTestUL_S] = deal(cell(size(TestSet,2),1));
for k = TestSet  % Run through the different sets of data for testing.
    XTestL{i} = [XTestL{i};DATA.Stained(k,TestLRep).DATA(NCells+1:NCells*2,FeatureInds)];   % Labeled features
    YTestL{i} = [YTestL{i};DATA.Stained(k,TestLRep).DATA(NCells+1:NCells*2,TargetInd)];     % Labeled targets
    XTestUL{i} = [XTestUL{i};DATA.Unstained(k,TestULRep).DATA(NCells+1:NCells*2,FeatureInds)];   % Unlabeled features
    XTestUL_S{i} = [XTestUL_S{i};DATA.Unstained(k,TestULRep).DATA(NCells+1:NCells*2,FeatureInds_S)];   % Unlabeled features for S
    i = i + 1;
end
end