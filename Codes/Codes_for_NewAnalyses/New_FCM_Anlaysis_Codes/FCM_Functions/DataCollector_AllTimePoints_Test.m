function [XTestL,YTestL,XTestUL] = DataCollector_AllTimePoints_Test(NCells,FeatureInds,TargetInd,TestSet,TestLRep,TestULRep,DATA)
i=1;
[XTestL,YTestL,XTestUL] = deal(cell(size(TestSet,2),1));
for k = TestSet  % Run through the different sets of data for testing. (cells 6001:9000)
    XTestL{i} = [XTestL{i};DATA.Stained(k,TestLRep).DATA(NCells*2+1:NCells*3-41,FeatureInds)];   % Labeled features (41 because of shortage on labeled data in timepoint 23)
    YTestL{i} = [YTestL{i};DATA.Stained(k,TestLRep).DATA(NCells*2+1:NCells*3-41,TargetInd)];     % Labeled targets
    XTestUL{i} = [XTestUL{i};DATA.Unstained(k,TestULRep).DATA(NCells*2+1:NCells*3-41,FeatureInds)];   % Unlabeled features
    i = i + 1;
end
end