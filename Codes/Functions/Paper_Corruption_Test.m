function [Stained,Unstained,KSDistance] = Paper_Corruption_Test(NCells,FeatureInds,SelectedTimePoint,StainedReplica,UnstainedReplica,DATA)
%Author: Mohammad Tanhaemami
%Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.
% This function finds possible corruptions in the features (labeled vs unlabeled) due to labeling
% Do the following for the selected time_point (timecount = time_point).
timepoint = SelectedTimePoint;
figure(8); clf; %Create figure, clear it.
for i=1:length(FeatureInds) %For the 11 features (13 total. 3 and 9 are targets.)
    try
        subplot(3,6,i)
        Stained=log(DATA.Stained(timepoint,StainedReplica).DATA(1:NCells,FeatureInds(i)));
        ksdensity(Stained(Stained~=Inf & (Stained~=-Inf))); %Plot a kernel density of stained data.
        hold on
        Unstained=log(DATA.Unstained(timepoint,UnstainedReplica).DATA(1:NCells,FeatureInds(i)));
        ksdensity(Unstained(Unstained~=Inf & (Unstained~=-Inf))); %Plot a kernel density of unstained data.
        title(DATA.FeatureNames(FeatureInds(i)),'FontSize',11)
        set(gca,'xlim',([min(min(Stained),min(Unstained)) max(max(Stained),max(Unstained))]))
    catch
    end
    %Create barchart for each timepoint with KS test method
    PL = DATA.Stained(timepoint,StainedReplica).DATA(1:NCells,FeatureInds(i));
    N_PL = length(PL);
    PUL = DATA.Unstained(timepoint,UnstainedReplica).DATA(1:NCells,FeatureInds(i));
    N_PUL = length(PUL);
    PLUL = sort([PUL;PL]);
    b_PL = zeros(length(PLUL),1);      % cumulative dist of 'PL'
    b_PUL = zeros(length(PLUL),1);      % cumulative dist of 'PUL'
    
    for count = 1:length(PLUL)
        b_PL(count) = sum(PL<PLUL(count))/N_PL;
        b_PUL(count) = sum(PUL<PLUL(count))/N_PUL;
    end
    KSDistance(i) = max(abs(b_PL-b_PUL));
end
lgnd = legend('Labeled','Unlabeled');
lgnd.Position = [0.83 0.41 0.06 0.05];
lgnd.FontSize = 7;
subplot(3,6,13:18);
bar(KSDistance);
hold on
ZZ = get(gca,'xlim'); %Get limits of the x axes of the current axes.
plot(ZZ,[0.5,0.5],'k--') %Ideal line of correlation.
BarchartInfo = get(gca);
xlabel('Features','FontSize',10)
ylabel('Distance','FontSize',10)
BarchartInfo.XAxis.TickLabels = DATA.FeatureNames(FeatureInds);
BarchartInfo.XAxis.TickLabelRotation = 45;
end
