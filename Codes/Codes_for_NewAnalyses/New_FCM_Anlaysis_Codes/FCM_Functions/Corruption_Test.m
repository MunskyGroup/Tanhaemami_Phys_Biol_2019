function [Stained,Unstained,KSDistance] = Corruption_Test(Step,NCells,FeatureInds,SelectedTimePoint,StainedReplica,UnstainedReplica,DATA)
close all; %Close all previous figures to avoide any conflict in numberieng figures in this section.

% Do the following for the selected time_point (timecount = time_point).
for timecount = SelectedTimePoint %Count from the first timepoint to the last.
    figure(timecount); clf; %Create figure, clear it.
    for i=1:length(FeatureInds) %For the 11 features (13 total. 3 and 9 are targets.)
        try
            subplot(3,6,i)
            Stained=log(DATA.Stained(timecount,StainedReplica).DATA(1:NCells,FeatureInds(i)));
            ksdensity(Stained(Stained~=Inf & (Stained~=-Inf))); %Plot a kernel density of stained data.
            hold on
            Unstained=log(DATA.Unstained(timecount,UnstainedReplica).DATA(1:NCells,FeatureInds(i)));
            ksdensity(Unstained(Unstained~=Inf & (Unstained~=-Inf))); %Plot a kernel density of unstained data.
            title(DATA.FeatureNames(FeatureInds(i)),'FontSize',11)
            set(gca,'xlim',([min(min(Stained),min(Unstained)) max(max(Stained),max(Unstained))]))
        catch
        end
        %Create barchart for each timepoint with KS test method
        PL = DATA.Stained(timecount,StainedReplica).DATA(1:NCells,FeatureInds(i));
        N_PL = length(PL);
        PUL = DATA.Unstained(timecount,UnstainedReplica).DATA(1:NCells,FeatureInds(i));
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
    BarchartInfo(timecount) = get(gca);
    xlabel('Features','FontSize',10)
    ylabel('Distance','FontSize',10)
    BarchartInfo(timecount).XAxis.TickLabels = DATA.FeatureNames(FeatureInds);
    BarchartInfo(timecount).XAxis.TickLabelRotation = 45;
    saveas(figure(timecount),['Figures/Step',num2str(Step),'_timepoint_',num2str(timecount),'.fig']); %Saving final figure
end
saveas(figure(timecount),['Figures/Step',num2str(Step),'_FINAL.pdf']); %Saving final figure as a *.pdf
end
%% Figure 2 - With KS test - Commented Out!!!
% %Here we redo figure 2 but with 2-sample Kolmogorov–Smirnov(KS) test.
% %This test determines whether stained and unstained data come from the same distribution.
% %First we gather stained and unstained data and plot empirical cumulative distribution functions of them to visualize the distributions.
% %Then we perform 2-sample KS test on each stained dataset with respect to their corresponsing unstained dataset.
% %KS-test results are saved in worksapce (95% CI).
% 
% %PL: Indicates probability distribution for labeled(stained) data
% %PUL: Indicates probability distribution for unlabeled(unstained) data
% %ECDF_L: Indicates empirical cumulative distribution function for PL
% %ECDF_UL: Indicates empirical cumulative distribution function for PUL
% 
% N_Cells = 3000;
% Feature_Inds = [1 2 4 5 6 7 8 10 11 12 13];
% for time_point = 1:13 %We have 13 time points
%     figure(time_point); clf; %Create figure, clear it.
%     for i=1:length(Feature_Inds) %For the 11 features (13 total. 3 and 9 are targets.)
%         try
%             subplot(3,4,i)
%             PL=DATA.Stained(time_point).DATA(1:N_Cells,Feature_Inds(i))/sum(DATA.Stained(time_point).DATA(1:N_Cells,Feature_Inds(i)));
%             [ECDF_L,PL_values] = ecdf(PL);
%             plot(PL_values,ECDF_L);
%             hold on
%             PUL=DATA.Unstained(time_point).DATA(1:N_Cells,Feature_Inds(i))/sum(DATA.Unstained(time_point).DATA(1:N_Cells,Feature_Inds(i)));
%             [ECDF_UL,PUL_values] = ecdf(PUL);
%             plot(PUL_values,ECDF_UL);
%             if length(ECDF_L)>length(ECDF_UL)
%                 ks_distance(i) = max(abs(ECDF_L(1:length(ECDF_UL))-ECDF_UL));
%             elseif length(ECDF_L)<length(ECDF_UL)
%                 ks_distance(i) = max(abs(ECDF_L-ECDF_UL(1:length(ECDF_L))));
%             else
%                 ks_distance(i) = max(abs(ECDF_L-ECDF_UL));
%             end
%             title(sprintf('Feature %d Max Diff = %1.4f',Feature_Inds(i),ks_distance(i)),'FontSize', 10)
%             xlim([0 0.001])
%             legend('Stained','Unstained','Location','SouthEast')
%             [ks_result(time_point).ks(i),ks_result(time_point).pvalue(i)] = kstest2(PL,PUL);
%         catch
%         end
%     end
%     %saveas(figure(time_point),sprintf('Fig2_KS_timepoint_%d',time_point));
% end
