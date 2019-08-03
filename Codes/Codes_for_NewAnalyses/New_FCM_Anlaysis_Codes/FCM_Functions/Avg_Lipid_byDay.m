function Lipid_Average = Avg_Lipid_byDay(TargetInd,TimePoint,TrainSet,ValidSet,ContainsReplicates,DATA)
%Plot lipid content averages for each day:
if strcmp(ContainsReplicates,'True')
    for m=1:size(DATA.Stained,2)
        figure(2);
        MU = [];
        err = [];
        [~,Date_Inds] = sort(DATA.Time);
        for k=Date_Inds
            MU = [MU;mean(DATA.Stained(k,m).DATA(:,TargetInd))];
            err = [err;std(DATA.Stained(k,m).DATA(:,TargetInd))];
        end
        kk = 1:size(MU,1);
        hold on
        Lipid_Average = errorbar(kk,MU,err,'--o');
        set(gca,'XTick',1:TimePoint);
        assignin('base','Average_Lipid_Content',MU);
    end
    plot(kk(TrainSet),MU(TrainSet),'o','MarkerFaceColor','r','MarkerSize',10);
    plot(kk(ValidSet),MU(ValidSet),'o','MarkerFaceColor','b','MarkerSize',10);
    xlabel('Time (days)');
    ylabel('Average Lipid Content (AUC)');
    legend('Replica 1','Replica 2','Replica 3','Training','Validation','Location','NorthWest');
    title(['Average Lipid Content for all ',num2str(m),' Replicates']);
else
    MU = [];
    err = [];
    figure(2);
    [~,Date_Inds] = sort(DATA.Time);
    for k=Date_Inds
        MU = [MU;mean(DATA.Stained(k).DATA(:,TargetInd))];
        err = [err;std(DATA.Stained(k).DATA(:,TargetInd))];
    end
    kk = 1:size(MU,1);
    hold on
    Lipid_Average = errorbar(kk,MU,err,'--o');
    set(gca,'XTick',1:TimePoint);
    assignin('base','Average_Lipid_Content',MU);
    plot(kk(TrainSet),MU(TrainSet),'o','MarkerFaceColor','r','MarkerSize',10);
    plot(kk(ValidSet),MU(ValidSet),'o','MarkerFaceColor','b','MarkerSize',10);
    xlabel('Time (days)');
    legend('Means','Training','Validation');
    title('Average Lipid Content');
end
saveas(figure(2),'Figures_FINAL/Average_Lipid_Content.pdf');
end