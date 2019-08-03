function [h] = Figure_Generator(Fig_Num,XtrainL,YtrainL,XtestL,YtestL,XtestUL,BINS)
%This function operates with KS test as default. For other metrics of error estimation the input for "Metric" in function "Simple_Linear_Regression" should change to other values.
if ~exist(['Results/Map',num2str(Fig_Num),'.mat'],'file') %Check whether the model (M and KS's) exists.
    [M] = Simple_Linear_Regression(Fig_Num,XtrainL,YtrainL,'true','KS');
    %Learn from training. Predict the same data.
    [KS_1] = Simple_Linear_Regression_Predict(M,XtrainL,YtrainL,'false','KS',BINS) %Perform a simple linear regression with a histogram of Ytest and Xtest*M
    %Learn from training. Predict with labeled testing features(X*M). Compare with labeled testing targets.
    [KS_2] = Simple_Linear_Regression_Predict(M,XtestL,YtestL,'false','KS',BINS) %Perform a simple linear regression with a histogram of Ytest and Xtest*M
    %Learn from training. Predict with unlabeled testing features(X*M). Compare with labeled testing targets.
    [KS_3] = Simple_Linear_Regression_Predict(M,XtestUL,YtestL,'false','KS',BINS) %Perform a simple linear regression with a histogram of Ytest and Xtest*M
    save(['Results/Map',num2str(Fig_Num),'.mat'],'M','KS_1','KS_2','KS_3'); %save data in Map1.mat in results folder
else
%     q=input('This Model Already Saved. Remake figures? (y/n)','s');
%     if strcmp(q,'y')
%         [M] = Simple_Linear_Regression(XtrainL,YtrainL,'true','KS'); %Perform a simple linear regression with a histogram of Ytest and Xtest*M
%         [KS_1] = Simple_Linear_Regression_Predict(M,XtrainL,YtrainL,'true','KS',BINS); %Perform a simple linear regression with a histogram of Ytest and Xtest*M
%         %Learn from training. Predict with labeled testing features(X*M). Compare with labeled testing targets.
%         [KS_2] = Simple_Linear_Regression_Predict(M,XtestL,YtestL,'true','KS',BINS); %Perform a simple linear regression with a histogram of Ytest and Xtest*M
%         %Learn from training. Predict with unlabeled testing features(X*M). Compare with labeled testing targets.
%         [KS_3] = Simple_Linear_Regression_Predict(M,XtestUL,YtestL,'true','KS',BINS); %Perform a simple linear regression with a histogram of Ytest and Xtest*M
%         save(['Results/Map',num2str(Fig_Num),'.mat'],'M','KS_1','KS_2','KS_3'); %save data in Map1.mat in results folder
%     else
        disp('Map exists in "Results" folder. Skipping figures');
        h=[];
        return
%     end
end
% figure;
%Saving Subplots 2, 3, and 4 of the final figure (subplot 1 is saved by Simple_Linear_Regression.m).
for k=0:3
    saveas(figure(2*k+1),['Figures/Fig',num2str(Fig_Num),'_',num2str(k+1),'.fig']);
end

%The 4 subplots of figure1 are saved separately. Showing above figures all in one.
% Load handels to the saved figures without showing them
Subplot1 = hgload(['Figures/Fig',num2str(Fig_Num),'_1.fig'],struct('Visible','off'));
Subplot2 = hgload(['Figures/Fig',num2str(Fig_Num),'_2.fig'],struct('Visible','off'));
Subplot3 = hgload(['Figures/Fig',num2str(Fig_Num),'_3.fig'],struct('Visible','off'));
Subplot4 = hgload(['Figures/Fig',num2str(Fig_Num),'_4.fig'],struct('Visible','off'));
% Prepare subplots
figure;
title_names = {'(a) Training - Labeled' '(b) Training - Labeled' '(c) Testing - Labeled' '(d) Testing - Unlabeled'};
xlabel_names = {'Measured' 'Measured Lipid Content (AUC)' 'Measured Lipid Content (AUC)' 'Measured Lipid Content (AUC)'};
ylabel_names = {'Predicted' '' '' ''};
for i=1:4 %Create subplots and set axes and title for each
    h(i)=subplot(2,2,i);
    SubplotsInfo(i) = get(gca);
    title(title_names{i},'FontSize',11)
    xlabel(xlabel_names{i},'FontSize',10)
    ylabel(ylabel_names{i},'FontSize',10)
end
for i=2:4 %Set axis limits for subplots 2,3,and 4
    h(i)=subplot(2,2,i);
    SubplotsInfo(i).XAxis.Limits = [0 3e6];
    SubplotsInfo(i).YAxis.Limits = [0 0.02];
    Distance = sprintf('KS_%d',i-1);
    text(1.1e6,0.011,sprintf('D(Measured,Predicted) = %1.4f',eval(Distance)),'FontSize',9);
end
% Paste figures on the subplots.
for i=1:length(h)
    copyobj(allchild(get(eval(sprintf('Subplot%d',i)),'CurrentAxes')),h(i));
end
%Add legend.
legend(h(1),'Training','Ideal','Location','SouthEast')
for i=2:length(h)
    legend(h(i),'Measured','Predicted','Location','NorthEast')
end

saveas(figure(12),['Figures/Fig',num2str(Fig_Num),'_FINAL.fig']); %Saving final figure
saveas(figure(12),['Figures/Fig',num2str(Fig_Num),'_FINAL.pdf']); %Saving final figure as a *.pdf
end