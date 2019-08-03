classdef FCSC_Log
    % Author: Mohammad Tanhaemami
    % Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.
    % This object is used for developing the regular strategy with respect to a single averaged modeling approach
    properties
        Step % Step of the calculation according to the main script
        Stage % Stage of computations (training or testing) used for generating figures
        BINS % Includes information regarding the bins used for the histograms
        XTrainL % Labeled training features
        XTrainUL % Unlabeled training features
        YTrainL % Labeled training targets
        XTestL % Labeled validation (or testing) features
        XTestUL % Unlabeled validation (or testing) features
        YTestL % Labeled validation (or testing) targets
        I % Vector indicating which features to be used as a result of applying the genetic algorithm
        MakeHistPlot % Logical to whether make histogram plots
        GenerateFigure % Logical to whether generate figures
    end
    
    properties(Constant)
    end
    
    properties(Dependent)
        RegressModel % Output of the object
    end
    %%
    methods
        function Model = get.RegressModel(obj)
            Model.Step = obj.Step;
            % Initial action: Convert input data values from linear to log scale.
            [Model.LogXTrainLAll,Model.LogYTrainLAll,Model.LogXTestLAll,Model.LogYTestLAll,Model.LogXTestULAll] = deal([]);
            
            %Training(concatenated):
            for i=1:size(obj.XTrainL,1)
                Model.LogXTrainL{i,1} = log(obj.XTrainL{i,1}); % Labeled training features
                Model.LogXTrainUL{i,1} = log(obj.XTrainUL{i,1}); % Unlabeled training features
                Model.LogYTrainL{i,1} = log(obj.YTrainL{i,1}); % Labeled training targets
                Model.LogXTrainLAll = [Model.LogXTrainLAll;Model.LogXTrainL{i,1}]; % Concatenate the labeled training features
                Model.LogYTrainLAll = [Model.LogYTrainLAll;Model.LogYTrainL{i,1}]; % Concatenate the labeled training targets
            end
            Model.MAll = Model.LogXTrainLAll(:,obj.I==1)\Model.LogYTrainLAll; % Regression coefficient %%%LOG-NORMAL ERRORS%%%
            Model.LogYPredTrainLAll = Model.LogXTrainLAll(:,obj.I==1)*Model.MAll; % Predicted labeled targets (concatenated)
            Model.CorrCoefTrainAll = corrcoef(Model.LogYPredTrainLAll,Model.LogYTrainLAll); % Compute correlation coefficient
            [~,~,Model.ErrTrainAll] = kstest2(Model.LogYTrainLAll,Model.LogYPredTrainLAll); % Compute the KS distance
            
            %Training(prediction for serparated timepoints):
            for i = 1:size(obj.XTrainL,1)
                Model.LogYPredTrainL{i,1} = Model.LogXTrainL{i,1}(:,obj.I==1)*Model.MAll; % Predict labeled targets indiviually with repect to the 3 time points
                [~,~,Model.ErrTrainL{i,1}] = kstest2(Model.LogYTrainL{i,1},Model.LogYPredTrainL{i,1}); % Compute the KS distance
            end
            
            %Validation(or Testing)(separated):
            for i=1:size(obj.XTestL,1) %i is the counter for timepoints
                Model.LogXTestL{i,1} = log(obj.XTestL{i,1}); % Labeled validation (or testing) features
                Model.LogYTestL{i,1} = log(obj.YTestL{i,1}); % Labeled validation (or testing) targets
                %Validation(labeled):
                Model.LogYPredTestL{i,1} = Model.LogXTestL{i,1}(:,obj.I==1)*Model.MAll; % Predict the labeled targets
                Model.CorrCoefTestL{i,1} = corrcoef(Model.LogYPredTestL{i,1},Model.LogYTestL{i,1}); % Correlation coefficients
                [~,~,Model.ErrTestL{i,1}] = kstest2(Model.LogYTestL{i,1},Model.LogYPredTestL{i,1}); % Compute the KS distance
                %Validation(unlabeled):
                Model.LogXTestUL{i,1} = log(obj.XTestUL{i,1});  % unlabeled validation (or testing) features
                Model.LogYPredTestUL{i,1} = Model.LogXTestUL{i,1}(:,obj.I==1)*Model.MAll; % Predict the unlabaled targets
                [~,~,Model.ErrTestUL{i,1}] = kstest2(Model.LogYTestL{i,1},Model.LogYPredTestUL{i,1}); % The KS distance
                %Combine timepoints:
                Model.LogXTestLAll = [Model.LogXTestLAll;Model.LogXTestL{i,1}]; % Concatenated labeled validation (or testing) features
                Model.LogYTestLAll = [Model.LogYTestLAll;Model.LogYTestL{i,1}]; % Concatenated labeled validation (or testing) targets
                Model.LogXTestULAll = [Model.LogXTestULAll;Model.LogXTestUL{i,1}]; % Concatenated unlabeled validation (or testing) features
            end
            %Validation(or Testing)(concatenated):
            Model.LogYPredTestLAll = Model.LogXTestLAll(:,obj.I==1)*Model.MAll; % Predict labeled validation (or testing) targets
            Model.CorrCoefTestLAll = corrcoef(Model.LogYPredTestLAll,Model.LogYTestLAll); % Correlation coefficients
            [~,~,Model.ErrTestLAll] = kstest2(Model.LogYTestLAll,Model.LogYPredTestLAll); % KS distance
            Model.LogYPredTestULAll = Model.LogXTestULAll(:,obj.I==1)*Model.MAll; % Predict unlabeled validation (or testing) targets
            [~,~,Model.ErrTestULAll] = kstest2(Model.LogYTestLAll,Model.LogYPredTestULAll); % KS distance
            
            %Histograms:
            %Histograms on top of each other for each step:
            if strcmp(obj.MakeHistPlot,'True')
                Model.HistTrain = FCSC_Log.Hist_Plot('Combined',Model.MAll,Model.LogXTrainLAll,Model.LogYTrainLAll,obj.I,obj.BINS);
                title('HistTrain');
                Model.HistTestL = FCSC_Log.Hist_Plot('SeparatedTest',Model.MAll,Model.LogXTestL,Model.LogYTestL,obj.I,obj.BINS);
                title('HistValidL');
                Model.HistTestUL = FCSC_Log.Hist_Plot('SeparatedTest',Model.MAll,Model.LogXTestUL,Model.LogYTestL,obj.I,obj.BINS);
                title('HistValidUL');
                Model.HistTrainLAll = FCSC_Log.Hist_Plot('Combined',Model.MAll,Model.LogXTrainLAll,Model.LogYTrainLAll,obj.I,obj.BINS);
                title('HistTrainAll');
                Model.HistTestLAll = FCSC_Log.Hist_Plot('Combined',Model.MAll,Model.LogXTestLAll,Model.LogYTestLAll,obj.I,obj.BINS);
                title('HistValidLAll');
                Model.HistTestULAll = FCSC_Log.Hist_Plot('Combined',Model.MAll,Model.LogXTestULAll,Model.LogYTestLAll,obj.I,obj.BINS);
                title('HistValidULAll');
                if strcmp(obj.GenerateFigure,'True')
                    Model.FigureSep = FCSC_Log.FigureGenerator_Sep(obj.Stage,size(obj.XTestL,1),obj.Step,Model.ErrTrainAll,Model.ErrTestL,Model.ErrTestUL);
                    Model.FigureComb = FCSC_Log.FigureGenerator_Comb(obj.Stage,size(obj.XTestL,1),obj.Step,Model.ErrTrainAll,Model.ErrTestLAll,Model.ErrTestULAll);
                end
            end
        end
    end
    %%
    methods(Static)
        
        function AvgKSDist = Compute_SepKS(XTrain,YTrain,XTest,YTest,I) % Calculates KS distance b/w measurements and predictions
            [SepKSDist,LogYPred] = deal({});
            [LogXTrainAll,LogYTrainAll] = deal([]);
            [LogXTest,LogYTest] = deal({});
            for j = 1:size(XTrain,1)
                LogXTrainAll = [LogXTrainAll;log(XTrain{j,1})];
                LogYTrainAll = [LogYTrainAll;log(YTrain{j,1})];
            end
            M = LogXTrainAll(:,I==1)\LogYTrainAll; % Computhe the model
            for j = 1:size(XTrain,1)
                LogXTest{j,1} = log(XTest{j,1}); % Validation (or testing) features
                LogYTest{j,1} = log(YTest{j,1}); % Validation (or testing) targets
                LogYPred{j,1} = LogXTest{j,1}(:,I==1)*M; % Predict targets
                [~,~,SepKSDist{j,1}] = kstest2(LogYTest{j,1},LogYPred{j,1}); % KS distance
            end
            SepKSDist = cell2mat(SepKSDist);
            AvgKSDist = mean(SepKSDist);
        end
        
        function Make_Hist_Plot = Hist_Plot(TimepointsFormat,M,X,Y,I,BINS)
            if strcmp(TimepointsFormat,'SeparatedTrain')
                for i = 1:size(X,1)
                    figure; clf;
                    A{i,1} = hist(Y{i,1},BINS)/sum(hist(Y{i,1},BINS)); %Take histogram values of measured target
                    B{i,1} = hist(X{i,1}(:,I==1)*M{i,1},BINS)/sum(hist(X{i,1}(:,I==1)*M{i,1},BINS)); %Take histogram values of predicted values
                    Make_Hist_Plot{i,1} = plot(BINS',[A{i,1};B{i,1}]','LineWidth',0.5); %Plot hist
                    set(gca,'fontsize',16)
                    legend('Measured 1','Predicted 1','Measured 2','Predicted 2','Measured 3','Predicted 3');
                end
            elseif strcmp(TimepointsFormat,'SeparatedTest')
                for i = 1:size(X,1)
                    figure; clf;
                    A{i,1} = hist(Y{i,1},BINS)/sum(hist(Y{i,1},BINS)); %Take histogram values of measured target
                    B{i,1} = hist(X{i,1}(:,I==1)*M,BINS)/sum(hist(X{i,1}(:,I==1)*M,BINS)); %Take histogram values of predicted values
                    Make_Hist_Plot{i,1} = plot(BINS',[A{i,1};B{i,1}]','LineWidth',0.5); %Plot hist
                    set(gca,'fontsize',16)
                    legend('Measured 1','Predicted 1','Measured 2','Predicted 2','Measured 3','Predicted 3');
                end
            elseif strcmp(TimepointsFormat,'Combined')
                figure; clf;
                A = hist(Y,BINS)/sum(hist(Y,BINS)); %Take histogram values of measured target
                B = hist(X(:,I==1)*M,BINS)/sum(hist(X(:,I==1)*M,BINS)); %Take histogram values of predicted values
                Make_Hist_Plot = plot(BINS',[A;B]','LineWidth',0.5); %Plot hist
                set(gca,'fontsize',16)
                legend('Measured','Predicted')
            end
        end
        
        function [h] = FigureGenerator_Sep(Stage,Num_Timepoints,Step_Num,KS1,KS2,KS3)
            Num_Figs = 2*Num_Timepoints+1;
            if strcmp(Stage,'TrainValid')
                for k=1:Num_Figs
                    saveas(figure(k),['Results/Figures/Step',num2str(Step_Num),'_',num2str(k),'.fig']);
                    Subplot(k) = hgload(['Results/Figures/Step',num2str(Step_Num),'_',num2str(k),'.fig'],struct('Visible','off'));
                end
                KSTrain = KS1;
                for i=1:Num_Timepoints
                    KSTestL(i) = KS2{i,1};
                    KSTestUL(i) = KS3{i,1};
                end
                KSTest = [KSTestL KSTestUL];
                %Prepare subplots:
                orient(figure,'landscape')
                titleTestL = repmat({'Testing - Labeled'},1,Num_Timepoints);
                titleTestUL = repmat({'Testing - Unlabeled'},1,Num_Timepoints);
                title_names = [titleTestL,titleTestUL];
                h(1) = subplot(3,Num_Timepoints,1);
                SubplotsInfo(1) = get(gca);
                title('Training');
                xlabel('Lipid Content (AUC)','FontSize',10)
                ylabel('Probability','FontSize',10)
                SubplotsInfo(1).XAxis.Limits = [0 20]; %Set axis limits for subplots
                SubplotsInfo(1).YAxis.Limits = [0 0.03];
                text(1,0.015,sprintf('D_{KS} = %1.4f',eval(num2str(KSTrain))),'FontSize',9);
                for i=Num_Timepoints+1:3*Num_Timepoints %Create subplots and set axes and title for each
                    h(i)=subplot(3,Num_Timepoints,i);
                    SubplotsInfo(i) = get(gca);
                    title(title_names{1,i-Num_Timepoints},'FontSize',10)
                    xlabel('Lipid Content (AUC)','FontSize',10)
                    ylabel('Probability','FontSize',10)
                    SubplotsInfo(i).XAxis.Limits = [0 20]; %Set axis limits for subplots
                    SubplotsInfo(i).YAxis.Limits = [0 0.03];
                    text(1,0.015,sprintf('D_{KS} = %1.4f',eval(num2str(KSTest(i-Num_Timepoints)))),'FontSize',9);
                end
                %Paste figures on the subplots:
                copyobj(allchild(get(Subplot(1),'CurrentAxes')),h(1));
                for i=2:Num_Figs
                    copyobj(allchild(get(Subplot(i),'CurrentAxes')),h(i+Num_Timepoints-1));
                end
                legend(h(1),'Measured','Predicted');
            elseif strcmp(Stage,'Test')
                for k=1:Num_Figs
                    saveas(figure(k),['Results/Figures/Test_Step',num2str(Step_Num),'_',num2str(k),'.fig']);
                    Subplot(k) = hgload(['Results/Figures/Test_Step',num2str(Step_Num),'_',num2str(k),'.fig'],struct('Visible','off'));
                end
                KSTrain = KS1;
                for i=1:Num_Timepoints
                    KSTestL(i) = KS2{i,1};
                    KSTestUL(i) = KS3{i,1};
                end
                
                %Figure for labeled testing plots:
                orient(figure,'landscape')
                suptitle('Testing - Labeled')
                h(1) = subplot(6,4,1);
                SubplotsInfo(1) = get(gca);
                title('Training');
                xlabel('Lipid Content (AUC)','FontSize',10)
                ylabel('Probability','FontSize',10)
                SubplotsInfo(1).XAxis.Limits = [0 20]; %Set axis limits for subplots
                SubplotsInfo(1).YAxis.Limits = [0 0.03];
                text(1,0.015,sprintf('D_{KS} = %1.4f',eval(num2str(KSTrain))),'FontSize',9);
                for i=5:21 %Create subplots and set axes and title for each
                    h(i)=subplot(6,4,i);
                    SubplotsInfo(i) = get(gca);
                    SubplotsInfo(i).XAxis.Limits = [0 20]; %Set axis limits for subplots
                    SubplotsInfo(i).YAxis.Limits = [0 0.03];
                    text(1,0.015,sprintf('D_{KS} = %1.4f',eval(num2str(KSTestL(i-4)))),'FontSize',9);
                end
                %Paste figures on the subplots:
                copyobj(allchild(get(Subplot(1),'CurrentAxes')),h(1));
                for i=2:18
                    copyobj(allchild(get(Subplot(i),'CurrentAxes')),h(i+3));
                end
                hL = legend(h(1),'Measured','Predicted');
                newPosition = [0.8 0.9 0.05 0.05];
                newUnits = 'normalized';
                set(hL,'Position', newPosition,'Units', newUnits);
                
                %Figure for unlabeled testing plots:
                orient(figure,'landscape')
                suptitle('Testing - Unlabeled')
                h(1) = subplot(6,4,1);
                SubplotsInfo(1) = get(gca);
                title('Training');
                xlabel('Lipid Content (AUC)','FontSize',10)
                ylabel('Probability','FontSize',10)
                %                 set(gca,'XScale','log'); %Set x-axis to log
                SubplotsInfo(1).XAxis.Limits = [0 20]; %Set axis limits for subplots
                SubplotsInfo(1).YAxis.Limits = [0 0.03];
                text(1,0.015,sprintf('D_{KS} = %1.4f',eval(num2str(KSTrain))),'FontSize',9);
                for i=5:21 %Create subplots and set axes and title for each
                    h(i)=subplot(6,4,i);
                    SubplotsInfo(i) = get(gca);
                    %                     set(gca,'XScale','log'); %Set x-axis to log
                    SubplotsInfo(i).XAxis.Limits = [0 20]; %Set axis limits for subplots
                    SubplotsInfo(i).YAxis.Limits = [0 0.03];
                    text(1,0.015,sprintf('D_{KS} = %1.4f',eval(num2str(KSTestUL(i-4)))),'FontSize',9);
                end
                % Paste figures on the subplots.
                copyobj(allchild(get(Subplot(1),'CurrentAxes')),h(1));
                for i=19:35
                    copyobj(allchild(get(Subplot(i),'CurrentAxes')),h(i-14));
                end
                hL = legend(h(1),'Measured','Predicted');
                newPosition = [0.8 0.9 0.05 0.05];
                newUnits = 'normalized';
                set(hL,'Position', newPosition,'Units', newUnits);
            end
        end
        
        function [h] = FigureGenerator_Comb(Stage,Num_Timepoints,Step_Num,KS1,KS2,KS3) % Generates figures
            if strcmp(Stage,'TrainValid')
                m = 7;
            elseif strcmp(Stage,'Test')
                m = 35;
            end
            for k=1:3
                saveas(figure(k+m),['Results/Figures/Step',num2str(Step_Num),'_',num2str(k+m),'.fig']);
                Subplot(k+m) = hgload(['Results/Figures/Step',num2str(Step_Num),'_',num2str(k+m),'.fig'],struct('Visible','off'));
            end
            KS = [KS1 KS2 KS3];
            %Prepare subplots
            orient(figure,'landscape')
            title_names = {'Training','Testing - Labeled','Testing - Unlabeled'};
            for i=1:3 %Create subplots and set axes and title for each
                h(i+m)=subplot(2,2,i);
                SubplotsInfo(i) = get(gca);
                title(title_names{i},'FontSize',10)
                xlabel('Lipid Content (AUC)','FontSize',10)
                ylabel('Probability','FontSize',10)
                %                 set(gca,'XScale','log'); %Set x-axis to log
                SubplotsInfo(i).XAxis.Limits = [0 20]; %Set axis limits for subplots
                SubplotsInfo(i).YAxis.Limits = [0 0.03];
                text(1,0.015,sprintf('D_{KS} = %1.4f',eval(num2str(KS(i)))),'FontSize',9);
            end
            %Paste figures on the subplots:
            for i=1:3
                copyobj(allchild(get(Subplot(i+m),'CurrentAxes')),h(i+m));
            end
            legend(h(3+m),'Measured','Predicted');
        end
        
    end
end
