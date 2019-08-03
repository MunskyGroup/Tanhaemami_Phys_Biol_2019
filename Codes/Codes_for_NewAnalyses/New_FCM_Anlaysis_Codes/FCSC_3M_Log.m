classdef FCSC_3M_Log
    %Author: Mohammad Tanhaemami
    %Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.
    % This object is used for our final weighted modeling multi-stage machine learning strategy
    properties % Input values as follows:
        Step % Step of the calculation according to the main script
        Stage % Stage of computations (training or testing) used for generating figures
        BINS % Includes information regarding the bins used for the histograms
        XTrainL % Labeled training features
        XTrainUL % Unlabeled training features
        XTrainUL_S % Unlabeled training features inclding information from main channels (to be used for secondary regression analyses) 
        YTrainL % Labeled training targets
        XTestL % Labeled validation (or testing) features
        XTestUL % Unlabeled validation (or testing) features
        XTestUL_S % Unlabeled validation (or testing) features inclding information from main channels (to be used for secondary regression analyses) 
        YTestL % Labeled validation (or testing) targets
        I % Vector indicating which features to be used as a result of applying the genetic algorithm
        M % Vector containing regression coefficients
        Q % The weight quotient
        IforS % Vector indicating which columns of the statistics matrix (S) or rows of the weight quotient (Q) to be used as a result of applying the genetic algorithm
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
            for i = 1:size(obj.XTrainL,1)
                Model.LogXTrainL{i,1} = log(obj.XTrainL{i,1}); % Labeled training features
                Model.LogXTrainUL{i,1} = log(obj.XTrainUL{i,1}); % Unlabeled training features
                Model.LogXTrainUL_S{i,1} = log(obj.XTrainUL_S{i,1}); % Unlabeled training data inclding information from main channels (to be used for secondary regression analyses)
                Model.LogYTrainL{i,1} = log(obj.YTrainL{i,1}); % Labled training targets
            end
            for i = 1:size(obj.XTestL,1)
                Model.LogXTestL{i,1} = log(obj.XTestL{i,1}); % Labled validation (or testing) features
                Model.LogXTestUL{i,1} = log(obj.XTestUL{i,1}); % Unlabled validation (or testing) features
                Model.LogXTestUL_S{i,1} = log(obj.XTestUL_S{i,1}); % Unlabeled validation (or testing) data inclding information from main channels (to be used for secondary regression analyses)
                Model.LogYTestL{i,1} = log(obj.YTestL{i,1}); % Labled validation (or testing) targets
            end
            % Predict training timepoints (separated):
            for i = 1:size(obj.XTrainL,1)
                % Unlabeled Cells:
                Model.STrUL(i,:) = [mean(Model.LogXTrainUL_S{i,1}),sqrt(var(Model.LogXTrainUL_S{i,1}))]; % Create the test statistics matrix (S)
                Model.AlphaTrUL(i,:) = Model.STrUL(i,obj.IforS==1)*obj.Q(obj.IforS==1,:); % Compute the weight as alpha = [statistics matrix (S)]*[weight quotient (Q)]
                Model.MwTrUL(:,i) = Model.AlphaTrUL(i,1)*obj.M{1,1} + Model.AlphaTrUL(i,2)*obj.M{2,1} + Model.AlphaTrUL(i,3)*obj.M{3,1}; % Compute the weighted model by using the weights nad regression coefficients weighte_model = alpha1*M1 + alpha2*M2 + alpha3*M3
                Model.LogYPredTrainUL{i,1} = Model.LogXTrainUL{i,1}*Model.MwTrUL(:,i); % Predict the unlabeled target values
                [~,~,Model.ErrTrainUL(i)] = kstest2(Model.LogYPredTrainUL{i,1},Model.LogYTrainL{i,1}); % Compute the KS distance between predicted and the measure labeled targets
                Model.LogYPredTrainL{i,1} = Model.LogXTrainL{i,1}*Model.MwTrUL(:,i); % Predict the labeled target values (we predict the labeled cells by applying the same weights as we used for unlabeld cells).
                [~,~,Model.ErrTrainL(i)] = kstest2(Model.LogYPredTrainL{i,1},Model.LogYTrainL{i,1}); % Compute the KS distance between predicted and the measure labeled targets
            end
            % Predict validation (or testing) timepoints (separated):
            for i = 1:size(obj.XTestL,1)
                % Unlabeled Cells:
                Model.STsUL(i,:) = [mean(Model.LogXTestUL_S{i,1}),sqrt(var(Model.LogXTestUL_S{i,1}))]; % Test statistics matrix (S)
                Model.AlphaTsUL(i,:) = Model.STsUL(i,obj.IforS==1)*obj.Q(obj.IforS==1,:); % Weights
                Model.MwTsUL(:,i) = Model.AlphaTsUL(i,1)*obj.M{1,1} + Model.AlphaTsUL(i,2)*obj.M{2,1} + Model.AlphaTsUL(i,3)*obj.M{3,1}; % weighted model
                Model.LogYPredTestUL{i,1} = Model.LogXTestUL{i,1}*Model.MwTsUL(:,i); % Predicted unlabeled targets
                [~,~,Model.ErrTestUL(i)] = kstest2(Model.LogYPredTestUL{i,1},Model.LogYTestL{i,1}); % KS distance
                % Labeled Cells:
                Model.LogYPredTestL{i,1} = Model.LogXTestL{i,1}*Model.MwTsUL(:,i); % Predict the labeled targets
                [~,~,Model.ErrTestL(i)] = kstest2(Model.LogYPredTestL{i,1},Model.LogYTestL{i,1}); % KS distance
            end
            
            %Histograms:
            %Histograms on top of each other for each step:
            if strcmp(obj.MakeHistPlot,'True')
                Model.HistTrainL = FCSC_3M_Log.Hist_Plot('Separated',Model.LogYTrainL,Model.LogYPredTrainL,obj.BINS);
                title('HistTrainL');
                Model.HistTrainUL = FCSC_3M_Log.Hist_Plot('Separated',Model.LogYTrainL,Model.LogYPredTrainUL,obj.BINS);
                title('HistTrainUL');
                Model.HistTestL = FCSC_3M_Log.Hist_Plot('Separated',Model.LogYTestL,Model.LogYPredTestL,obj.BINS);
                title('HistValid(Test)L');
                Model.HistTestUL = FCSC_3M_Log.Hist_Plot('Separated',Model.LogYTestL,Model.LogYPredTestUL,obj.BINS);
                title('HistValid(Test)UL');
                if strcmp(obj.GenerateFigure,'True')
                    Model.FigureSep = FCSC_3M_Log.FigureGenerator_Sep(obj.Stage,size(obj.XTestL,1),obj.Step,Model.ErrTrainL,Model.ErrTrainUL,Model.ErrTestL,Model.ErrTestUL);
                end
            end
        end
    end
    %%
    methods(Static)
        
        function WKS = Compute_WeightedKS_S(S,Q,LogXTestUL,LogYTestL,M,I) % This is the objfunc for GA on S
            a = S(:,I==1)*Q(I==1,:); % [weights (alpha)] = [statistics matrix (S)] * [weight quotient (Q)]
            WM = a(1)*M{1,1} + a(2)*M{2,1} + a(3)*M{3,1}; % [weighted model (WM)] = alpha1*M1 + alpha2*M2 + alpha3*M3
            LogYPred = LogXTestUL*WM; % Predict with regression
            [~,~,WKS] = kstest2(LogYPred,LogYTestL); % KS distance between predicted and measured
        end
        
        function KSDist = Compute_KS(XTrain,YTrain,XTest,YTest,I) % Calculates the KS distance b/w measurements and predictions for the GA
            LogXTrain = log(XTrain); % Get the log values
            LogYTrain = log(YTrain);
            LogXTest = log(XTest);
            LogYTest = log(YTest);
            M = LogXTrain(:,I==1)\LogYTrain; % Compute the regression coefficient (the model)
            LogYPred = LogXTest(:,I==1)*M; % Predict the targets
            [~,~,KSDist] = kstest2(LogYPred,LogYTest); % Compute the KS distsnce
        end
        
        function Make_Hist_Plot = Hist_Plot(TimepointsFormat,Y,YPred,BINS) % Function to plot histograms of the results
            if strcmp(TimepointsFormat,'Separated')
                for i = 1:size(Y,1)
                    figure; clf;
                    A{i,1} = hist(Y{i,1},BINS)/sum(hist(Y{i,1},BINS)); % Take histogram values of measured target
                    B{i,1} = hist(YPred{i,1},BINS)/sum(hist(YPred{i,1},BINS)); % Take histogram values of predicted values
                    Make_Hist_Plot{i,1} = stairs(BINS',[A{i,1};B{i,1}]','LineWidth',0.5); % Plot histograms
                    set(gca,'fontsize',16)
                    legend('Measured 1','Predicted 1','Measured 2','Predicted 2','Measured 3','Predicted 3');
                end
            elseif strcmp(TimepointsFormat,'Combined')
                figure; clf;
                A = hist(Y,BINS)/sum(hist(Y,BINS)); %Take histogram values of measured target
                B = hist(YPred,BINS)/sum(hist(YPred,BINS)); %Take histogram values of predicted values
                Make_Hist_Plot = stairs(BINS',[A;B]','LineWidth',0.5); %Plot hist
                set(gca,'fontsize',16)
                legend('Measured','Predicted')
            end
        end
        
        function [h] = FigureGenerator_Sep(Stage,Num_Timepoints,Step_Num,KS1,KS2,KS3,KS4) % Function to generate figures of the results
            if strcmp(Stage,'TrainValid')
                Num_Figs = 4*Num_Timepoints; % Determines number of figures to be generated
                for k = 1:Num_Figs
                    saveas(figure(k),['Results/Figures/Step',num2str(Step_Num),'_',num2str(k),'.fig']);
                    Subplot(k) = hgload(['Results/Figures/Step',num2str(Step_Num),'_',num2str(k),'.fig'],struct('Visible','off')); % Get the figure information
                end
                KSTrainTest = [KS1 KS2 KS3 KS4];
                %Prepare subplots:
                orient(figure,'landscape')
                titleTrainL = repmat({'Training - Labeled'},1,Num_Timepoints);
                titleTrainUL = repmat({'Training - Unlabeled'},1,Num_Timepoints);
                titleTestL = repmat({'Validation - Labeled'},1,Num_Timepoints);
                titleTestUL = repmat({'Validation - Unlabeled'},1,Num_Timepoints);
                title_names = [titleTrainL,titleTrainUL,titleTestL,titleTestUL];
                for i = 1:4*Num_Timepoints %Create subplots and set axes and title for each
                    h(i) = subplot(4,Num_Timepoints,i);
                    SubplotsInfo(i) = get(gca); % Get axis information
                    title(title_names{1,i},'FontSize',10); %Title
                    xlabel('Lipid Content (AUC)','FontSize',10); % X-axis label
                    ylabel('Probability','FontSize',10); % Y-axis label
                    SubplotsInfo(i).XAxis.Limits = [5 25]; % Set axis limits for subplots
                    SubplotsInfo(i).YAxis.Limits = [0 0.08];
                    text(6,0.04,sprintf('D_{KS} = %1.4f',eval(num2str(KSTrainTest(i)))),'FontSize',9); % Display the KS distance on plot
                end
                % Paste figures on the subplots:
                for i = 1:Num_Figs
                    copyobj(allchild(get(Subplot(i),'CurrentAxes')),h(i));
                end
                legend(h(1),'Measured','Predicted');
            elseif strcmp(Stage,'Test')
                Num_Figs = 2*Num_Timepoints+6;
                for k=1:Num_Figs
                    saveas(figure(k),['Results/Figures/Test_Step',num2str(Step_Num),'_',num2str(k),'.fig']);
                    Subplot(k) = hgload(['Results/Figures/Test_Step',num2str(Step_Num),'_',num2str(k),'.fig'],struct('Visible','off'));
                end
                KSTestL = KS3;
                KSTestUL = KS4;
                
                % Figure for labeled testing plots:
                orient(figure,'landscape')
                suptitle('Testing - Labeled')
                for i = 1:17 % Create subplots and set axes and title for each
                    h(i) = subplot(5,4,i);
                    SubplotsInfo(i) = get(gca);
                    xlabel('Lipid Content (AUC)','FontSize',10);
                    ylabel('Probability','FontSize',10);
                    SubplotsInfo(i).XAxis.Limits = [5 25]; %Set axis limits for subplots
                    SubplotsInfo(i).YAxis.Limits = [0 0.08];
                    text(6,0.04,sprintf('D_{KS} = %1.4f',eval(num2str(KSTestL(i)))),'FontSize',9);
                end
                % Paste figures on the subplots:
                for i=1:17
                    copyobj(allchild(get(Subplot(i+6),'CurrentAxes')),h(i));
                end
                % Set the legend information and position:
                hL = legend(h(4),'Measured','Predicted');
                newPosition = [0.8 0.9 0.05 0.05];
                newUnits = 'normalized';
                set(hL,'Position', newPosition,'Units', newUnits);
                
                % Figure for unlabeled testing plots:
                orient(figure,'landscape');
                suptitle('Testing - Unlabeled');
                for i=1:17 %Create subplots and set axes and title for each
                    h(i)=subplot(5,4,i);
                    SubplotsInfo(i) = get(gca);
                    xlabel('Lipid Content (AUC)','FontSize',10);
                    ylabel('Probability','FontSize',10);
                    SubplotsInfo(i).XAxis.Limits = [5 25]; %Set axis limits for subplots
                    SubplotsInfo(i).YAxis.Limits = [0 0.08];
                    text(2,0.02,sprintf('D_{KS} = %1.4f',eval(num2str(KSTestUL(i)))),'FontSize',9);
                end
                % Paste figures on the subplots.
                for i=1:17
                    copyobj(allchild(get(Subplot(i+23),'CurrentAxes')),h(i));
                end
                % Set the legend information and position:
                hL = legend(h(4),'Measured','Predicted');
                newPosition = [0.8 0.9 0.05 0.05];
                newUnits = 'normalized';
                set(hL,'Position', newPosition,'Units', newUnits);
            end
        end
    end
end