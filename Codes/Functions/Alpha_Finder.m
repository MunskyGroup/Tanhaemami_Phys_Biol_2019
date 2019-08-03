function Alphas = Alpha_Finder(DATA,I,FeatureInds,FeatureInds_S,FeaturesType,PrevBestVecS)
%Author: Mohammad Tanhaemami
%Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.

% This function performs analysis to develop a weighted model by performing a secondary regression analysis
% This function performs the following:
%       Computes the models with respect to 3 training time points (M1, M2, and M3)
%       Computes the weight quotient (Q) by taking random samples from the training and validation data
%       Find the best columns of the test staticsics matrix (S), or equivaletly the best rows of the Q via using the GA
Alphas.NCellsL = 6000;
Alphas.NCellsUL = 4500;
TargetInd = 3;
TrainSet = [2,10,23];
ValidSet = [1,11,22];
if strcmp(FeaturesType,'Linear')
    Alphas.FeatursType = 'Linear Features';
elseif strcmp(FeaturesType,'Quadratic')
    Alphas.FeatursType = 'Quadratic Features';
end
[XTrainL,XTrainUL,XTrainUL_S,YTrainL,XValidUL,XValidUL_S,YValidL] = Alpha_CollectData(Alphas.NCellsL,Alphas.NCellsUL,FeatureInds,FeatureInds_S,FeaturesType,TargetInd,TrainSet,ValidSet,DATA);
for i = 1:size(XTrainL,1)
    Alphas.LogXTrainL{i,1} = log(XTrainL{i,1});
    Alphas.LogXTrainUL{i,1} = log(XTrainUL{i,1});
    Alphas.LogXTrainUL_S{i,1} = log(XTrainUL_S{i,1});
    Alphas.LogYTrainL{i,1} = log(YTrainL{i,1});
    Alphas.LogXValidUL{i,1} = log(XValidUL{i,1});
    Alphas.LogXValidUL_S{i,1} = log(XValidUL_S{i,1});
    Alphas.LogYValidL{i,1} = log(YValidL{i,1});
end

%% Compute models for the 3 training timepoints
Alphas.M = cell(3,1);
for i = 1:size(Alphas.LogXTrainL,1)
    Alphas.M{i,1} = zeros(size(Alphas.LogXTrainL{1,1},2),1);
    Alphas.M{i,1}(I(:,i)==1) = Alphas.LogXTrainL{i,1}(:,I(:,i)==1)\Alphas.LogYTrainL{i,1};
end

%% Sample from the data and compute alpha, Q, and S
% M = alpha1*M1 + alpha2*M2 + alpha3*M3
% S: The test statistic defined as S = [mean , covariance]
% Q: The quantity that maps alpha to S as: alpha = S * Q ------> Q = S\alpha
options = optimset('Display','off');
nrep = 1;
for i = 1:nrep*6
    if i <= nrep
        rng(i);
        R1 = ceil(6000*rand(3000,1));
        R2 = ceil(6000*rand(3000,1));
        XUL = Alphas.LogXTrainUL{1,1}(R1,:);
        XUL_S = Alphas.LogXTrainUL_S{1,1}(R1,:);
        YL = Alphas.LogYTrainL{1,1}(R2,:);
    elseif i <= 2*nrep
        rng(i);
        R1 = ceil(6000*rand(3000,1));
        R2 = ceil(6000*rand(3000,1));
        XUL = Alphas.LogXTrainUL{2,1}(R1,:);
        XUL_S = Alphas.LogXTrainUL_S{2,1}(R1,:);
        YL = Alphas.LogYTrainL{2,1}(R2,:);
    elseif i <= 3*nrep
        rng(i);
        R1 = ceil(6000*rand(3000,1));
        R2 = ceil(6000*rand(3000,1));
        XUL = Alphas.LogXTrainUL{3,1}(R1,:);
        XUL_S = Alphas.LogXTrainUL_S{3,1}(R1,:);
        YL = Alphas.LogYTrainL{3,1}(R2,:);
    elseif i <= 4*nrep
        rng(i);
        R1 = ceil(6000*rand(3000,1));
        R2 = ceil(6000*rand(3000,1));
        XUL = Alphas.LogXValidUL{1,1}(R1,:);
        XUL_S = Alphas.LogXValidUL_S{1,1}(R1,:);
        YL = Alphas.LogYValidL{1,1}(R2,:);
    elseif i <= 5*nrep
        rng(i);
        R1 = ceil(6000*rand(3000,1));
        R2 = ceil(6000*rand(3000,1));
        XUL = Alphas.LogXValidUL{2,1}(R1,:);
        XUL_S = Alphas.LogXValidUL_S{2,1}(R1,:);
        YL = Alphas.LogYValidL{2,1}(R2,:);
    elseif i <= 6*nrep
        rng(i);
        R1 = ceil(6000*rand(3000,1));
        R2 = ceil(6000*rand(3000,1));
        XUL = Alphas.LogXValidUL{3,1}(R1,:);
        XUL_S = Alphas.LogXValidUL_S{3,1}(R1,:);
        YL = Alphas.LogYValidL{3,1}(R2,:);
    end
    old_fval = inf;
    Objfun = @(a)Alpha_Compute_WeightedKS(a,XUL,YL,Alphas.M);
    for k = 1:100
        if k == 1
            a0 = [0 -9 -9];
        elseif k == 2
            a0 = [-9 0 -9];
        elseif k == 3
            a0 = [-9 -9 0];
        elseif k == 4
            a0 = [-9 0 0];
        elseif k == 5
            a0 = [0 -9 0];
        elseif k == 6
            a0 = [0 0 -9];
        else
            a0 = randn(1,3)-2;
        end
        [a,fval] = fminsearch(Objfun,a0,options);
        if fval < old_fval
            Alphas.alpha(i,:) = exp(a);
            Alphas.alphaError(i,1) = fval;
            old_fval = fval;
        end
    end
    Alphas.S(i,:) = [mean(XUL_S),sqrt(var(XUL_S))];
end
Alphas.Q = Alphas.S\Alphas.alpha;

%% Which columns of S (or rows of Q) should be selected:
nrep = 1;
rng('default');
for i = 1:nrep*6
    if i <= nrep
        rng(6+i);
        R1 = 6000 + ceil(3000*rand(3000,1));
        R2 = 6000 + ceil(3000*rand(3000,1));
        Alphas.XUL_is{i,1} = Alphas.LogXTrainUL{1,1}(R1,:);
        Alphas.XUL_S_is{i,1} = Alphas.LogXTrainUL_S{1,1}(R1,:);
        Alphas.YL_is{i,1} = Alphas.LogYTrainL{1,1}(R2,:);
    elseif i <= 2*nrep
        rng(6+i);
        R1 = 6000 + ceil(3000*rand(3000,1));
        R2 = 6000 + ceil(3000*rand(3000,1));
        Alphas.XUL_is{i,1} = Alphas.LogXTrainUL{2,1}(R1,:);
        Alphas.XUL_S_is{i,1} = Alphas.LogXTrainUL_S{2,1}(R1,:);
        Alphas.YL_is{i,1} = Alphas.LogYTrainL{2,1}(R2,:);
    elseif i <= 3*nrep
        rng(6+i);
        R1 = 6000 + ceil(3000*rand(3000,1));
        R2 = 6000 + ceil(3000*rand(3000,1));
        Alphas.XUL_is{i,1} = Alphas.LogXTrainUL{3,1}(R1,:);
        Alphas.XUL_S_is{i,1} = Alphas.LogXTrainUL_S{3,1}(R1,:);
        Alphas.YL_is{i,1} = Alphas.LogYTrainL{3,1}(R2,:);
    elseif i <= 4*nrep
        rng(6+i);
        R1 = 6000 + ceil(3000*rand(3000,1));
        R2 = 6000 + ceil(3000*rand(3000,1));
        Alphas.XUL_is{i,1} = Alphas.LogXValidUL{1,1}(R1,:);
        Alphas.XUL_S_is{i,1} = Alphas.LogXValidUL_S{1,1}(R1,:);
        Alphas.YL_is{i,1} = Alphas.LogYValidL{1,1}(R2,:);
    elseif i <= 5*nrep
        rng(6+i);
        R1 = 6000 + ceil(3000*rand(3000,1));
        R2 = 6000 + ceil(3000*rand(3000,1));
        Alphas.XUL_is{i,1} = Alphas.LogXValidUL{2,1}(R1,:);
        Alphas.XUL_S_is{i,1} = Alphas.LogXValidUL_S{2,1}(R1,:);
        Alphas.YL_is{i,1} = Alphas.LogYValidL{2,1}(R2,:);
    elseif i <= 6*nrep
        rng(6+i);
        R1 = 6000 + ceil(3000*rand(3000,1));
        R2 = 6000 + ceil(3000*rand(3000,1));
        Alphas.XUL_S_is{i,1} = Alphas.LogXValidUL_S{3,1}(R1,:);
        Alphas.YL_is{i,1} = Alphas.LogYValidL{3,1}(R2,:);
    end
    Alphas.S_is(i,:) = [mean(Alphas.XUL_S_is{i,1}),sqrt(var(Alphas.XUL_S_is{i,1}))];
end

% Apply GA to select the best columns of S:
objfS = @(I)Alpha_Compute_KS_S(Alphas.S_is,Alphas.Q,Alphas.XUL_is,Alphas.YL_is,Alphas.M,I);
NColumns = size(Alphas.S_is,2);

if isempty(PrevBestVecS) == 1
    Alphas.initpop = 1*(rand(50*NColumns,NColumns)<0.1);
else
    Alphas.initpop = repmat(PrevBestVecS,50*NColumns,1);
end

options = optimoptions('ga','PopulationType','bitstring','UseParallel',true,'Display','off',...
    'InitialPopulationMatrix',Alphas.initpop,'PopulationSize',50*NColumns,'MutationFcn',{@mutationuniform, 0.5/NColumns});
Alphas.IvecBest_is = ga(objfS,NColumns,options);

end
