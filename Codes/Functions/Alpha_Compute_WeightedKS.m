function KS = Alpha_Compute_WeightedKS(a,XUL,YL,M)
%Author: Mohammad Tanhaemami
%Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.

% X & Y should be in log scale
a = exp(a);
WM = a(1)*M{1,1} + a(2)*M{2,1} + a(3)*M{3,1};
YPred = XUL*WM;
[~,~,KS] = kstest2(YPred,YL);
end