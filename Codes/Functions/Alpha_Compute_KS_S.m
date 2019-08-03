function SKS = Alpha_Compute_KS_S(SS,QS,XULS,YLS,MS,IS)
%Author: Mohammad Tanhaemami
%Organization: Dr. Brian Munsky Research Group, Department of Chemical and Biological Engineering, Colorado State University.

aS = SS(:,IS==1)*QS(IS==1,:);
YPredS = cell(0);
for i = 1:size(XULS,1)
    WMS(:,i) = aS(i,1)*MS{1,1} + aS(i,2)*MS{2,1} + aS(i,3)*MS{3,1};
    YPredS{i,1} = XULS{i,1}*WMS(:,i);
    [~,~,KSS(i)] = kstest2(YPredS{i,1},YLS{i,1});
end
SKS = mean(KSS)+ 0.001*sum(IS);
end