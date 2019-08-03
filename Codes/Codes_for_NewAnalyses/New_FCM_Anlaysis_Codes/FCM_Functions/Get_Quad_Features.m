function [XQ] = Get_Quad_Features(X)
%This function takes values from any given X matrix and generates a new
%matrix with quadratic values. The resulting matrix XQ can be used for
%performing linear regression (simple or other) analyses.
N_data = size(X,1);      %Record number of datapoints.
N_features = size(X,2);  %Record number of features.
XQ = zeros(N_data,N_features+N_features*(N_features+1)/2); %Initiate quadratic matrix XQ.
XQ(:,1:N_features) = X; %Values  for XQ are the same as X for columns 1 to N_features.
k = N_features; %Start counter from N_features to move forward to the final column of XQ.
for i=1:N_features
    for j = i:N_features
        k=k+1; %increase counter by one.
        %Mupltiply each element to generate quadaratic elements and assign
        %them to the remaining columns of XQ until the final column:
        XQ(:,k) = X(:,i).*X(:,j);
    end
end