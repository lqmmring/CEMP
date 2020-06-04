function RC = Cal_RC(X, labels)
%==========================================================================
% FUNCTION: RC = Cal_RC(X, label) 
% DESCRIPTION: A function for computing ratio cut 
%              for a clustering result
%
% INPUTS:  X = a dataset, rows of X correspond to observations; columns
%              correspond to variables (exclude class labels!!)
%     labels = cluster labels from a clustering result (N-by-1 vector)
%
% OUTPUT: RC = ratio cut score
%==========================================================================
% copyright (c) 2019 Liu & Wang
%==========================================================================

C = unique(labels); %all clusters
k = length(C); %number of clusters
RC = 0; %initialize compactness

for i=1:k %for each cluster
    ind1 = find(labels==C(i)); %find data point members for the i-th cluster
    ind2 = find(labels~=C(i)); %find data point members not for the i-th cluster
    nk = length(ind1);
    if nk <= 1 %singleton cluster
        RC = RC + 0;
    else
        sum_d = sum(sum(pdist2(X(ind1,:),X(ind2,:))))/2;
        RC = RC + (sum_d/nk);
        
    end
end
% fprintf('RCָ�꣺%f   \n',RC);
end