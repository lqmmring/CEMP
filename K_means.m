function [idx, center, m] = K_means(data, m, MaxIter, metrices)
% function [idx, center, m] = eff_kmeans(data, m, MaxIter);
if ~exist('metrices', 'var') || isempty(metrices)
        metrices = 'sqdist';
end
[n, dim] = size(data);
dex = randperm(n);
center = data(dex(1:m),:);
similarity = zeros(n,n);
for i = 1:MaxIter;
    nul = zeros(m,1);
    for c =1:m
        aa=center(c,:)';
        for d = 1:n
            bb=data(d,:)';
            eval(['similarity(c,d)=',metrices,'(aa, bb);'])
        end
    end
    [xx, idx]=min(similarity);
    for j = 1:m;
        dex = find(idx == j);
        l = length(dex);
        cltr = data(dex,:);
%         center(j,:) = mean(cltr);
        if l > 1;
            center(j,:) = mean(cltr);
        elseif l == 1;
            center(j,:) = cltr;
        else
            nul(j) = 1;
        end;
    end;
    dex = find(nul == 0);
    m = length(dex);
    center = center(dex,:);
end;

end