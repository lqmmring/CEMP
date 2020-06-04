function d=persondist(a,b)
% SQDIST - computes person colleration distance matrix
%          computes a rectangular matrix of pairwise distances
% between points in A (given in columns) and points in B

% NB: very fast implementation taken from Roland Bunschoten

fenzi = sum(a.*b)-(sum(a)*sum(b))/length(a);
fenmu = sqrt((sum(a.*a)-sum(a)^2/length(a))*(sum(b.*b)-sum(b)^2/length(a))); 
d = fenzi/fenmu;
end