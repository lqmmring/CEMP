%APCLUSTER Affinity Propagation Clustering (Frey/Dueck, Science 2007) 
% [idx,netsim,dpsim,expref]=APCLUSTER(s,p) clusters data, using a set  
% of real-valued pairwise data point similarities as input. Clusters  
% are each represented by a cluster center data point (the "exemplar").  
% The method is iterative and searches for clusters so as to maximize  
% an objective function, called net similarity. 
%  
% For N data points, there are potentially N^2-N pairwise similarities;  
% this can be input as an N-by-N matrix 's', where s(i,k) is the  
% similarity of point i to point k (s(i,k) needn?t equal s(k,i)).  In  
% fact, only a smaller number of relevant similarities are needed; if  
% only M similarity values are known (M < N^2-N) they can be input as  
% an M-by-3 matrix with each row being an (i,j,s(i,j)) triple. 
%  
% APCLUSTER automatically determines the number of clusters based on  
% the input preference 'p', a real-valued N-vector. p(i) indicates the  
% preference that data point i be chosen as an exemplar. Often a good  
% choice is to set all preferences to median(s); the number of clusters  
% identified can be adjusted by changing this value accordingly. If 'p'  
% is a scalar, APCLUSTER assumes all preferences are that shared value. 
%  
% The clustering solution is returned in idx. idx(j) is the index of  
% the exemplar for data point j; idx(j)==j indicates data point j  
% is itself an exemplar. The sum of the similarities of the data points to  
% their exemplars is returned as dpsim, the sum of the preferences of  
% the identified exemplars is returned in expref and the net similarity  
% objective function returned is their sum, i.e. netsim=dpsim+expref. 
%  
% 	[ ... ]=apcluster(s,p,'NAME',VALUE,...) allows you to specify  
% 	  optional parameter name/value pairs as follows: 
%  
%   'maxits'     maximum number of iterations (default: 1000) 
%   'convits'    if the estimated exemplars stay fixed for convits  
%          iterations, APCLUSTER terminates early (default: 100) 
%   'dampfact'   update equation damping level in [0.5, 1).  Higher  
%        values correspond to heavy damping, which may be needed  
%        if oscillations occur. (default: 0.9) 
%   'plot'       (no value needed) Plots netsim after each iteration 
%   'details'    (no value needed) Outputs iteration-by-iteration  
%      details (greater memory requirements) 
%   'nonoise'    (no value needed) APCLUSTER adds a small amount of  
%      noise to 's' to prevent degenerate cases; this disables that. 
%  
% Copyright (c) B.J. Frey & D. Dueck (2006). This software may be  
% freely used and distributed for non-commercial purposes. 
%          (RUN APCLUSTER WITHOUT ARGUMENTS FOR DEMO CODE) 
function [idx,netsim,dpsim,expref]=apcluster(s,p,varargin); 
start = clock; 
% Handle arguments to function 
if nargin<2 error('Too few input arguments');  %必须输入至少两个参数，即相似度矩阵和聚类倾向程度
else 
    maxits=1000; convits=100; lam=0.9; plt=0; details=0; nonoise=0;     %设置其它参数
    i=1; 
    while i<=length(varargin)     
        if strcmp(varargin{i},'plot')   %画出，收敛曲线
            plt=1; i=i+1;          
        elseif strcmp(varargin{i},'details') 
            details=1; i=i+1;        
		elseif strcmp(varargin{i},'sparse') 
% 			[idx,netsim,dpsim,expref]=apcluster_sparse(s,p,varargin{:}); 
			fprintf('''sparse'' argument no longer supported; see website for additional software\n\n'); 
			return; 
        elseif strcmp(varargin{i},'nonoise') 
            nonoise=1; i=i+1; 
        elseif strcmp(varargin{i},'maxits')    %最大迭代次数
            maxits=varargin{i+1};   
            i=i+2; 
            if maxits<=0 error('maxits must be a positive integer'); end; 
        elseif strcmp(varargin{i},'convits')   %收敛迭代间隔
            convits=varargin{i+1}; 
            i=i+2; 
            if convits<=0 error('convits must be a positive integer'); end; 
        elseif strcmp(varargin{i},'dampfact')  %阻尼系数，防止数据震荡
            lam=varargin{i+1}; 
            i=i+2; 
            if (lam<0.5)||(lam>=1) 
                error('dampfact must be >= 0.5 and < 1'); 
            end; 
        else i=i+1; 
        end; 
    end; 
end; 
if lam>0.9   %阻尼系数最好不要大于0.9
    fprintf('\n*** Warning: Large damping factor in use. Turn on plotting\n'); 
    fprintf('    to monitor the net similarity. The algorithm will\n'); 
    fprintf('    change decisions slowly, so consider using a larger value\n'); 
    fprintf('    of convits.\n\n'); 
end; 
 
% Check that standard arguments are consistent in size
%检查相似矩阵，需要二维数据，参数P必须是个标量或矢量
if length(size(s))~=2 error('s should be a 2D matrix'); 
elseif length(size(p))>2 error('p should be a vector or a scalar');   %
elseif size(s,2)==3 %检查，为形成相似矩阵的数据，是否完整
    tmp=max(max(s(:,1)),max(s(:,2))); 
    if length(p)==1 N=tmp; else N=length(p); end; 
    if tmp>N 
        error('data point index exceeds number of data points'); 
    elseif min(min(s(:,1)),min(s(:,2)))<=0 
        error('data point indices must be >= 1'); 
    end; 
elseif size(s,1)==size(s,2) 
    N=size(s,1); 
    if (length(p)~=N)&&(length(p)~=1) 
        error('p should be scalar or a vector of size N'); 
    end; 
else error('s must have 3 columns or be square'); end; 
 
% Construct similarity matrix   构建相似矩阵
if N>3000 
    fprintf('\n*** Warning: Large memory request. Consider activating\n'); 
    fprintf('    the sparse version of APCLUSTER.\n\n'); 
end; 
if size(s,2)==3 && size(s,1)~=3,    %针对的是未形成相似矩阵的数据
    S=-Inf*ones(N,N,class(s));  
    for j=1:size(s,1), S(s(j,1),s(j,2))=s(j,3); end; 
else S=s; 
end; 
 
if S==S', symmetric=true; else symmetric=false; end; 
realmin_=realmin(class(s)); realmax_=realmax(class(s)); 
 
% In case user did not remove degeneracies from the input similarities, 
% avoid degenerate solutions by adding a small amount of noise to the 
% input similarities       
if ~nonoise 
    rns=randn('state'); randn('state',0); 
    S=S+(eps*S+realmin_*100).*rand(N,N); 
    randn('state',rns); 
end; 
 
% Place preferences on the diagonal of S   把聚类倾向系数填入相似矩阵的对角线
if length(p)==1 for i=1:N S(i,i)=p; end; 
else for i=1:N S(i,i)=p(i); end; 
end; 
 
% Numerical stability -- replace -INF with -realmax   用最小实数代替非数字
n=find(S<-realmax_); if ~isempty(n), warning('-INF similarities detected; changing to -REALMAX to ensure numerical stability'); S(n)=-realmax_; end; clear('n'); 
if ~isempty(find(S>realmax_,1)), error('+INF similarities detected; change to a large positive value (but smaller than +REALMAX)'); end; 
 
 
% Allocate space for messages, etc   给结果分配空间
dS=diag(S); A=zeros(N,N,class(s)); R=zeros(N,N,class(s)); t=1; 
if plt, netsim=zeros(1,maxits+1); end; 
if details 
    idx=zeros(N,maxits+1); 
    netsim=zeros(1,maxits+1);  
    dpsim=zeros(1,maxits+1);  
    expref=zeros(1,maxits+1);  
end; 
 
% Execute parallel affinity propagation updates 
e=zeros(N,convits); dn=0; i=0; 
if symmetric, ST=S; else ST=S'; end; % saves memory if it's symmetric 
while ~dn 
    i=i+1;  
 
    % Compute responsibilities 
	%%A=A'; R=R';
    T=A+R;
    A=A';  
    R=R'; 
	for ii=1:N, 
		old = R(:,ii); 
		AS = A(:,ii) + ST(:,ii); [Y,I]=max(AS); AS(I)=-Inf; 
		[Y2,I2]=max(AS); 
		R(:,ii)=ST(:,ii)-Y; 
		R(I,ii)=ST(I,ii)-Y2; 
		R(:,ii)=(1-lam)*R(:,ii)+lam*old; % Damping 
        R(R(:,ii)>realmax_,ii)=realmax_; 
	end; 
	A=A'; R=R'; 
 
    % Compute availabilities 
	for jj=1:N, 
		old = A(:,jj); 
		Rp = max(R(:,jj),0); Rp(jj)=R(jj,jj); 
		A(:,jj) = sum(Rp)-Rp; 
		dA = A(jj,jj); A(:,jj) = min(A(:,jj),0); A(jj,jj) = dA; 
		A(:,jj) = (1-lam)*A(:,jj) + lam*old; % Damping 
	end; 
	 
    % Check for convergence 
    E=((diag(A)+diag(R))>0); e(:,mod(i-1,convits)+1)=E; K=sum(E); 
    if i>=convits || i>=maxits, 
        se=sum(e,2); 
        unconverged=(sum((se==convits)+(se==0))~=N); 
        if (~unconverged&&(K>0))||(i==maxits) dn=1; end; 
    end; 
 
    % Handle plotting and storage of details, if requested 
    if plt||details 
        if K==0 
            tmpnetsim=nan; tmpdpsim=nan; tmpexpref=nan; tmpidx=nan; 
        else 
            I=find(E); notI=find(~E); [tmp c]=max(S(:,I),[],2); c(I)=1:K; tmpidx=I(c); 
            tmpdpsim=sum(S(sub2ind([N N],notI,tmpidx(notI)))); 
            tmpexpref=sum(dS(I)); 
            tmpnetsim=tmpdpsim+tmpexpref; 
        end; 
    end; 
    if details 
        netsim(i)=tmpnetsim; dpsim(i)=tmpdpsim; expref(i)=tmpexpref; 
        idx(:,i)=tmpidx; 
    end; 
    if plt, 
        netsim(i)=tmpnetsim; 
		figure(234); 
        plot(((netsim(1:i)/10)*100)/10,'r-'); xlim([0 i]); % plot barely-finite stuff as infinite 
        xlabel('# Iterations'); 
        ylabel('Fitness (net similarity) of quantized intermediate solution'); 
%         drawnow;  
    end; 
end; % iterations 
I=find((diag(A)+diag(R))>0); K=length(I); % Identify exemplars 
if K>0 
    [tmp c]=max(S(:,I),[],2); c(I)=1:K; % Identify clusters 
    % Refine the final set of exemplars and clusters and return results 
    for k=1:K ii=find(c==k); [y j]=max(sum(S(ii,ii),1)); I(k)=ii(j(1)); end; notI=reshape(setdiff(1:N,I),[],1); 
    [tmp c]=max(S(:,I),[],2); c(I)=1:K; tmpidx=I(c); 
	tmpdpsim=sum(S(sub2ind([N N],notI,tmpidx(notI)))); 
	tmpexpref=sum(dS(I)); 
	tmpnetsim=tmpdpsim+tmpexpref; 
else 
    tmpidx=nan*ones(N,1); tmpnetsim=nan; tmpexpref=nan; 
end; 
if details 
    netsim(i+1)=tmpnetsim; netsim=netsim(1:i+1); 
    dpsim(i+1)=tmpdpsim; dpsim=dpsim(1:i+1); 
    expref(i+1)=tmpexpref; expref=expref(1:i+1); 
    idx(:,i+1)=tmpidx; idx=idx(:,1:i+1); 
else 
    netsim=tmpnetsim; dpsim=tmpdpsim; expref=tmpexpref; idx=tmpidx; 
end; 
if plt||details 
    fprintf('\nNumber of exemplars identified: %d  (for %d data points)\n',K,N); 
    fprintf('Net similarity: %g\n',tmpnetsim); 
    fprintf('  Similarities of data points to exemplars: %g\n',dpsim(end)); 
    fprintf('  Preferences of selected exemplars: %g\n',tmpexpref); 
    fprintf('Number of iterations: %d\n\n',i); 
	fprintf('Elapsed time: %g sec\n',etime(clock,start)); 
end; 
if unconverged 
	fprintf('\n*** Warning: Algorithm did not converge. Activate plotting\n'); 
	fprintf('    so that you can monitor the net similarity. Consider\n'); 
	fprintf('    increasing maxits and convits, and, if oscillations occur\n'); 
	fprintf('    also increasing dampfact.\n\n'); 
end;
