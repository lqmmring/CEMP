function [Offspring, Offspring_Vel] = Operator(c1,c2,CA,W,Particle,Parent,Pbest,Gbest,Parent_Vel)
% function [Offspring, Offspring_Vel] = Operator(lower,upper,W,Particle,Parent,Pbest,Gbest,Parent_Vel,Parameter)
% <operator> <real>
% Particle swarm optimization in MPSO/D
% c1   ---   2 --- Parameter in updating particle's velocity
% c2   ---   2 --- Parameter in updating particle's velocity
% CR   --- 0.5 --- Parameter CR in differental evolution
% F    --- 0.5 --- Parameter F in differental evolution
% proM ---   1 --- The expectation of number of bits doing mutation 
% disM ---  20 --- The distribution index of polynomial mutation

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
%     if nargin < 8
%         [c1,c2,CR,F,proM,disM,CA] = deal(Parameter{:});
%     else
%         [c1,c2,CR,F,proM,disM,CA] = deal(1.42,1.63,0.5,0.5,1,20,0.5);
%     end
    ParticleDec = Particle(:,Parent);
    PbestDec = Particle(:,Pbest);
    GbestDec = Particle(:,Gbest);
    [N,D] = size(ParticleDec);
    
    Y_Index = zeros(size(Parent));
    Y_Index(Parent==Pbest) = 1;
    Y_Index(Parent==Gbest) = -1;
    Y_Index(Parent==Gbest &Parent==Pbest & rand(size(Parent))>0.5) = 1;
    Y_Index = repmat(Y_Index,N,1);
    
%     Lower = repmat(lower,N,D);
%     Upper = repmat(upper,N,D);
    r1 = repmat(rand(N,1),1,D);
    r2 = repmat(rand(N,1),1,D);
    ParticleVel = W.*Parent_Vel + c1.*r1.*(-1-Y_Index) + c2.*r2.*(1-Y_Index);
    L_index = Y_Index + ParticleVel;
    Y2_Index = zeros(size(Y_Index));
    Y2_Index(L_index>CA)=1;
    Y2_Index(L_index<-CA)=-1;
    ParticleDec(Y2_Index==1)=GbestDec(Y2_Index==1);
    ParticleDec(Y2_Index==-1)=PbestDec(Y2_Index==-1);
    ParticleDec(Y2_Index==0)=randi(max(max(Particle)),size(ParticleDec(Y2_Index==0)));
    Offspring = ParticleDec;
    Offspring_Vel = ParticleVel;
%     %% Particle swarm optimization
%     Lower = repmat(lower,N,D);
%     Upper = repmat(upper,N,D);
%     DoPSO = repmat(rand(N,1)<0.5,1,D);
%     
%     r1 = repmat(rand(N,1),1,D);
%     r2 = repmat(rand(N,1),1,D);
%     OffVel = ParticleVel;
%     OffDec = ParticleDec;
%     OffVel(DoPSO) = W.*ParticleVel(DoPSO) + c1.*r1(DoPSO).*(PbestDec(DoPSO)-ParticleDec(DoPSO)) + c2.*r2(DoPSO).*(GbestDec(DoPSO)-ParticleDec(DoPSO));
%     OffDec(DoPSO) = ParticleDec(DoPSO) + OffVel(DoPSO);
%     % Set the infeasible decision variables to the value of their parents
%     Invalid = OffDec < Lower | OffDec > Upper;
%     OffDec(Invalid) = ParticleDec(Invalid);
%     
%     %% DE
%     Site = ~DoPSO & rand(N,D)<CR;
%     OffDec(Site) = ParticleDec(Site) + F.*(GbestDec(Site)-PbestDec(Site));
%     % Set the infeasible decision variables to boundary values
%     OffDec = max(min(OffDec,Upper),Lower);
% 
%     %% Polynomial mutation
%     Site  = rand(N,D) < proM/D;
%     mu    = rand(N,D);
%     temp  = Site & mu<=0.5;
%     OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
%                    (1-(OffDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
%     temp = Site & mu>0.5; 
%     OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
%                    (1-(Upper(temp)-OffDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
%     Offspring = OffDec;
%     Offspring_Vel = OffVel;
end