function [Results] = CEMP(path_data,index_of_similarity,CA,C1,C2,T_index,evaluation)
%% ==========================================================================
% % FUNCTION: CEMP
% % DESCRIPTION: A pipline for multi-view clustering via multiobjective
% % particle swarm optimization

% % ==========================================================================
% % % % copyright (c) 2019 Q.M Liu & G.H Wang & X.D Zhao
% % ==========================================================================
    tic
    %% Input data and data pre-processing
    Data=load(path_data);
    similarity={'ED','PC','SC','PC_ED','SC_ED','PC_SC','PC_SC_ED'};
    % index_of_similarity=4;

    %% 输入数据，数据预处理
    X=Data.in_X;
    X_norm=normalizeData(X);
    num_Cluster = length(unique(Data.true_labs)); % the number of clusters in the final clustering (using in consensus functions)
    [n,m]=size(X);

    labelParicles=[];
    for i = 2: 20
        option.algorithm='nmfrule';
        [A,Y]=nmf(X_norm',i,option); 
        Cl = kmeans(Y',num_Cluster);
        switch index_of_similarity
            case 1
                CL1 = kmeans(Y',num_Cluster);
                labelParicles = [labelParicles CL1];
                if i ==20
                    fprintf('SIMILARITY: %s',similarity{index_of_similarity});
                end
            case 2
                [CL2,~,~] = K_means(Y',num_Cluster, 500, 'persondist');
                labelParicles = [labelParicles CL2'];
                if i ==20
                    fprintf('SIMILARITY: %s',similarity{index_of_similarity});
                end
            case 3
                [CL3,~,~] = K_means(Y',num_Cluster, 500, 'spermandist');
                labelParicles = [labelParicles CL3'];
                if i ==20
                    fprintf('SIMILARITY: %s',similarity{index_of_similarity});
                end
            case 4
                CL1 = kmeans(Y',num_Cluster);
                [CL2,~,~] = K_means(Y',num_Cluster, 500, 'persondist');
                labelParicles = [labelParicles CL1 CL2'];
                if i ==20
                    fprintf('SIMILARITY: %s',similarity{index_of_similarity});
                end
            case 5
                CL1 = kmeans(Y',num_Cluster);
                [CL3,~,~] = K_means(Y',num_Cluster, 500, 'spermandist');
                labelParicles = [labelParicles CL1 CL3'];
                if i ==20
                    fprintf('SIMILARITY: %s',similarity{index_of_similarity});
                end
            case 6
                [CL2,~,~] = K_means(Y',num_Cluster, 500, 'persondist');
                [CL3,~,~] = K_means(Y',num_Cluster, 500, 'spermandist');
                labelParicles = [labelParicles CL2' CL3'];
                if i ==20
                    fprintf('SIMILARITY: %s',similarity{index_of_similarity});
                end
            case 7
                CL1 = kmeans(Y',num_Cluster);
                [CL2,~,~] = K_means(Y',num_Cluster, 500, 'persondist');
                [CL3,~,~] = K_means(Y',num_Cluster, 500, 'spermandist');
                labelParicles = [labelParicles CL1 CL2' CL3'];
                if i ==20
                    fprintf('SIMILARITY: %s',similarity{index_of_similarity});
                end
        end %end of switch
    end %end of for i =2:20

    %% Generate the weight vectors
    M=3;
    N=size(labelParicles,2);
    [W,N] = UniformPoint(N,M);
    W = W./repmat(sqrt(sum(W.^2,2)),1,size(W,2));
           %T计算方法1：
    if T_index == 1
        T = ceil(N/3);
        if T<3
            T=3;
        end
    else
        T=4;
    end
    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);

    %% Generate random population
    Lower = 0;
    Upper = 1;
    Population1 = labelParicles;
    [tn,tm]=size(Population1);

    Value1 = [];
    for i = 1:tm
        value_cp(i)=valid_compactness(X_norm,Population1(:,i));
        value_Dea(i)=Devation(X_norm, Population1(:,i));
        [Sinter]=valid_sumsqures(X_norm, Population1(:,i), num_Cluster);
        value_sum(i)=sum(sum(Sinter));
        TEMP = [value_cp(i) value_Dea(i) value_sum(i)];
        Value1=[Value1; TEMP];
    end

    Z1 = min(Value1,[],1);
    PopObj = Value1 - repmat(Z1,size(Population1,2),1);
    true_PopObj  = PopObj*W';
    [~,P]  = max(true_PopObj,[],2);

    Population = labelParicles(:,P(1:N));
    Value = Value1(P(1:N),:);
    Z = min(Value,[],1);
    W_plus=W;

    %% Optimization
    B_plus=B;
    %evaluation = 100; %the number of iteration
    %% Paramaters intilization
    Offspring_Vel = zeros(size(Population));
    POP_best=zeros(1,evaluation);
    bestPOP=zeros(size(Population,1),evaluation);
    GAL_best=zeros(1,evaluation);
    gbestPOP=zeros(size(Population,1),evaluation);

    disp('START EVALUTION');
    all_cluster_results1=[];
    all_cluster_results2=[];
    Results=[];
    %% Evaluation start
    for evaluated = 1:evaluation
        if mod(evaluated,10) == 0
            fprintf('SIMILARITY: %s, T= %d, C1 = %d, C2 = %d, CA = %d,evaluation: %d\n',similarity{index_of_similarity},T,C1,C2,CA, evaluated);
        end
        [Parent,Pbest,Gbest] = MatingSelection(Value,B_plus,W_plus,Z);
        w_1 = 0.4 + (evaluation-evaluated/evaluation)*0.5;
        [Offspring, Offspring_Vel] = Operator(C1,C2,CA,w_1,Population,Parent,Pbest,Gbest,Offspring_Vel);%fixed C1, C2, and CA
        all_cluster_results1(:,:,evaluated)=Offspring;
        [fn,fm]=size(Offspring);
        Offspring_Value = [];
        Offspring_Value_with_label = [];
        for i = 1:fm
            % objective functions without labels
            Offspring_value_cp(i)=valid_compactness(X_norm,Offspring(:,i));
            Offspring_value_Dea(i)=Devation(X_norm, Offspring(:,i));
            [Sinter]=valid_sumsqures(X_norm, Offspring(:,i), num_Cluster);
            Offspring_value_sum(i)=sum(sum(Sinter));
            TEMP = [Offspring_value_cp(i) Offspring_value_Dea(i) Offspring_value_sum(i)];
            Offspring_Value=[Offspring_Value; TEMP];

            % objective functions with labels
            Offspring_NMI(i)=Cal_NMI(Data.true_labs, Offspring(:,i));
            Offspring_RI(i)=RandIndex(Data.true_labs, Offspring(:,i));
            TEMP1 = [Offspring_NMI(i) Offspring_RI(i)];
            Offspring_Value_with_label = [Offspring_Value_with_label; TEMP1];
        end
        all_Offspring_Value(:,:,evaluated)=Offspring_Value;
        all_Offspring_Value_with_label(:,:,evaluated)=Offspring_Value_with_label;
        Z = min([Z;Offspring_Value],[],1);

        %计算每次最好的个体和全局最好个体
        [POP_best(evaluated),index]=max(Offspring_Value_with_label(:,1));
        bestPOP(:,evaluated) = Offspring(:,index);
        if evaluated == 1
            GAL_best(evaluated)=POP_best(evaluated);
            gbestPOP(:,evaluated)=bestPOP(:,evaluated);
        elseif GAL_best(evaluated-1)<POP_best(evaluated)
            GAL_best(evaluated)=POP_best(evaluated);
            gbestPOP(:,evaluated)=bestPOP(:,evaluated);
        else
            GAL_best(evaluated)=GAL_best(evaluated-1);
            gbestPOP(:,evaluated)=gbestPOP(:,evaluated-1);
        end

        Population = Classification(num_Cluster,[Population Offspring],[Value; Offspring_Value],W_plus,Z);
        all_cluster_results2(:,:,evaluated)=Population;
        toc
    end %end of evaluated
    [Results.best,Results.index]=max(GAL_best);
    Results.withoutlabel=all_Offspring_Value;
    Results.withlabel=all_Offspring_Value_with_label;
    Results.cluster=all_cluster_results2;
    Results.finalCluster=gbestPOP(:,Results.index);
    [rs,jc,fm,f1]=evalution(X, Data.true_labs, Results.finalCluster);
    nmi=Cal_NMI(Data.true_labs, Results.finalCluster);
    ri=RandIndex(Data.true_labs, Results.finalCluster);
    Results.evalution=[nmi,ri,rs,jc,fm,f1];
    toc
