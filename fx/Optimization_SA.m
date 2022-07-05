
%% Function: Connectivity optimization of Hopfield Network using a Simulated Annealing algorithm
function [PMAX,DataOutput] = Optimization_SA(N,c,PAT,parameters,e)
%% INPUTS:
% N          = network dimension
% c          = network connectivity
% PAT        = matrix containing random patterns
% parameters = [initial guess, Delta]
% e	         = defines noise bias to optimize (e=0 -> epsilon = 0; e = 1 -> epsilon = p/c)

%% OUTPUTS:
% PMAX                  = maximum number of patterns an optimized network can efficiently store
% DataOutput.pmax       = memories that the network can handle
% DataOutput.PAT        = set of random patterns used
% DataOutput.PreLoading = values of preloaded memories
% DataOutput.CMSq       = Conectivity matrix optimized for max(DataOutput.pmax)

    %% First Step: DataOutputerate variables
    DataOutput                            = struct;
    DataOutput.pmax                       = [];
    %Patterns to use
    if length(PAT)<2
        DataOutput.PAT = sign(rand(N,N)-0.5);
    else
        DataOutput.PAT = PAT;
    end
    
    %% Second Step: begins optimization
    Delta = parameters(2); lb = 0; ub = N;
    
    DataOutput.pmax(1)       = parameters(1) - Delta;  
    DataOutput.PreLoading(1) = parameters(1) - Delta;
    DataOutput.CMSq          = zeros(N,N,'single');
    %Useful matrices
    Matrix            = zeros(N,N-1,'single');
    J                 = zeros(N,N,'single');
    %Indices for lower and upper triangular parts of matrices
    LowerTriJ         = find(tril(ones(size(J)),-1));   
    UpperTriJ         = find(triu(ones(size(J)),1));
    LowerTriMatrix    = find(tril(ones(size(Matrix)),-1));
    UpperTriMatrix    = find(triu(ones(size(Matrix))));
    iter = 1;
    %% Deinfe upper and lower limits of binary search
    % Initial guess is given in parameters(1). Bounds are initially set to
    % [0,N]. If the optimized network is capable of loading the initial
    % guess number of memories, then it is set as a lower bound and search
    % of upper bound continues with constant Delta steps (=parameters(2)) 
    % till the network can't handle that amount of memories. On the other 
    % hand, if the optimized network can not handle the initial amount of 
    % memories, then it is set as an upper bound and search of lower bound 
    % continues with constant Delta steps
    
    Num2Sto = DataOutput.PreLoading(1); 
    W       = (DataOutput.PAT(1:Num2Sto,:))'*(DataOutput.PAT(1:Num2Sto,:))/c;  
    while ~all([lb,ub]-[0,N]) 
            
        DataOutput.PreLoading(iter+1) = DataOutput.PreLoading(iter) + Delta;
        DataOutput.pmax(iter+1)       = DataOutput.PreLoading(iter+1);
        Num2Sto                       = DataOutput.PreLoading(iter+1);  
        
        min_load = min(Num2Sto-Delta,Num2Sto); max_load = max(Num2Sto-Delta,Num2Sto);
        W        = W + sign(Delta)*(DataOutput.PAT((min_load+1):max_load,:))'*(DataOutput.PAT((min_load+1):max_load,:))/c;        
        iter = iter + 1;
        
        epsilon  = e*(Num2Sto)/c;
        Matrix   = Optimize_Connectivity(N,c,Num2Sto,W,DataOutput.PAT,epsilon,Matrix);
        
        %Test memory capacity
        J(LowerTriJ) = Matrix(LowerTriMatrix);
        J(UpperTriJ) = Matrix(UpperTriMatrix);
        J       = J.*W; %Jij = Cij*Wij %Full Weight Matrix
        
        DataOutput.pmax(iter) = Test_LoadingCapacity(J,DataOutput.PAT,Num2Sto);
        
        disp(['Preload = ',num2str(DataOutput.PreLoading(iter)),'; Pmax = ',num2str(DataOutput.pmax(iter))])
        
        if DataOutput.PreLoading(iter) == DataOutput.pmax(iter)
            DataOutput.CMSq(LowerTriJ) = Matrix(LowerTriMatrix);
            DataOutput.CMSq(UpperTriJ) = Matrix(UpperTriMatrix);
            lb                         = DataOutput.PreLoading(iter);
            Delta                      = abs(Delta);
        else
            ub                         = DataOutput.PreLoading(iter);
            Delta                      = -abs(Delta);
        end

    end 
    
    %% Start Binary Search
    % Given the initial upper and lower bounds, start a binary search
    Delta = round(Delta/2);
    while ub - lb > 1    

        DataOutput.PreLoading(iter+1) = DataOutput.PreLoading(iter) + Delta;
        DataOutput.pmax(iter+1)       = DataOutput.PreLoading(iter+1);
        Num2Sto                       = DataOutput.PreLoading(iter+1);         
        
        min_load = min(Num2Sto-Delta,Num2Sto); max_load = max(Num2Sto-Delta,Num2Sto);
        W        = W + sign(Delta)*(DataOutput.PAT((min_load+1):max_load,:))'*(DataOutput.PAT((min_load+1):max_load,:))/c;        
        iter = iter + 1;

        epsilon  = e*(Num2Sto)/c;
        Matrix   = Optimize_Connectivity(N,c,Num2Sto,W,DataOutput.PAT,epsilon,Matrix);

        %Test memory capacity
        J(LowerTriJ) = Matrix(LowerTriMatrix);
        J(UpperTriJ) = Matrix(UpperTriMatrix);
        J            = J.*W; %Jij = Cij*Wij %Full Weight Matrix
        
        DataOutput.pmax(iter) = Test_LoadingCapacity(J,DataOutput.PAT,Num2Sto);

        disp(['Preload = ',num2str(DataOutput.PreLoading(iter)),'; Pmax = ',num2str(DataOutput.pmax(iter))])

        if DataOutput.PreLoading(iter) == DataOutput.pmax(iter)
            DataOutput.CMSq(LowerTriJ)  = Matrix(LowerTriMatrix);
            DataOutput.CMSq(UpperTriJ)  = Matrix(UpperTriMatrix);
            lb                          = DataOutput.PreLoading(iter);
            Delta                       = round((ub-lb)/2);
        else
            ub                          = DataOutput.PreLoading(iter);
            Delta                       = -round((ub-lb)/2);
        end

    end        

    index = logical(DataOutput.pmax - DataOutput.PreLoading == 0);
    PMAX  = max(DataOutput.pmax(index));
    disp([c,N])
end