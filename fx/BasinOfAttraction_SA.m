%% Basin Of Attraction - Simulated Annealing Optimization
function [BOA,C] = BasinOfAttraction_SA(N,c,Pattern,epsilon,evaluate,Initial_Error_Perc)

    W = Pattern'*Pattern/c;
    C                 = zeros(N,N,'single');
    Matrix            = zeros(N,N-1,'single');
    %Indices for lower and upper triangular parts of matrices
    LowerTriJ         = find(tril(ones(size(W)),-1));   
    UpperTriJ         = find(triu(ones(size(W)),1));
    LowerTriMatrix    = find(tril(ones(size(Matrix)),-1));
    UpperTriMatrix    = find(triu(ones(size(Matrix))));
    
    Matrix = Optimize_Connectivity(N,c,size(W,1),W,Pattern,epsilon,Matrix);
    
    C(LowerTriJ) = Matrix(LowerTriMatrix);
    C(UpperTriJ) = Matrix(UpperTriMatrix);
    %Estimate Basin Of Attraction for C
    
    BOA  = zeros(1,length(Initial_Error_Perc)); 
    for perc = 1:length(Initial_Error_Perc)            
        [~, BOA(perc)] = overlap_err(W.*C,Pattern,length(Pattern(:,1)),evaluate,Initial_Error_Perc(perc));
        if BOA(perc) == 0
            break
        end
    end
end