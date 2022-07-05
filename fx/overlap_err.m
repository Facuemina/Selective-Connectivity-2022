%% Overlap Error
% Returns overlap (m) and fraction of memories retrived with an err*100%
% percentage of initial error after a 100 iteration process
function [m, frac] = overlap_err(J,Pattern,Num2Sto,evaluate,err)

    N = size(J,2);
    %W = transpose(Pattern(1:Num2Sto,:))*Pattern(1:Num2Sto,:);       
    %J = C .*W/c; %Jij = Cij*Wij %Full Weight Matrix
    
    if evaluate == 0
        
        error       = round(N*err);
        Neurons     = transpose(Pattern(1:Num2Sto,:));
        
        if error >= 1 
            init_err            = randperm(N,error);
            Neurons(init_err,:) = -Neurons(init_err,:);
        end
        
        m            = zeros(Num2Sto,100,'single');
        m(:,1)       = sum(Pattern(1:Num2Sto,:).*Neurons',2)/N;
        
        for i = 2:100
            Neurons  = sign(J*Neurons);
            m(:,i)   = sum(Pattern(1:Num2Sto,:).*Neurons',2)/N;
        end

        frac = length(find(m(:,end)>0.7))/Num2Sto;
        
    else
        error    = round(N*err);
        index    = randperm(Num2Sto,evaluate);
        Neurons  = transpose(Pattern(index,:));
        
        if error >= 1
            init_err            = randperm(N,error);
            Neurons(init_err,:) = -Neurons(init_err,:);
        end
        
        m        = zeros(evaluate,100,'single');
        m(:,1)   = sum(Pattern(index,:).*Neurons',2)/N;
        
        for i = 2:100
            Neurons  = sign(J*Neurons);
            m(:,i)   = sum(Pattern(index,:).*Neurons',2)/N;
        end

        frac = length(find(m(:,end)>0.7))/evaluate;
    end
end