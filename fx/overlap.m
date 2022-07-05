function [m, frac] = overlap(Pattern,Num2Sto,c,C,evaluate)

    N = size(Pattern,2);
    W = transpose(Pattern(1:Num2Sto,:))*Pattern(1:Num2Sto,:);       
    J = C .*W/c; %Jij = Cij*Wij %Full Weight Matrix
    
    if evaluate == 0
        
        Neurons     = transpose(Pattern(1:Num2Sto,:));
        m           = zeros(Num2Sto,100,'single');
        m(:,1)      = sum(Pattern(1:Num2Sto,:).*Neurons',2)/N;
        for i = 2:100
            Neurons = sign(J*Neurons);
            m(:,i)  = sum(Pattern(1:Num2Sto,:).*Neurons',2)/N;
        end

        frac = length(find(m(:,end)>0.7))/Num2Sto;
        
    else
        index    = randperm(Num2Sto,evaluate);
        Neurons  = transpose(Pattern(1:index,:));
        m        = zeros(Num2Sto,100,'single');
        m(:,1)   = sum(Pattern(1:evaluate,:).*Neurons',2)/N;
        
        for i = 2:100
            Neurons = sign(J*Neurons);
            m(:,i)  = sum(Pattern(1:evaluate,:).*Neurons',2)/N;
        end
        frac = length(find(m(:,end)>0.7))/evaluate;
    end
end