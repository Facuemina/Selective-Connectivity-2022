function p = Test_LoadingCapacity(J,Patterns,Num2Sto)

    Neurons = transpose(Patterns(1:Num2Sto,:));
    N       = length(Neurons);
    m       = zeros(Num2Sto,100,'single');

    for i = 1:100
        Heff     = J*Neurons;
        Neurons = sign(Heff);
        m(:,i)   = sum(Patterns(1:Num2Sto,:).*Neurons',2)/N;
    end
    p = sum(m(:,end)>0.7);
end