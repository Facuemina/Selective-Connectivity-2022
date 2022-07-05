function Matrix = Optimize_Connectivity(N,c,Num2Sto,W,PAT,epsilon,Matrix)
%% Create Noise Matrices and optimize
%Define Initial Temperature / Compute Noise Matrices and optimize
    row = 1;
    indices  = find(1:N ~= row);   
    Omega    = W(row,indices).*PAT(1:Num2Sto,indices) - PAT(1:Num2Sto,row)/c;   
    T0       = Find_initial_T(N-1,c,Omega,1000,0,0.8,epsilon,PAT(1:Num2Sto,row));
    %Start Optimization
    parfor row = 1:N                
        indices = find(1:N ~= row);
        Omega   = W(row,indices).*PAT(1:Num2Sto,indices) - PAT(1:Num2Sto,row)/c;   
        [Matrix(row,:), ~, ~] = SimulatedAnnealing(N-1,c,T0,10^(-10),0.99,300,800,Omega, PAT(1:Num2Sto,row),epsilon,0);
    end
end