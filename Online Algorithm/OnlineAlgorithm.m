function [C,MeanConnectivity,StdConnectivity,PEFF,iterations] = OnlineAlgorithm(cinit,stop,Patterns,DeltaP,e)
%% Start Online Algorithm

%Num2Sto = initial pattern load
%C       = initial random connectivity matrix
[Num2Sto,C] = StorageCapacity(size(Patterns,2),cinit,0,Patterns);

iterations = 1;
MeanConnectivity = cinit; StdConnectivity = 0; PEFF = Num2Sto;

while PEFF(end)>=Num2Sto%0.9    
    disp(Num2Sto)
    Num2Sto = Num2Sto+DeltaP;
    
    [C,peff,CMean,Cstd] = online_connectivity_optimization(mean(sum(C,2)),C,Patterns(1:Num2Sto,:),e,stop,1);%0.9
    
    iterations       = [iterations,length(peff)];
    MeanConnectivity = [MeanConnectivity,CMean];
    StdConnectivity  = [StdConnectivity,Cstd];
    PEFF             = [PEFF,peff];    
end

%% Online Connectivity Optimization function
function [C,peff,Cmean,Cstd] = online_connectivity_optimization(c,C,Patterns,e,stop,f)    
    %% Initialize parameters and useful matrices
    N = size(Patterns,2); %Number of neurons
    P = size(Patterns,1); %Number of stored patterns
    
    %Initialize random connectivity matrix if none is provided
    if C == 0
        C = zeros(N,N);
        x = [ones(1,c),zeros(1,N-c-1)];
        for k = 1:N
            C(k,1:N ~= k) = x(randperm(N-1));
        end
    end

    % HS: PxN matrix containing the total local field each neuron receives when
    % retriving exactly a pattern, multiplied by the sign that neuron has in
    % that pattern   
    W  = Patterns'*Patterns;
    HS = (Patterns*(W.*C')).*Patterns;
      
    epsilon = e(1)*P + e(2);
    E = sqrt(sum((HS-epsilon-e(3)*sum(C,2)').^2,1)); %Energy to minimize

    Cmean = 1:stop; %Cmean will store mean pre-synaptic connectivity
    %Cmean(1) = mean(sum(C,2));
    Cstd = Cmean; %Cstd(1) = std(sum(C,2)); %Cstd will store the std of the pre-synaptic connectivity
    peff = [0]; %peff will store the number of patterns the network can retreive (from the P stored)
    
    i = 0;    
    rows = 1:N;
    
    %% Start iteration
    while peff(end) < P*f

        i=i+1; 

            for it_pret = 1:10
                
                %For each neuron (i=1,...,N), select randomly one
                %pre-synaptic neuron (no matter if it's connected or not)
                random_index = randi(N-1,1,N); random_index = random_index + (random_index>=rows); %presynaptic neruon subindex
                R = (random_index-1)*N + rows; %neuron index for connectivity matrix                
                
                %Compute the difference in each element of the HS matrix
                %after changing the connectivity state of each randomly
                %selected pre-synaptic neuron
                DeltaHS = ((2*(0.5-C(R)).*W(R)).*Patterns(:,random_index)).*Patterns;
                
                %Compute the change in energy if the connectivity change is made 
                DeltaE = sqrt(sum((HS+DeltaHS-epsilon-e(3)*(sum(C,2)' + 2*(0.5-C(R)))).^2,1)) - E;
                
                %Keep those changes that minimize local energy
                chgs      = (DeltaE < 0); %boolean: 1 if change is accepted and 0 if not
                E(chgs)   = E(chgs) + DeltaE(chgs); 
                HS(:,chgs)= DeltaHS(:,chgs) + HS(:,chgs);
                C(R(chgs))= 1-C(R(chgs));

            end    
            
        Cmean(i)= mean(sum(C,2));
        Cstd(i) = std(sum(C,2));
        peff(i) = Test_LoadingCapacity(C.*W./sum(C,2),Patterns,P);
        
        %If mean connectivity does not change or iterations reach maximum
        %given by "stop", break optimization
        if isequal(Cmean(end)*ones(1,stop),Cmean(end-stop+1:end)) ||i>10000
            break
        end

    end
    
    %clean arrays
    if i<stop
        Cmean = Cmean(1:i);
        Cstd  = Cstd(1:i);
    end
end
end