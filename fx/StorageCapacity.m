%% Storage Capacity
% Computes the Storage Capacity of a network with N nodes and c
% connections per node. If C == 0, then de function generates a random
% conectivity matrix with c presynaptic connections per neuron. If Patterns == 0, then
% the function generates N random patterns to be loaded to the network.
% The function starts loading p0 = max(1,round(c*0.14)-1) patterns and adds
% Delta = p0 patterns in every iteration. Whenever the network can't handle
% that memory load, Delta is reduced by a half and the loop restarts
% from the last fully handled memory load. The loops continue untill
% Delta = 0 (or the difference between the upper limit and lower limit is
% 1). 
% Criteria: given a set of p patterns. The network is able to retrive them
% all if the overlap bewtween the corresponding attractor to that memory
% and the actual network state is more than 70%
%%
function [pmax,C] = StorageCapacity(N,c,C,Patterns) 
%#Codegen
p            = 1;
Norm         = c;
%%
m1 = 1;
%% Patterns
if Patterns == 0
    if N>=1000
        Patterns = sign(rand(N,round(N*0.5))-0.5);
    else
        Patterns = sign(rand(N,N)-0.5);
    end
end
%% Connectivity
if C == 0
    C = zeros(N,N,'single');
    v = [ones(1,c),zeros(1,N-1-c)];
    for i = 1:N
        C(i,1:N ~= i) = v(randperm(N-1));
    end
end

%% Pre-Store patterns to the network
Delta = max(1,round(c*0.14)-1);
lb = 1; ub = N;
while Delta >= 1 && ub-lb > 1

    while all(m1>0.7)   
        lb = p;

        p  = p + Delta;

        W       = transpose(Patterns(1:p,:))*Patterns(1:p,:);       
        J       = C .*W/Norm; %Jij = Cij*Wij %Full Weight Matrix
        Neurons = transpose(Patterns(1:p,:));
        h = 6;
        m = zeros(p,100,'single');

        for i = 1:h
            Heff     = J*Neurons;
            Neurons = sign(Heff);
            m(:,h)   = sum(Patterns(1:p,:).*Neurons',2)/N;
        end

        while any(abs(m(:,h)-m(:,h-5))>0)
            h = h+1;
            Heff     = J*Neurons;
            Neurons = sign(Heff);
            m(:,h)   = sum(Patterns(1:p,:).*Neurons',2)/(N);
            if any(m(:,h)<0.7) || h>100
                break
            end 
        end
        m1 = m(:,h);
    end

    ub = p;
    p  = p - Delta;
    Delta = round(Delta/2 - 0.1);
    m1=1;

end
pmax = p;