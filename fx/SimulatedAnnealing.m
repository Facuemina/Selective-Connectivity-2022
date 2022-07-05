function [xinit, fval, sol] = SimulatedAnnealing(N,c,Ti,Tf,alpha,preterm,stop,Omega,Pat,epsilon,xinit)

if xinit == 0
    xinit = [ones(1,c),zeros(1,N-c)];
    xinit = xinit(randperm(N));
end

E       = norm(Omega*xinit' - Pat*epsilon);
randvec = 1:N;
WOmega  = Omega'*Omega;
Signal  = 2*Pat'*Omega*epsilon;
sol     = 1:1:stop;
i       = stop;

T = Ti;
while T > Tf

    i=i+1; 

    for it_pret = 1:preterm

        R = randvec(randperm(N,2));
        if xinit(R(1))-xinit(R(2)) == 0
            continue
        else
            DeltaX = [xinit(R(2))-xinit(R(1)),xinit(R(1))-xinit(R(2))];
        end  

        DeltaW = WOmega(:,R)*DeltaX';
        DeltaE = sqrt(DeltaX*DeltaW(R) + 2*(xinit*DeltaW) + E^2 - Signal(R)*DeltaX') - E;

        if DeltaE < 0 
            xinit(R) = xinit(R) + DeltaX;
            E        = E + DeltaE;
        elseif rand<exp(-DeltaE/T)
            xinit(R) = xinit(R) + DeltaX;
            E        = E + DeltaE;
        end    

    end

    sol(i) = E;    
    if isequal(sol(end-stop+1:end),ones(1,stop)*sol(end)) || E < 10^(-4)
        break
    end
    T = alpha*T;
end

fval = sol(end);
sol = sol(stop+1:end);
end