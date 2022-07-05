function [T0] = Find_initial_T(N,c,Omega,iterations,xinit,prob,epsilon,Pat)

if xinit == 0
    xinit = [ones(1,c,'single'),zeros(1,N-c,'single')];
    xinit = xinit(randperm(N));
end

E      = norm(Omega*xinit' - Pat*epsilon);
WOmega = Omega'*Omega;
DeltaE = zeros(1,iterations);
Signal = 2*Pat'*Omega*epsilon;

for it = 1:iterations
    R          = randperm(N,2);
    DeltaX     = [xinit(R(2))-xinit(R(1)),xinit(R(1))-xinit(R(2))];
    DeltaW     = WOmega(:,R)*DeltaX';
    DeltaE(it) = sqrt(DeltaX*DeltaW(R) + 2*(xinit*DeltaW) + E^2 - Signal(R)*DeltaX') - E;

    xinit(R) = xinit(R) + DeltaX;
    E        = E + DeltaE(it);      
end
T0    = -max(abs(DeltaE))/log(prob);
end