function [DeltaE,T0,ErrT0] = Copy_of_Test_T0(N,c,Omega,iterations,xinit,prob,mod,epsilon,Pat)

if xinit == 0
    xinit = [ones(1,c,'single'),zeros(1,N-c,'single')];
    xinit = xinit(randperm(N));
end

E      = norm(Omega*xinit' - Pat*epsilon);
DeltaX = [-1,1];
WOmega = Omega'*Omega;
DeltaE = zeros(1,iterations);
Signal  = 2*Pat'*Omega*epsilon;

for it = 1:iterations

    active = find(xinit);
    particle = active(randi(length(active)));
    exchange = find(xinit<mod);
    exchange(find(exchange == particle)) = [];
    exchange = exchange(randi(length(exchange)));

    R = [particle, exchange]; %R = sort(R,'ascend');    

    DeltaW     = WOmega(:,R)*DeltaX';
    DeltaE(it) = sqrt(DeltaX*DeltaW(R) + 2*(xinit*DeltaW) + E^2 - Signal(R)*DeltaX') - E;

    xinit(R) = xinit(R) + DeltaX;
    E        = E + DeltaE(it);      

end
T0 = -mean(abs(DeltaE))/log(prob);
ErrT0 = abs(std(abs(DeltaE))/log(prob));
end