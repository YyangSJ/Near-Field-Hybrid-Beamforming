function [channel_SW,phi,theta,r] = Channel_realization(d,lambda,L,Nr,Nt)
r0=unifrnd(5,100);

theta=unifrnd(-1,1,L-1,1);
phi=unifrnd(-1,1,L-1,1);
r=unifrnd(5,100,L-1,1);
r2=unifrnd(5,100,L-1,1);
z=1/sqrt(2)*(normrnd(0,1,L-1,1)+1i*normrnd(0,1,L-1,1));
channel_SW=1/sqrt(2)*(normrnd(0,1)+1i*normrnd(0,1))*sqrt(Nr*Nt)*Channel_LoS(Nt,Nr,r0,d,lambda);

for l=1:L-1
    channel_SW=channel_SW+sqrt(Nr*Nt/L)*z(l)*SW(phi(l),r2(l),d,lambda,Nr)*SW(theta(l),r(l),d,lambda,Nt)';
end

end

