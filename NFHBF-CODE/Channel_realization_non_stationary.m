function [channel_SW,phi,theta,r] = Channel_realization_non_stationary(d,lambda,L,Nr,Nt,P,vn)

theta=unifrnd(-1,1,L,1); 
phi=unifrnd(-1,1,L,1); 
r=unifrnd(5,100,L,1);
r2=unifrnd(5,100,L,1);
z=1/sqrt(2)*(normrnd(0,1,L,1)+1i*normrnd(0,1,L,1));
 
r0=unifrnd(5,100);
channel_SW=1/sqrt(2)*(normrnd(0,1)+1i*normrnd(0,1))*sqrt(Nr*Nt)*Channel_LoS(Nt,Nr,r0,d,lambda);
%channel_SW=zeros(Nr,Nt);
for l=1:L-1
    v=ones(Nt,1); 
    NVR=randperm(P,vn);
    for i=1:vn
          v((NVR(i)-1)*Nt/P+1:NVR(i)*Nt/P)=zeros(Nt/P,1);
    end 
    channel_SW=channel_SW+sqrt(Nr*Nt/L)*z(l)*SW(phi(l),r2(l),d,lambda,Nr)*(SW(theta(l),r(l),d,lambda,Nt).*v)';
end

end

