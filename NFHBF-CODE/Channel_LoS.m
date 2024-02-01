function H_LoS = Channel_LoS(Nt,Nr,r0,d,lambda)

for nr=1:Nr
    for nt=1:Nt
        r(nr,nt)=sqrt((nr-nt)^2*d^2+r0^2);
        H_LoS(nr,nt)=1/r(nr,nt)*exp(1i*2*pi/lambda*r(nr,nt));
    end
end
H_LoS=H_LoS/sqrt(trace(H_LoS'*H_LoS)); %gain normalization
end

