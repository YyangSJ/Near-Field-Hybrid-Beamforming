function [SUPP,Gain,Lambda] = LC_SOLS_MMV(Y,A,P,M,Nsub,d,lambda) 

[N,R]=size(Y);
G=size(A,2);
Lambda=[];SUPP=[];
for m=1:M
    supp_all=[];
    for p=1:P
        for g=1:G
            if m==1
                PHI_G=A(:,g)
            else
                PHI_G=[SUPP((p-1)*Nsub+1:p*Nsub,:),A(:,g)];
            end
            Res(g)=trace(((eye(Nsub)-PHI_G*((PHI_G'*PHI_G)\PHI_G'))*Y((p-1)*Nsub+1:p*Nsub,:))'...
                *((eye(Nsub)-PHI_G*((PHI_G'*PHI_G)\PHI_G'))*Y((p-1)*Nsub+1:p*Nsub,:)));
        end
        [~,indm]=min(Res);
        supp=A(:,indm);
        supp_all=[supp_all;supp]
    end
    
    SUPP=[SUPP, supp_all];
end
Gain=((SUPP'*SUPP)\SUPP')*Y;
end
