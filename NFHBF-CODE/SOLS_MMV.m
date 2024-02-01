function [SUPP,Gain,Lambda] = SOLS_MMV(Y,A,P,M,Nsub,d,lambda,I)

angle_P=unifrnd(-1,1,P,1);
PHI=[];
for p=1:P
    ANG=repmat(PW(angle_P(p),d,lambda,Nsub),1,size(A,2));
    PHI=[PHI;ANG];
end
%PHI=zeros(Nsub*P,size(A,2));
[N,R]=size(Y);
G=size(PHI,2);
Lambda=[];supp=[];SUPP=[];
%I=12
for l=1:M
    for i=1:I
    for p=1:P
        PHI((p-1)*Nsub+1:p*Nsub,:)=A;
       
            for g=1:G
%                 if ismember(g,Lambda)
%                     g=g+1;
%                 end
                %PHI_G=PHI(:,[Lambda,g]);
                PHI_G=[SUPP,PHI(:,g)];
               % Res(g)=trace((PHI_G'*Y)'*(PHI_G'*Y));
                Res(g)=trace(((eye(N)-PHI_G*((PHI_G'*PHI_G)\PHI_G'))*Y)'*((eye(N)-PHI_G*((PHI_G'*PHI_G)\PHI_G'))*Y));
            end
    
        [~,indm]=min(Res);
        supp=PHI(:,indm);     
        PHI((p-1)*Nsub+1:p*Nsub,:)=repmat(A(:,indm),1,size(A,2));
        Lambda=[Lambda,indm];
    end
    end
    SUPP=[SUPP, supp];
end
    Gain=((SUPP'*SUPP)\SUPP')*Y;
end
