function [SUPP,Gain,Lambda] = SOMP_MMV(Y,A,P,M,Nsub,d,lambda,I)
angle_P=unifrnd(-1,1,P,1);
PHI=[];
for p=1:P
    ANG=repmat(PW(angle_P(p),d,lambda,Nsub),1,size(A,2));
    PHI=[PHI;ANG];
end
%PHI=zeros(Nsub*P,size(A,2));
[N,R]=size(Y);
Y_Res=Y;
G=size(PHI,2);
Lambda=[];supp=[];SUPP=[];

for l=1:M
    for i=1:I
        for p=1:P
            PHI((p-1)*Nsub+1:p*Nsub,:)=A;
            
            for g=1:G
                %                 if ismember(g,Lambda)
                %                     g=g+1;
                %                 end
                tep(g)=sum(abs(PHI(:,g)'*Y_Res));
                %     Res(g)=trace(((eye(N)-PHI_G*((PHI_G'*PHI_G)\PHI_G'))*Y)'*((eye(N)-PHI_G*((PHI_G'*PHI_G)\PHI_G'))*Y));
            end
            
            [~,indm]=max(tep);
            supp=PHI(:,indm);
            PHI((p-1)*Nsub+1:p*Nsub,:)=repmat(A(:,indm),1,size(A,2));
            Lambda=[Lambda,indm];
        end
    end
    SUPP=[SUPP, supp];
    Y_Res=Y-SUPP*((SUPP'*SUPP)\SUPP')*Y;
end 

Gain=((SUPP'*SUPP)\SUPP')*Y; 
end
