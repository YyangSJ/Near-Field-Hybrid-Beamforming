clear all;
close all;
warning off;
%% Some notes can be found in file 'NFHBF_Fig_DFTvsPolar.m'.
rng(0)
f=28e9;
c=3e8;
lambda=c/f;
d=lambda/2;
Nr=8;Nt=128;
P=8;
Nsub=Nt/P;
M=4;
L=4;
Ns=4;
R_sample=0.01:0.01:0.2;
r_sample=1./R_sample;
Realization=200
SNR=-30:5:20;
I=8  ;
vn=5;
CS_NCHT_1=expandMatrix([1,1;1,-1],4);
CS_NCHT_h1= n_product(I_line(Nsub),n_product(blkdiag(I_line(Nsub/2),I_prime(Nsub/2)),n_product(blkdiag(I_line(Nsub/4),I_prime(Nsub/4),...
    I_prime(Nsub/4),I_line(Nsub/4)),CS_NCHT_1,2),4),8);
CS_NCHT_2=expandMatrix([1,1;1i,-1i],4);
CS_NCHT_h2= n_product(I_line(Nsub),n_product(blkdiag(I_line(Nsub/2),I_prime(Nsub/2)),n_product(blkdiag(I_line(Nsub/4),I_prime(Nsub/4),...
    I_prime(Nsub/4),I_line(Nsub/4)),CS_NCHT_2,2),4),8);
Q=zeros(Nsub,Nsub);
for n1=1:Nsub
    n2=bin2dec(flip(dec2bin(n1-1,log2(Nsub))))+1;
    Q(n1,n2)=1;
end

CS_SCHT_h1=Q*CS_NCHT_h1;
CS_SCHT_h2=Q*CS_NCHT_h2;
CS_SCHT_4bit=CS_SCHT_h1;
CS_SCHT_5bit=[CS_SCHT_h1 CS_SCHT_h2];

A_4bit=ULA_FFDic(d,lambda,Nsub,-1+2/(2^4):2/(2^4):1);
A_5bit=ULA_FFDic(d,lambda,Nsub,-1+2/(2^5):2/(2^5):1);
for reali=1:Realization
    reali
    % H=Channel_realization(d,lambda,L,Nr,Nt);
    H=Channel_realization_non_stationary(d,lambda,L,Nr,Nt,P,vn); 
    [~,Sig,V]=svd(H); 
    OV=V(:,1:Ns);
    %% Angle-Stepwise HBF-DFT
    
    [SUPP_4bit,Gain_4bit,~] = SOLS_MMV(OV,A_4bit,P,M,Nsub,d,lambda,I);
    [SUPP_5bit,Gain_5bit,~] = SOLS_MMV(OV,A_5bit,P,M,Nsub,d,lambda,I);
    Gain_h=Gain_4bit/sqrt(trace((SUPP_4bit*Gain_4bit)'*(SUPP_4bit*Gain_4bit)))*sqrt(Ns);%power normalization
    F_ASP_DFT_4bit=SUPP_4bit*Gain_h;
    Gain_h=Gain_5bit/sqrt(trace((SUPP_5bit*Gain_5bit)'*(SUPP_5bit*Gain_5bit)))*sqrt(Ns);%power normalization
    F_ASP_DFT_5bit=SUPP_5bit*Gain_h;
    
    
    %% LC-AS-HBF-DFT
    
    [SUPP_4bit,Gain_4bit,~] = LC_SOLS_MMV(OV,A_4bit,P,M,Nsub,d,lambda);
    [SUPP_5bit,Gain_5bit,~] = LC_SOLS_MMV(OV,A_5bit,P,M,Nsub,d,lambda); 
    Gain_h=Gain_4bit/sqrt(trace((SUPP_4bit*Gain_4bit)'*(SUPP_4bit*Gain_4bit)))*sqrt(Ns);%power normalization
    F_LC_ASP_DFT_4bit=SUPP_4bit*Gain_h;
    Gain_h=Gain_5bit/sqrt(trace((SUPP_5bit*Gain_5bit)'*(SUPP_5bit*Gain_5bit)))*sqrt(Ns);%power normalization
    F_LC_ASP_DFT_5bit=SUPP_5bit*Gain_h;
    %% Angle-Stepwise HBF-SCHT
    
    [SUPP_4bit,Gain_4bit,~] = SOLS_MMV(OV,CS_SCHT_4bit,P,M,Nsub,d,lambda,I);
    [SUPP_5bit,Gain_5bit,~] = SOLS_MMV(OV,CS_SCHT_5bit,P,M,Nsub,d,lambda,I);
    Gain_h=Gain_4bit/sqrt(trace((SUPP_4bit*Gain_4bit)'*(SUPP_4bit*Gain_4bit)))*sqrt(Ns);%power normalization
    F_ASP_SCHT_4bit=SUPP_4bit*Gain_h;
    Gain_h=Gain_5bit/sqrt(trace((SUPP_5bit*Gain_5bit)'*(SUPP_5bit*Gain_5bit)))*sqrt(Ns);%power normalization
    F_ASP_SCHT_5bit=SUPP_5bit*Gain_h;
    
    %% LC-AS-HBF-DFT-SCHT
    
    [SUPP_4bit,Gain_4bit,~] = LC_SOLS_MMV(OV,CS_SCHT_4bit,P,M,Nsub,d,lambda);
    [SUPP_5bit,Gain_5bit,~] = LC_SOLS_MMV(OV,CS_SCHT_5bit,P,M,Nsub,d,lambda); 
    Gain_h=Gain_4bit/sqrt(trace((SUPP_4bit*Gain_4bit)'*(SUPP_4bit*Gain_4bit)))*sqrt(Ns);%power normalization
    F_LC_ASP_SCHT_4bit=SUPP_4bit*Gain_h;
    Gain_h=Gain_5bit/sqrt(trace((SUPP_5bit*Gain_5bit)'*(SUPP_5bit*Gain_5bit)))*sqrt(Ns);%power normalization
    F_LC_ASP_SCHT_5bit=SUPP_5bit*Gain_h;
    
    for i=1:length(SNR)
        snr=10^(SNR(i)/10); 
        SE_ASP_DFT_4bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ASP_DFT_4bit)'*H*F_ASP_DFT_4bit));
        SE_ASP_DFT_5bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ASP_DFT_5bit)'*H*F_ASP_DFT_5bit));
        SE_ASP_LC_DFT_4bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_LC_ASP_DFT_4bit)'*H*F_LC_ASP_DFT_4bit));
        SE_ASP_LC_DFT_5bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_LC_ASP_DFT_5bit)'*H*F_LC_ASP_DFT_5bit));
        
        SE_ASP_SCHT_4bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ASP_SCHT_4bit)'*H*F_ASP_SCHT_4bit));
        SE_ASP_SCHT_5bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ASP_SCHT_5bit)'*H*F_ASP_SCHT_5bit));
        SE_ASP_LC_SCHT_4bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_LC_ASP_SCHT_4bit)'*H*F_LC_ASP_SCHT_4bit));
        SE_ASP_LC_SCHT_5bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_LC_ASP_SCHT_5bit)'*H*F_LC_ASP_SCHT_5bit));
    end
end 
SE_ASP_DFT_4bit_mean=mean(SE_ASP_DFT_4bit);
SE_ASP_DFT_5bit_mean=mean(SE_ASP_DFT_5bit);
SE_ASP_LC_DFT_4bit_mean=mean(SE_ASP_LC_DFT_4bit);
SE_ASP_LC_DFT_5bit_mean=mean(SE_ASP_LC_DFT_5bit);
SE_ASP_SCHT_4bit_SCHT_mean=mean(SE_ASP_SCHT_4bit);
SE_ASP_SCHT_5bit_SCHT_mean=mean(SE_ASP_SCHT_5bit);
SE_ASP_LC_SCHT_4bit_mean=mean(SE_ASP_LC_SCHT_4bit);
SE_ASP_LC_SCHT_5bit_mean=mean(SE_ASP_LC_SCHT_5bit);
co1= [0, 161, 241]/255;
co2=[29, 191, 151]/255
co3= [70, 158, 180]/255
co4=[253,185,106]/255
co5=[214,64,78]/255
figure 
plot(SNR, abs(SE_ASP_DFT_4bit_mean), 'sk--', 'linewidth', 1.1, 'markerfacecolor', co5,'markersize', 7.5)
hold on
plot(SNR, abs(SE_ASP_DFT_5bit_mean), 'sk-', 'linewidth', 1.1, 'markerfacecolor', co5,'markersize', 7.5)
hold on
plot(SNR, abs(SE_ASP_SCHT_4bit_SCHT_mean), '^k--', 'linewidth', 1.1, 'markerfacecolor', co4,'markersize', 7.5)
hold on
plot(SNR, abs(SE_ASP_SCHT_5bit_SCHT_mean), '^k-', 'linewidth', 1.1, 'markerfacecolor', co4,'markersize', 7.5)
hold on
plot(SNR, abs(SE_ASP_LC_DFT_4bit_mean), 'dk--', 'linewidth', 1.1, 'markerfacecolor', co1,'markersize', 7.5)
hold on
plot(SNR, abs(SE_ASP_LC_DFT_5bit_mean), 'dk-', 'linewidth', 1.1, 'markerfacecolor', co1,'markersize', 7.5)
hold on
plot(SNR, abs(SE_ASP_LC_SCHT_4bit_mean), '>k--', 'linewidth', 1.1, 'markerfacecolor', co3,'markersize', 7.5)
hold on
plot(SNR, abs(SE_ASP_LC_SCHT_5bit_mean), '>k-', 'linewidth', 1.1, 'markerfacecolor', co3,'markersize', 7.5)

grid on
lgh=legend('SJG-HBF-I, $Q=16$, $4$-bit PS','SJG-HBF-I, $Q=32$, $5$-bit PS',...
    'SJG-HBF-II, $Q=16$, $2$-bit PS','SJG-HBF-II, $Q=32$, $2$-bit PS',...
    'SIG-HBF-I, $Q=16$, $4$-bit PS','SIG-HBF-I, $Q=32$, $5$-bit PS',...
    'SIG-HBF-II, $Q=16$, $2$-bit PS','SIG-HBF-II, $Q=32$, $2$-bit PS');
set(lgh,'interpreter','latex', 'fontsize', 14);
xlabel('SNR [dB]','interpreter','latex','fontsize',14)
ylabel('Spectral Efficiency [bit/s/Hz]','interpreter','latex','fontsize',14)