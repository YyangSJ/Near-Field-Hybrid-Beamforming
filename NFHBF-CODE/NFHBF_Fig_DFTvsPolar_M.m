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
MM=4:1:10;
L=4;
Ns=4;
R_sample=0.01:0.01:0.2;
r_sample=1./R_sample;
Realization=200
SNR=10;
I=8; 
vn=5;

for reali=1:Realization
    reali
     H=Channel_realization(d,lambda,L,Nr,Nt);
    %H=Channel_realization_non_stationary(d,lambda,L,Nr,Nt,P,vn);
    [~,Sig,V]=svd(H);
    for i=1:length(MM)
        M=MM(i);
        snr=10^(SNR/10);
        
        %% Optimal Fully-Digital BF
        
        OV=V(:,1:Ns);
        OV=OV/sqrt(trace(OV'*OV))*sqrt(Ns);  %power normalization
        SE_Opt(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*OV)'*H*OV));
        
        %% Angle-Distance HBF
        G_theta=Nt/8;
        theta_sample_4=-1+2/G_theta:2/G_theta:1;
        OFRF_4=ULA_NFDic(d,lambda,Nt,theta_sample_4,r_sample);
        G_theta=Nt/2;
        theta_sample_6=-1+2/G_theta:2/G_theta:1;
        OFRF_6=ULA_NFDic(d,lambda,Nt,theta_sample_6,r_sample);
        
        [Lambda,GAIN]=OLS_MMV(OV,OFRF_4,M);
        FRF_EST_4=OFRF_4(:,Lambda);
        FBB_EST_4=GAIN;
        FBB_EST_4=FBB_EST_4/sqrt(trace((FRF_EST_4*FBB_EST_4)'*(FRF_EST_4*FBB_EST_4)))*sqrt(Ns); %power normalization
        F_ADP_4=FRF_EST_4*FBB_EST_4;
        
        [Lambda,GAIN]=OLS_MMV(OV,OFRF_6,M);
        FRF_EST_6=OFRF_6(:,Lambda);
        FBB_EST_6=GAIN;
        FBB_EST_6=FBB_EST_6/sqrt(trace((FRF_EST_6*FBB_EST_6)'*(FRF_EST_6*FBB_EST_6)))*sqrt(Ns); %power normalization
        F_ADP_6=FRF_EST_6*FBB_EST_6;
        
        SE_ADP_4(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ADP_4)'*H*F_ADP_4));
        SE_ADP_6(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ADP_6)'*H*F_ADP_6));
        
        %% Angle-Stepwise HBF
        A_4bit=ULA_FFDic(d,lambda,Nsub,-1+2/(2^4):2/(2^4):1);
        
        A_6bit=ULA_FFDic(d,lambda,Nsub,-1+2/(2^6):2/(2^6):1);
 
        [SUPP_4bit,Gain_4bit,~] = SOLS_MMV(OV,A_4bit,P,M,Nsub,d,lambda,I); 
        [SUPP_6bit,Gain_6bit,~] = SOLS_MMV(OV,A_6bit,P,M,Nsub,d,lambda,I); 
        
        
        [SUPP_OMP_4bit,Gain_omp_4bit,~] = SOMP_MMV(OV,A_4bit,P,M,Nsub,d,lambda,I);
        [SUPP_OMP_6bit,Gain_omp_6bit,~] = SOMP_MMV(OV,A_6bit,P,M,Nsub,d,lambda,I);
        
        
        
        Gain_h=Gain_4bit/sqrt(trace((SUPP_4bit*Gain_4bit)'*(SUPP_4bit*Gain_4bit)))*sqrt(Ns);%power normalization
        F_ASP_4bit=SUPP_4bit*Gain_h;
        Gain_h=Gain_6bit/sqrt(trace((SUPP_6bit*Gain_6bit)'*(SUPP_6bit*Gain_6bit)))*sqrt(Ns);%power normalization
        F_ASP_6bit=SUPP_6bit*Gain_h;
        Gain_h=Gain_omp_4bit/sqrt(trace((SUPP_OMP_4bit*Gain_omp_4bit)'*(SUPP_OMP_4bit*Gain_omp_4bit)))*sqrt(Ns);%power normalization
        F_ASP_omp_4bit=SUPP_OMP_4bit*Gain_h;
        Gain_h=Gain_omp_6bit/sqrt(trace((SUPP_OMP_6bit*Gain_omp_6bit)'*(SUPP_OMP_6bit*Gain_omp_6bit)))*sqrt(Ns);%power normalization
        F_ASP_omp_6bit=SUPP_OMP_6bit*Gain_h;
        SE_ASP_4bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ASP_4bit)'*H*F_ASP_4bit));
        SE_ASP_6bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ASP_6bit)'*H*F_ASP_6bit));
        SE_ASP_omp_5bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ASP_omp_4bit)'*H*F_ASP_omp_4bit));
        SE_ASP_omp_7bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ASP_omp_6bit)'*H*F_ASP_omp_6bit));
        
        %% LC-AS-HBF 
        
        [SUPP_4bit,Gain_4bit,~] = LC_SOLS_MMV(OV,A_4bit,P,M,Nsub,d,lambda); 
        [SUPP_6bit,Gain_6bit,~] = LC_SOLS_MMV(OV,A_6bit,P,M,Nsub,d,lambda);
        
        Gain_h=Gain_4bit/sqrt(trace((SUPP_4bit*Gain_4bit)'*(SUPP_4bit*Gain_4bit)))*sqrt(Ns);%power normalization
        F_ASP_4bit=SUPP_4bit*Gain_h;
        Gain_h=Gain_6bit/sqrt(trace((SUPP_6bit*Gain_6bit)'*(SUPP_6bit*Gain_6bit)))*sqrt(Ns);%power normalization
        F_ASP_6bit=SUPP_6bit*Gain_h;
        SE_LCASP_4bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ASP_4bit)'*H*F_ASP_4bit));
        SE_LCASP_6bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ASP_6bit)'*H*F_ASP_6bit));
    end
end
SE_Opt_mean=mean(SE_Opt);
SE_ADP_4_mean=mean(SE_ADP_4);
SE_ADP_6_mean=mean(SE_ADP_6);
SE_ASP_4bit_mean=mean(SE_ASP_4bit);
SE_ASP_6bit_mean=mean(SE_ASP_6bit);
SE_LCASP_4bit_mean=mean(SE_LCASP_4bit);
SE_LCASP_6bit_mean=mean(SE_LCASP_6bit); 
co1= [0, 161, 241]/255;
co2=[29, 191, 151]/255
co3= [70, 158, 180]/255
co4=[253,185,106]/255
co5=[214,64,78]/255
figure
plot(MM, abs(SE_Opt_mean), 'ok-', 'linewidth', 1.1, 'markerfacecolor', co2,'markersize', 7.5)
hold on
plot(MM, abs(SE_ADP_4_mean), '^k--', 'linewidth', 1.1, 'markerfacecolor', co4,'markersize', 7.5)
hold on
plot(MM, abs(SE_ADP_6_mean), '^k-', 'linewidth', 1.1, 'markerfacecolor', co4,'markersize', 7.5)
hold on
plot(MM, abs(SE_ASP_4bit_mean), 'sk--', 'linewidth', 1.1, 'markerfacecolor', co5,'markersize', 7.5)
hold on
plot(MM, abs(SE_ASP_6bit_mean), 'sk-', 'linewidth', 1.1, 'markerfacecolor', co5,'markersize', 7.5)
hold on
plot(MM, abs(SE_LCASP_4bit_mean), 'dk--', 'linewidth', 1.1, 'markerfacecolor', co1,'markersize', 7.5)
hold on
plot(MM, abs(SE_LCASP_6bit_mean), 'dk-', 'linewidth', 1.1, 'markerfacecolor', co1,'markersize', 7.5) 

grid on
 lgh=legend('Fully-Digital Beamforming','DG-HBF, $Q=16$, $16$-bit PS','DG-HBF, $Q=64$, $18$-bit PS'...
     ,'SJG-HBF-I, $Q=16$, $4$-bit PS','SJG-HBF-I, $Q=64$, $6$-bit PS','SIG-HBF-I, $Q=16$, $4$-bit PS','SIG-HBF-I, $Q=64$, $6$-bit PS');
set(lgh,'interpreter','latex', 'fontsize', 14); 
xlabel('Number of RF chians $M$','interpreter','latex','fontsize',14)
ylabel('Spectral Efficiency [bit/s/Hz]','interpreter','latex','fontsize',14) 
 axis([4,10,10,38])