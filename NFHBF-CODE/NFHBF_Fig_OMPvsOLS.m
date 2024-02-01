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
Realization=2
SNR=-30:5:20;
I=8  ;
vn=5;

A_4bit=ULA_FFDic(d,lambda,Nsub,-1+2/(2^4):2/(2^4):1);
A_6bit=ULA_FFDic(d,lambda,Nsub,-1+2/(2^6):2/(2^6):1);
for reali=1:Realization
    reali
    H=Channel_realization(d,lambda,L,Nr,Nt);
    %H=Channel_realization_non_stationary(d,lambda,L,Nr,Nt,P,vn);
    
    [~,Sig,V]=svd(H);
    OV=V(:,1:Ns);
    OV=OV/sqrt(trace(OV'*OV))*sqrt(Ns);  %power normalization
    %% Angle-Stepwise HBF-DFT
    G_theta=Nt/8;
    theta_sample_4=-1+2/G_theta:2/G_theta:1;
    OFRF_4=ULA_NFDic(d,lambda,Nt,theta_sample_4,r_sample);
    G_theta=Nt/2;
    theta_sample_6=-1+2/G_theta:2/G_theta:1;
    OFRF_6=ULA_NFDic(d,lambda,Nt,theta_sample_6,r_sample);
%     [Lambda,GAIN]=OLS_MMV(OV,OFRF_4,M);
%     FRF_EST_OLS_4=OFRF_4(:,Lambda);
%     FBB_EST_OLS_4=GAIN;
%     FBB_EST_OLS_4=FBB_EST_OLS_4/sqrt(trace((FRF_EST_OLS_4*FBB_EST_OLS_4)'*(FRF_EST_OLS_4*FBB_EST_OLS_4)))*sqrt(Ns); %power normalization
%     F_ADP_OLS_4=FRF_EST_OLS_4*FBB_EST_OLS_4;
    [Lambda,GAIN]=OLS_MMV(OV,OFRF_6,M);
    FRF_EST_OLS_6=OFRF_6(:,Lambda);
    FBB_EST_OLS_6=GAIN;
    FBB_EST_OLS_6=FBB_EST_OLS_6/sqrt(trace((FRF_EST_OLS_6*FBB_EST_OLS_6)'*(FRF_EST_OLS_6*FBB_EST_OLS_6)))*sqrt(Ns); %power normalization
    F_ADP_OLS_6=FRF_EST_OLS_6*FBB_EST_OLS_6;
    
%     [SUPP_OMP_4bit,Gain_omp_4bit] = cs_somp(OV,OFRF_4,size(OFRF_4,2),M) ;
%     FRF_EST_OMP_4=OFRF_4(:,SUPP_OMP_4bit);
%     FBB_EST_OMP_4=Gain_omp_4bit(SUPP_OMP_4bit,:);
%     FBB_EST_OMP_4=FBB_EST_OMP_4/sqrt(trace((FRF_EST_OMP_4*FBB_EST_OMP_4)'*(FRF_EST_OMP_4*FBB_EST_OMP_4)))*sqrt(Ns); %power normalization
%     F_ADP_OMP_4=FRF_EST_OMP_4*FBB_EST_OMP_4;
    [SUPP_OMP_6bit,Gain_omp_6bit] = cs_somp(OV,OFRF_6,size(OFRF_6,2),M) ;
    FRF_EST_OMP_6=OFRF_6(:,SUPP_OMP_6bit);
    FBB_EST_OMP_6=Gain_omp_6bit(SUPP_OMP_6bit,:);
    FBB_EST_OMP_6=FBB_EST_OMP_6/sqrt(trace((FRF_EST_OMP_6*FBB_EST_OMP_6)'*(FRF_EST_OMP_6*FBB_EST_OMP_6)))*sqrt(Ns); %power normalization
    F_ADP_OMP_6=FRF_EST_OMP_6*FBB_EST_OMP_6;
    
    
%     [SUPP_4bit,Gain_4bit,~] = SOLS_MMV(OV,A_4bit,P,M,Nsub,d,lambda,I);
%     FBB_EST_SJGOLS_4=Gain_4bit/sqrt(trace((SUPP_4bit*Gain_4bit)'*(SUPP_4bit*Gain_4bit)))*sqrt(Ns); %power normalization
%     F_SJG_OLS_4=SUPP_4bit*FBB_EST_SJGOLS_4;
%     [SUPP2,Gain2,~] = SOMP_MMV(OV,A_4bit,P,M,Nsub,d,lambda);
%     FBB_EST_SJGOMP_4=Gain2/sqrt(trace((SUPP2*Gain2)'*(SUPP2*Gain2)))*sqrt(Ns); %power normalization
%     F_SJG_OMP_4=SUPP2*FBB_EST_SJGOMP_4;
    [SUPP_6bit,Gain_6bit,~] = SOLS_MMV(OV,A_6bit,P,M,Nsub,d,lambda,I);
    FBB_EST_SJGOLS_6=Gain_6bit/sqrt(trace((SUPP_6bit*Gain_6bit)'*(SUPP_6bit*Gain_6bit)))*sqrt(Ns); %power normalization
    F_SJG_OLS_6=SUPP_6bit*FBB_EST_SJGOLS_6;
    [SUPP2,Gain2,~] = SOMP_MMV(OV,A_6bit,P,M,Nsub,d,lambda,I);
    FBB_EST_SJGOMP_6=Gain2/sqrt(trace((SUPP2*Gain2)'*(SUPP2*Gain2)))*sqrt(Ns); %power normalization
    F_SJG_OMP_6=SUPP2*FBB_EST_SJGOMP_6;
    for i=1:length(SNR)
        snr=10^(SNR(i)/10);
        
        SE_Opt(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*OV)'*H*OV));
%                SE_ADP_OLS_4bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ADP_OLS_4)'*H*F_ADP_OLS_4));
%         SE_ADP_OMP_4bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ADP_OMP_4)'*H*F_ADP_OMP_4));
%         SE_SJG_OMP_4bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_SJG_OMP_4)'*H*F_SJG_OMP_4)); 
%         SE_SJG_OLS_4bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_SJG_OLS_4)'*H*F_SJG_OLS_4));
        SE_ADP_OLS_6bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ADP_OLS_6)'*H*F_ADP_OLS_6));
        SE_ADP_OMP_6bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_ADP_OMP_6)'*H*F_ADP_OMP_6));
 
        SE_SJG_OMP_6bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_SJG_OMP_6)'*H*F_SJG_OMP_6));
        SE_SJG_OLS_6bit(reali,i)=log2(det(eye(Ns)+1/Ns*snr*(H*F_SJG_OLS_6)'*H*F_SJG_OLS_6));
    end
end
SE_Opt_mean=mean(SE_Opt);

% SE_ADP_OLS_4bit_mean=mean(SE_ADP_OLS_4bit);
% SE_SJG_OLS_4bit_mean=mean(SE_SJG_OLS_4bit);
% SE_ADP_OMP_4bit_mean=mean(SE_ADP_OMP_4bit);
% SE_SJG_OMP_4bit_mean=mean(SE_SJG_OMP_4bit);

SE_ADP_OLS_6bit_mean=mean(SE_ADP_OLS_6bit);
SE_SJG_OLS_6bit_mean=mean(SE_SJG_OLS_6bit);
SE_ADP_OMP_6bit_mean=mean(SE_ADP_OMP_6bit);
SE_SJG_OMP_6bit_mean=mean(SE_SJG_OMP_6bit);
co1= [0, 161, 241]/255;
co2=[29, 191, 151]/255
co3= [70, 158, 180]/255
co4=[253,185,106]/255
co5=[214,64,78]/255
figure
plot(SNR, abs(SE_Opt_mean), 'ok-', 'linewidth', 1.1, 'markerfacecolor', co2,'markersize', 7.5)
hold on
% plot(SNR, abs(SE_ADP_OLS_4bit_mean), 'sk--', 'linewidth', 1.1, 'markerfacecolor', co5,'markersize', 7.5)
% hold on
plot(SNR, abs(SE_ADP_OLS_6bit_mean), 'sk-', 'linewidth', 1.1, 'markerfacecolor', co5,'markersize', 7.5)
hold on
% plot(SNR, abs(SE_ADP_OMP_4bit_mean), '^k--', 'linewidth', 1.1, 'markerfacecolor', co4,'markersize', 7.5)
% hold on
plot(SNR, abs(SE_ADP_OMP_6bit_mean), '^k-', 'linewidth', 1.1, 'markerfacecolor', co4,'markersize', 7.5)
hold on
% plot(SNR, abs(SE_SJG_OLS_4bit_mean), '<k--', 'linewidth', 1.1, 'markerfacecolor', co1,'markersize', 7.5)
% hold on
plot(SNR, abs(SE_SJG_OLS_6bit_mean), '<k-', 'linewidth', 1.1, 'markerfacecolor', co1,'markersize', 7.5)
hold on
% plot(SNR, abs(SE_SJG_OMP_4bit_mean), '>k--', 'linewidth', 1.1, 'markerfacecolor', co3,'markersize', 7.5)
% hold on
plot(SNR, abs(SE_SJG_OMP_6bit_mean), '>k-', 'linewidth', 1.1, 'markerfacecolor', co3,'markersize', 7.5)


grid on
lgh=legend('Fully-Digital Beamforming','DG-HBF, OLS, $Q=64$',...
    'DG-HBF, OMP, $Q=64$','SJG-HBF, OLS, $Q=64$',...
    'SJG-HBF, OMP, $Q=64$');
% lgh=legend('Fully-Digital Beamforming','DG-HBF, OLS, $Q=16$','DG-HBF, OLS, $Q=64$',...
%     'DG-HBF, OMP $Q=16$','DG-HBF, OMP, $Q=64$','SJG-HBF, OLS, $Q=16$','SJG-HBF, OLS, $Q=64$',...
%     'SJG-HBF, OMP $Q=16$','SJG-HBF, OMP, $Q=64$');
set(lgh,'interpreter','latex', 'fontsize', 14);
xlabel('SNR [dB]','interpreter','latex','fontsize',14)
ylabel('Spectral Efficiency [bit/s/Hz]','interpreter','latex','fontsize',14)