function ULA_Dic = ULA_FFDic(d,lambda,N,theta_sample)
ULA_Dic=[];

for theta=theta_sample 
        G_s=PW(theta,d,lambda,N);
        ULA_Dic=[ULA_Dic G_s]; 
end
end

