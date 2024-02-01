function ULA_Dic = ULA_NFDic(d,lambda,N,theta_sample,r_sample)
 
ULA_Dic=[];

for theta=theta_sample
    for r=r_sample
        r_s=[];G_s=[];
        G_s=SW(theta,r,d,lambda,N);
        ULA_Dic=[ULA_Dic G_s];
    end
end
end