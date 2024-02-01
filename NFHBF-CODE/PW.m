function g_s = PW(phi,d,lambda,N)  
 
    for n=1:N  
        G_s(n)=1/sqrt(N)*exp(-1i*2*pi/lambda*(-(n-1)*d*phi)); 
    end
 
g_s=G_s(:);
 

end

