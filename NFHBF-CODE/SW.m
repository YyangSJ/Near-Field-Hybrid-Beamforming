function g_s = SW(theta,r,d,lambda,N)
 
    for n=1:N
        r_s(n)=sqrt(r^2+((n-1)*d)^2-2*r*(n-1)*d*theta);   
        G_s(n)=1/sqrt(N)*exp(-1i*2*pi/lambda*(r_s(n)-r)); 
    end
 
g_s=G_s(:);
end

