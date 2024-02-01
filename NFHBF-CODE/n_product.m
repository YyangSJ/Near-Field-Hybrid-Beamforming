function WI = n_product(I,W,n)
%N_PRODUCT 此处显示有关此函数的摘要
%   此处显示详细说明
[N1,N2]=size(W);
WI=zeros(N1,N2);
for i=1:n:N1
    for j=1:n:N2
        WI(i:i+n-1,j:j+n-1)=W(i:i+n-1,j:j+n-1)*I(i:i+n-1,i:i+n-1);
    end
end
end

