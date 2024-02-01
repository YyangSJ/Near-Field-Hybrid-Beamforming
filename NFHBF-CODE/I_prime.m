function I_prime = I_prime(N)
%I_PRIME 此处显示有关此函数的摘要
%   此处显示详细说明
I_prime=blkdiag(eye(N/2),eye(N/4),-eye(N/4));
end

