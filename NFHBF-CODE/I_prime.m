function I_prime = I_prime(N)
%I_PRIME �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
I_prime=blkdiag(eye(N/2),eye(N/4),-eye(N/4));
end

