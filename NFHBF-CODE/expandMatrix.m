function W_expanded = expandMatrix(W, n)
    % �ݹ麯��������չ���� W
    % W - ��ʼ����
    % n - ���ƴ���

    if n == 1
        % ������������س�ʼ����
        W_expanded = W;
    else
        % �ݹ���ã���չ����
        W_prev = expandMatrix(W, n - 1);
        W_expanded = [W_prev, W_prev; W_prev, -W_prev];
    end
end