function W_expanded = expandMatrix(W, n)
    % 递归函数用于扩展矩阵 W
    % W - 初始矩阵
    % n - 递推次数

    if n == 1
        % 基本情况，返回初始矩阵
        W_expanded = W;
    else
        % 递归调用，扩展矩阵
        W_prev = expandMatrix(W, n - 1);
        W_expanded = [W_prev, W_prev; W_prev, -W_prev];
    end
end