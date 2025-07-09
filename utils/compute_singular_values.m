%%% 计算H矩阵的奇异值 %%%
function [sigma, sigma2] = compute_singular_values(H)
    % compute_singular_values 计算给定矩阵H的奇异值及其平方
    %
    % 输入:
    %   H - 信道矩阵 (复数矩阵)
    %
    % 输出:
    %   sigma  - 奇异值向量，降序排列
    %   sigma2 - 奇异值平方向量

    % 计算奇异值
    [~, S, ~] = svd(H, 'econ'); % Debug: compute singular values of H
    sigma = diag(S);
    
    % 计算奇异值平方
    sigma2 = sigma.^2;
    
end
