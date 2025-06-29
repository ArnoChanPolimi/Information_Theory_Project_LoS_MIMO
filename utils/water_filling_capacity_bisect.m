%%% Water Filling %%%%%%
function [C, p_opt, sigma2, mu] = water_filling_capacity_bisect(H, P, N0, B)
    [~, S, ~] = svd(H);
    sigma = diag(S);            % 信道奇异值
    sigma2 = sigma.^2;          % 奇异值平方，对应信道增益（特征值）
    n = length(sigma2);
    
    % 关键修正：计算实际噪声功率
    noise_power = N0 * B;  % 噪声功率 = 功率谱密度 × 带宽
    
    % 修正求和函数：使用噪声功率而非功率谱密度
    power_sum = @(mu) sum(max(0, mu - noise_power ./ sigma2));
    
    % 二分法参数（修正边界计算）
    lower = min(noise_power ./ sigma2);
    upper = P + max(noise_power ./ sigma2);
    tol = 1e-10;
    max_iter = 2000;
    mu = (lower + upper)/2;  % 初始值
    
    % 二分法求mu
    for iter = 1:max_iter
        mu = (lower + upper)/2;
        total_power = power_sum(mu);
        
        if abs(total_power - P) < tol
            break;
        elseif total_power > P
            upper = mu;
        else
            lower = mu;
        end
    end
    
    % 计算功率分配（使用噪声功率）
    p_opt = max(0, mu - noise_power ./ sigma2);
    
    % 修正容量计算公式（使用正确的信噪比）
    C = sum(B * log2(1 + p_opt .* sigma2 / noise_power)); % 输出 C 单位为 bit/s
end
