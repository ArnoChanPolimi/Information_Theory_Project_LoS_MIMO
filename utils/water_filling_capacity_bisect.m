%%% Water Filling %%%%%%
% function [C, p_opt, sigma2, mu, Ka] = water_filling_capacity_bisect(H, energy, P, N0, B)
function [C, p_opt, sigma2, mu, Ka] = water_filling_capacity_bisect(H, P, N0, B)
    % H = [1 0; 0 1]; % de-bug用途
    Nt = 2;

    % [~, S, ~] = svd(H, 'econ');
    [~, S, ~] = svd(H);
    sigma = diag(S);            % 信道奇异值
    % disp(sigma);
    Ka = max(sigma)/min(sigma); % 条件数，越小越好，最小是1
    % disp(Ka);
    sigma2 = sigma.^2;          % 奇异值平方，对应信道增益（特征值）

    n = length(sigma2);
    
    % 关键修正：计算实际噪声功率
    noise_power = N0 * B;  % 噪声功率 = 功率谱密度 × 带宽
    % noise_power = noise_power/energy; % 对噪声也进行归一化 难崩！
    N_eq = noise_power ./ sigma2;   % 等效噪声：N0*B / |h_i|^2
    
    % 修正求和函数：使用噪声功率而非功率谱密度
    power_sum = @(mu) sum(max(0, mu - N_eq));
    
    % 二分法参数设置（修正边界）
    lower_bound = 0;                % 允许关闭弱信道
    upper_bound = P/2 + max(N_eq);    % 合理上界
    tol = 1e-12;
    max_iter = 2000;
    mu = (lower_bound + upper_bound) / 2;
    
    % 二分法求最优注水水平 mu
    for iter = 1:max_iter
        total_power = power_sum(mu);
        
        if abs(total_power - P) < tol
            break;
        elseif total_power > P
            upper_bound = mu;
        else
            lower_bound = mu;
        end
        mu = (lower_bound + upper_bound) / 2;
    end

    % 计算最优功率分配
    p_opt = max(0, mu - N_eq); % p_i = max(0, mu - N_eq) with N_eq = B*N0 / sigma_i^2
    
    epsilon = 1e-4;  % 误差容忍范围，可以根据实际调整大小
    if p_opt(1) < 0 || p_opt(2) < 0
        fprintf('p_opt(1) = %.4f, p_opt(2) = %.4f < 0\n', p_opt(1), p_opt(2));
    end
    if abs(p_opt(1) + p_opt(2) - P) > epsilon
        fprintf('Warning: WF distribution ERROR!!!! ---P is %.4fW, but p_opt(1) = %.4f, p_opt(2) = %.4f\n', P, p_opt(1), p_opt(2));
    end

    % 计算总容量（修正带宽位置）
    SNR = p_opt .* sigma2 / noise_power; % 各子信道信噪比
    % C = 0;
    % for i = 1 : 2
    %     C = C + B * log2(1 + p_opt(i)*sigma2(i)/noise_power);
    % end
    
    C = B * sum(log2(1 + SNR));         % 总容量 (bit/s)

    %% debug:
    
    epsilon_sigma2 = 1e-19;
    if abs(sigma2(1)-sigma2(2)) < epsilon_sigma2
        fprintf('Now, the p_opt(1) = %.4f W, p_opt(2) = %.4f W\n', p_opt(1), p_opt(2));
        fprintf('!!!equal!!!---sigma2(1)= %.2ef, fsigma2(2)=%.2ef, C_WF = %.3e (bits/s)\n', sigma2(1), sigma2(2), C);
    end
    disp('sigma = ')
    disp(sigma);
    disp('奇异值 sigma2:');
    disp(sigma2);
    disp('等效噪声 N_eq:');
    disp(N_eq);
    disp('水位线 mu:');
    disp(mu);
    disp('功率分配 p_opt:');
    disp(p_opt);
    disp('总功率和:');
    disp(sum(p_opt));

end



    % % 二分法求mu
    % for iter = 1:max_iter
    %     mu = (lower + upper)/2;
    %     total_power = power_sum(mu);
    % 
    %     if abs(total_power - P) < tol
    %         break;
    %     elseif total_power > P
    %         upper = mu;
    %     else
    %         lower = mu;
    %     end
    % end