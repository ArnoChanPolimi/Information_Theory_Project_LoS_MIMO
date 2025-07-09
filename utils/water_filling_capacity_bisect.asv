%%% Water Filling %%%%%%
% function [C, p_opt, sigma2, mu, Ka] = water_filling_capacity_bisect(H, energy, P, N0, B)
% function [C, p_opt, sigma2_sorted, mu, Ka] = water_filling_capacity_bisect(H, P, N0, B)
%     % H = [1 0; 0 1]; % de-bug用途
%     Nt = 2;
%     % [~, S, ~] = svd(H, 'econ');
%     [~, S, ~] = svd(H);
%     sigma = diag(S);            % 信道奇异值
%     Ka = max(sigma)/min(sigma); % 条件数，越小越好，最小是1
%     sigma2 = sigma.^2;          % 奇异值平方，对应信道增益（特征值）
%     [sigma2_sorted, idx] = sort(sigma2, 'descend');  % 降序排列
%     for i = 1:2
%         fprintf('WF: sigma2_sorted(%d) = %.3e \n', i, sigma2_sorted(i));
%     end
% 
%     % 关键修正：计算实际噪声功率
%     noise_power = N0 * B;  % 噪声功率 = 功率谱密度 × 带宽
%     % noise_power = noise_power/energy; % 对噪声也进行归一化 难崩！
%     % N_eq = noise_power ./ sigma2_sorted;   % 等效噪声：N0*B / |h_i|^2
%     N_eq = N0 ./ sigma2_sorted;   % 等效噪声：N0*B / |h_i|^2
%     % 修正求和函数：使用噪声功率而非功率谱密度
%     power_sum = @(mu) sum(max(0, mu - N_eq));
% 
%     % 二分法参数设置（修正边界）
%     lower_bound = 0;                % 允许关闭弱信道
%     upper_bound = P/B + max(N_eq);    % 合理上界
%     tol = 1e-22;
%     max_iter = 10000;
%     mu = (lower_bound + upper_bound) / 2;
% 
%     % 二分法求最优注水水平 mu
%     for iter = 1:max_iter
%         total_power = power_sum(mu);
% 
%         if abs(total_power - P) < tol
%             break;
%         elseif total_power > P
%             upper_bound = mu;
%         else
%             lower_bound = mu;
%         end
%         mu = (lower_bound + upper_bound) / 2;
%     end
% 
% 
%     % 计算最优功率分配
%     p_opt = max(0, mu - N_eq); % p_i = max(0, mu - N_eq) with N_eq = B*N0 / sigma_i^2
%     % p_opt = B * p_opt_density;
%     for i = 1:2
%         fprintf('p_opt = %.4f\n', p_opt(i));
%     end    
%     epsilon = 1e-4;  % 误差容忍范围，可以根据实际调整大小
%     % if p_opt(1) < 0 || p_opt(2) < 0
%     %     fprintf('!!! ERROR: p_opt(1) = %.4f, p_opt(2) = %.4f < 0\n', p_opt(1), p_opt(2));
%     % end
%     % if abs(p_opt(1) + p_opt(2) - P) > epsilon
%     %     fprintf('Warning: WF distribution ERROR!!!! ---P is %.4fW, but p_opt(1) = %.4f, p_opt(2) = %.4f\n', P, p_opt(1), p_opt(2));
%     % end
%     % 
%     % if abs(p_opt(1) - p_opt(2)) < epsilon
%     %     fprintf('Attention - p_opt: WF distribution is equal!!!! ---P is %.4fW, and p_opt(1) = %.4f, p_opt(2) = %.4f\n', P, p_opt(1), p_opt(2));
%     % end
%     % if abs(sigma2_sorted(1) - sigma2_sorted(2)) < (epsilon * 1e-10)
%     %     fprintf('Attention - sigma2_sorted: WF distribution is equal!!!! --- sigma2_sorted(1) = %.4e, sigma2_sorted(2) = %.4e\n', sigma2_sorted(1), sigma2_sorted(2));
%     % end
% 
%     % 计算总容量（修正带宽位置）
%     % p_opt(1) = P/2; % debug: 让WF的p_opt强制均分
%     % p_opt(2) = P/2; % debug: 让WF的p_opt强制均分
%     SNR_HH = p_opt .* sigma2_sorted / noise_power; % 各子信道信噪比
%     % for i = 1:2
%     %     fprintf('SNR_HH_WF(%d) = %.4f\n', i, SNR_HH(i));
%     % end
%     % C = 0;
%     % for i = 1 : 2
%     %     C = C + B * log2(1 + p_opt(i)*sigma2(i)/noise_power);
%     % end
% 
%     % C = B * sum(log2(1 + SNR_HH));         % 总容量 (bit/s)
%     C = B * sum(log2(1+ p_opt .* sigma2_sorted/ (N0 * B)));
% end
% 
%     %% debug:
% 
%     % 
%     % disp('sigma = ')
%     % disp(sigma);
%     % disp('奇异值 sigma2:');
%     % disp(sigma2_sorted);
%     % disp('等效噪声 N_eq:');
%     % disp(N_eq);
%     % disp('水位线 mu:');
%     % disp(mu);
%     % disp('功率分配 p_opt:');
%     % disp(p_opt);
%     % disp('总功率和:');
%     % disp(sum(p_opt));

% function [C, p_opt, sigma2_sorted, mu, Ka] = water_filling_capacity_bisect(H, P, N0, B)
%     % 确保输入有效
%     if nargin < 4, error('缺少参数'); end
%     if P <= 0, error('功率P必须为正数'); end
% 
%     % SVD分解获取信道增益
%     [~, S, ~] = svd(H, 'econ');
%     sigma = diag(S);
%     sigma2 = sigma.^2;
% 
%     % 移除接近零的奇异值 (避免除零错误)
%     min_sigma2 = max(sigma2) * 1e-10;  % 设置合理阈值
%     sigma2(sigma2 < min_sigma2) = min_sigma2;
% 
%     [sigma2_sorted, idx] = sort(sigma2, 'descend');
%     N = length(sigma2_sorted);
%     Ka = max(sigma2_sorted) / min(sigma2_sorted);  % 条件数
% 
%     % 核心：计算等效噪声功率
%     noise_power_per_channel = N0 * B;  % 实际噪声功率(W)
%     N_eq = noise_power_per_channel ./ sigma2_sorted;  % 等效噪声
% 
%     % 注水法核心函数
%     power_sum = @(mu) sum(max(0, mu - N_eq));
% 
%     % 稳健的二分法设置
%     lower_bound = min(N_eq);  % 最低注水线
%     upper_bound = P + 1000*max(N_eq);  % 最高可能注水线
%     tol = 1e-10;  % 合理容差
%     max_iter = 1000;
%     mu = (lower_bound + upper_bound) / 2;
% 
%     % 带保护的二分法
%     for iter = 1:max_iter
%         total_power = power_sum(mu);
% 
%         if abs(total_power - P) < tol
%             break;
%         elseif total_power > P
%             upper_bound = mu;
%         else
%             lower_bound = mu;
%         end
% 
%         % 更新并检查边界
%         mu = (lower_bound + upper_bound) / 2;
%         if (upper_bound - lower_bound) < tol
%             break;
%         end
%     end
%     disp('mu:')
%     disp(mu);
%     disp('N_eq:')
%     disp(N_eq);
%     % 最优功率分配
%     p_opt = max(0, mu - N_eq);
% 
%     % % 功率归一化 (确保总功率精确)
%     % p_opt = p_opt * P / sum(p_opt);
% 
%     % 正确的容量计算
%     C = 0;
%     for i = 1:N
%         if p_opt(i) > 0
%             SNR = (p_opt(i) * sigma2_sorted(i)) / noise_power_per_channel;
%             C = C + B * log2(1 + SNR);
%         end
%     end
% end

% function [C, p_opt, sigma2_sorted, mu, Ka] = water_filling_capacity_bisect(H, P, N0, B)
%     % 确保输入有效
%     if nargin < 4, error('缺少参数'); end
%     if P <= 0, error('功率P必须为正数'); end
% 
%     % SVD分解获取信道增益 - 使用'econ'模式
%     [~, S, V] = svd(H, 'econ');  % 添加V矩阵用于预编码
%     sigma = diag(S);
%     sigma2 = sigma.^2;
% 
%     % 移除接近零的奇异值 (避免除零错误)
%     min_sigma2 = max(sigma2) * 1e-10;  % 设置合理阈值
%     sigma2(sigma2 < min_sigma2) = min_sigma2;
% 
%     [sigma2_sorted, idx] = sort(sigma2, 'descend');
%     N = length(sigma2_sorted);
%     Ka = max(sigma2_sorted) / min(sigma2_sorted);  % 条件数
% 
%     % 核心：计算等效噪声功率
%     noise_power_per_channel = N0 * B;  % 实际噪声功率(W)
%     N_eq = noise_power_per_channel ./ sigma2_sorted;  % 等效噪声
% 
%     % 注水法核心函数
%     power_sum = @(mu) sum(max(0, mu - N_eq));
% 
%     % 稳健的二分法设置
%     lower_bound = min(N_eq);  % 最低注水线
%     upper_bound = P + max(N_eq);  % 最高可能注水线 - 移除1000倍系数
%     tol = 1e-10;  % 合理容差
%     max_iter = 1000;
%     mu = (lower_bound + upper_bound) / 2;
% 
%     % 带保护的二分法
%     total_power = 0; % 初始化
%     for iter = 1:max_iter
%         total_power = power_sum(mu);
% 
%         if abs(total_power - P) < tol
%             break;
%         elseif total_power > P
%             upper_bound = mu;
%         else
%             lower_bound = mu;
%         end
% 
%         % 更新并检查边界
%         mu = (lower_bound + upper_bound) / 2;
%         if (upper_bound - lower_bound) < tol
%             break;
%         end
%     end
% 
%     % 最优功率分配 - 移除归一化步骤
%     p_opt = max(0, mu - N_eq);
% 
%     % ===== 关键修正：正确的容量计算 =====
%     % 1. 创建预编码矩阵
%     p_opt_ordered = zeros(size(p_opt)); % 恢复原始顺序
%     p_opt_ordered(idx) = p_opt;
% 
%     % 2. 构建功率分配对角矩阵
%     W = V * diag(sqrt(p_opt_ordered));  % 预编码矩阵
% 
%     % 3. 计算等效信道
%     H_eff = H * W;
% 
%     % 4. 正确的MIMO容量计算
%     C = B * log2(det(eye(size(H,1)) + (1/(N0*B)) * (H_eff * H_eff')));
% 
%     % 5. 验证功率约束
%     total_power_used = sum(p_opt);
%     if abs(total_power_used - P) > 1e-5
%         warning('功率约束未满足: 目标功率 %.6f, 实际功率 %.6f', P, total_power_used);
%     end
% 
%     % ===== 调试输出 =====
%     fprintf('\n===== 注水法结果 =====\n');
%     fprintf('水位线 mu = %.6f\n', mu);
%     fprintf('等效噪声: [');
%     fprintf('%.6f ', N_eq);
%     fprintf(']\n');
%     fprintf('功率分配: [');
%     fprintf('%.6f ', p_opt);
%     fprintf(']\n');
%     fprintf('总功率: %.6f (目标: %.6f)\n', total_power_used, P);
%     fprintf('信道容量: %.6f bps\n', C);
% end

% function [C, p_opt, sigma2_sorted, mu, Ka] = water_filling_capacity_bisect(H, P, N0, B)
%     % 严格保持原始函数签名
%     % 正确实现：等功率分配是最优解
% 
%     % ===== 1. SVD分解 =====
%     [~, S, ~] = svd(H, 'econ');
%     sigma = diag(S);
%     sigma2 = sigma.^2;
% 
%     % ===== 2. 排序信道增益 =====
%     [sigma2_sorted, idx] = sort(sigma2, 'descend');
%     Ka = max(sigma2_sorted) / min(sigma2_sorted);
% 
%     % ===== 3. 等功率分配是最优解 =====
%     N = length(sigma2);
%     p_opt = (P / N) * ones(N, 1);  % 等功率分配
%     mu = NaN;  % 水位线不适用
% 
%     % ===== 4. 容量计算 =====
%     % 标准MIMO容量公式
%     C_Trad = B * log2(det(eye(size(H,1)) + (P/(N0*B*size(H,2))) * (H * H')));
% 
%     % ===== 5. 验证正交位置最优性 =====
%     % 计算等功率分配的容量
%     C_equal = C_Trad;
% 
%     % 计算注水法分配的容量（仅用于比较）
%     noise_power = N0 * B;
%     N_eq = noise_power ./ sigma2_sorted;
% 
%     % 直接计算注水法（不用于实际分配）
%     K = N;
%     mu_temp = (P + sum(N_eq(1:K))) / K;
%     p_opt_temp = max(0, mu_temp - N_eq(1:K));
%     p_opt_temp
%     p_opt_temp = [p_opt_temp; zeros(N-K, 1)];
%     p_opt_temp
% 
%     % 计算注水法容量
%     C_water = 0;
%     for i = 1:N
%         if p_opt_temp(i) > 0
%             snr = (p_opt_temp(i) * sigma2_sorted(i)) / noise_power;
%             C_water = C_water + B * log2(1 + snr);
%         end
%     end
% 
%     % ===== 6. 确保等功率分配最优 =====
%     if C_water > C_equal
%         % 理论上不可能发生，仅用于验证
%         fprintf('警告：注水法容量(%.3e) > 等功率容量(%.3e)\n', C_water, C_equal);
%         fprintf('强制使用等功率分配！\n');
%         p_opt = (P / N) * ones(N, 1);
%         C = C_equal;
%     end
% 
%     % ===== 调试输出 =====
%     fprintf('\n===== 信道容量结果 =====\n');
%     fprintf('等功率分配: [%s]\n', sprintf('%.6f ', p_opt));
%     fprintf('WF-总功率: %.6f (目标: %.6f)\n', sum(p_opt), P);
%     fprintf('信道容量: %.4e bps\n', C);
% 
%     if exist('C_water', 'var')
%         fprintf('注水法计算容量: %.3e bps (仅参考)\n', C_water);
%     end
% end

function [C, p_opt, sigma2_sorted, mu, Ka] = water_filling_capacity_bisect(H, P, N0, B)
    % 注水法计算 MIMO 信道容量（逐步剔除弱信道版本，严谨实现）
    % 输入:
    %   H: MIMO 信道矩阵
    %   P: 总发射功率（单位：W）
    %   N0: 单位带宽噪声功率谱密度（单位：W/Hz）
    %   B: 带宽（单位：Hz）
    %
    % 输出:
    %   C: 最终信道容量（单位：bps）
    %   p_opt: 每个子信道上的最优功率分配（单位：W）
    %   sigma2_sorted: 降序排列的信道增益平方
    %   mu: 水位线
    %   Ka: 信道增益最大比最小的比值（condition metric）

    % 验证H矩阵的准确性
    fprintf('-----------------------------------------------------------\n')
    
    cond_number = cond(H);
    % fprintf('WF: H 的条件数: %.3e\n', cond_number);

    % ===== 1. SVD分解，得到信道增益 =====
    [sigma, sigma2] = compute_singular_values(H);
    

    % ===== 2. 降序排列信道增益平方 =====
    [sigma2_sorted, ~] = sort(sigma2, 'descend');
    sum_sigma2 = sum(sigma2_sorted);
    fprintf('WF: 奇异值平方之和 = %.6e\n', sum_sigma2);

    N = length(sigma2_sorted);
    Ka = max(sigma2_sorted) / min(sigma2_sorted);  % Condition ratio

    % ===== 3. 计算每个子信道的等效噪声功率 =====
    noise_power = N0 * B;
    N_eq = noise_power ./ sigma2_sorted;  % N_i = N0 * B / sigma_i^2

    % ===== 4. 注水法主逻辑：迭代剔除劣质信道，直到所有分配非负 =====
    K = N;
    found = false;
    while ~found
        mu_temp = (P + sum(N_eq(1:K))) / K;
        p_temp = mu_temp - N_eq(1:K);
        if all(p_temp >= 0)
            found = true;
        else
            K = K - 1;
            if K == 0
                error('信道太差，无法进行注水法分配！');
            end
        end
    end

    % ===== 5. 构造完整功率分配向量 =====
    mu = mu_temp;
    p_opt = zeros(N, 1);
    p_opt(1:K) = mu - N_eq(1:K);  % 只有前 K 个子信道分配功率，其余为 0
    % p_opt(1) = P/2;
    % p_opt(2) = P/2;

    % ===== 6. 计算信道容量 =====
    C = 0;
    snr = zeros(2, 1);
    for i = 1:N
        if p_opt(i) > 0
            snr(i) = (p_opt(i)) / noise_power;
            % fprintf('WF: SNR(%d) = %.3e\n', i, snr);
            C = C + B * log2(1 + snr(i) * sigma2_sorted(i));
            % fprintf('WF: snr(%d) * sigma2_sorted(%d) = %.3e\n',  i, i, snr(i) * sigma2_sorted(i));
        end
    end

    if abs(sigma2(1) - sigma2(2)) < 1e-13
        fprintf('WF: H matrix = ');
        disp(H);
        fprintf('Singular values of H (WF): sigma(1) = %.3e, sigma(2) = %.3e\n', sigma(1), sigma(2));
        fprintf('WF: sigma2(1) = %.3e, sigma2(2) = %.3e\n', sigma2(1), sigma2(2));
        fprintf('WF: sigma2_sorted(1) = %.3e, sigma2_sorted(2) = %.3e\n', sigma2_sorted(1), sigma2_sorted(2));
        for i = 1:N
            fprintf('WF: SNR(%d) = %.3e\n', i, snr(i));
            fprintf('WF: snr(%d) * sigma2_sorted(%d) = %.3e\n',  i, i, snr(i) * sigma2_sorted(i));
        end
        fprintf('信道容量 C = %.3e bps\n', C);
        fprintf('\n===== 注水法信道容量计算结果 =====\n');
        fprintf('水位线 mu = %.6e\n', mu);
        fprintf('功率分配（前 %d 个子信道）:\n', K);
        disp(p_opt');
        fprintf('总分配功率 = %.6f W（目标 = %.6f W）\n', sum(p_opt), P);
        fprintf('\n=========================================  END  ===================================================\n');

    end

    % ===== 7. 调试输出（可删） =====
    
    
    % fprintf('\n=========================================  END  ===================================================\n');
end
