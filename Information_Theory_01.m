% LoS MIMO 2x2 容量计算示范代码
% 输入参数 a,b (天线间距米), D (距离米), f (频率Hz), P (发射功率), N0 (噪声功率)
%%%%%%%%%%%% ---Information Theory --- %%%%%%%%%%%%%%

%%
global c
c = 3e8;
%% over
% plot_capacity_vs_spacing_symmetric_geometry(); % 计算对称情况下信道容量随a的变化
function plot_capacity_vs_spacing_symmetric_geometry()
% plot_capacity_vs_spacing 绘制对称天线间距 a=b 情况下的 2x2 LoS MIMO 容量变化曲线

    % 参数设置
    global c
    
    a_vals = linspace(0.5, 10, 100);  % 天线间距从0.5到10米
    b_vals = a_vals;                  % 对称情况，b = a
    D = 10e3;                        % 链路长度（18 km）
    f = 18e9;                         % 载波频率（18 GHz）
    lambda = c/f;
    P = 1;                            % 发射功率（归一化）
    N0 = 4e-21;                        % 噪声功率
    B = 120e6;                           % 信道带宽

    % 初始化结果向量
    capacity_results = zeros(length(a_vals),1);

    % 遍历每个间距，计算容量
    for i = 1:length(a_vals)
        a = a_vals(i);
        b = b_vals(i);
        H = create_H_matrix(a, b, D, f);
        capacity_results(i) = los_mimo_capacity(a, b, D, f, P, N0, B);
    end

    % 绘图
    figure;
    plot(a_vals, capacity_results, '-o', 'LineWidth', 1.5);
    xlabel('Antenna Spacing a = b (m)');
    ylabel('Capacity (bits/s)');
    title('2x2 LoS MIMO Capacity vs Symmetric Antenna Spacing');
    grid on;

end

function H = create_H_matrix(a, b, D, f)
% 计算路径距离 d_ij
    
    global c;
    d11 = D;
    d12 = sqrt(D^2 + b^2);
    d21 = sqrt(D^2 + a^2);
    d22 = sqrt(D^2 + (a - b)^2);
    lambda = c / f; % 波长
    
    % 计算相位 phi_ij
    phi11 = 2*pi*d11/lambda;
    phi12 = 2*pi*d12/lambda;
    phi21 = 2*pi*d21/lambda;
    phi22 = 2*pi*d22/lambda;
    
    % 构造信道矩阵 H
    H = [ (1/d11)*exp(-1j * phi11),   (1/d12)*exp(-1j * phi12);
          (1/d21)*exp(-1j * phi21),   (1/d22)*exp(-1j * phi22) ];
end
%% 计算信道容量 over
function C = los_mimo_capacity(H, P, N0, B)
    global c        
    
    % 计算容量
    Nt = 2; % 发射天线数
    SNR = P/(B*N0); % 信噪比
    
    % 计算容量 C = log2(det(I + (P/N0)*H*H^H))
    C = B * real(log2(det(eye(Nt) + (SNR/Nt)*(H*H'))));
end


%% a和b变化，其他固定，得到容量C-over
% plot_LoS_MIMO_capacity() % 计算不同a和b的情况下的信道容量
function plot_LoS_MIMO_capacity()
    global c
    if isempty(c)
        c = 3e8; % 光速 m/s
    end

    % 固定参数
    f = 18e9;           % 载波频率 18 GHz
    lambda = c / f;     % 波长
    D = 10e3;           % 链路距离 10 km
    P = 1;              % 归一化功率
    N0 = 1e-10;          % 噪声功率
    B = 28e6;             % 信道带宽

    % 设置采样范围（单位米）
    a_range = [0.5, 30]; % 发射端天线间距范围
    b_range = [0.5, 30]; % 接收端天线间距范围
    num_points = 1000;    % 采样点数

    % 生成天线间距向量
    a_vals = linspace(a_range(1), a_range(2), num_points);
    b_vals = linspace(b_range(1), b_range(2), num_points);

    Capacity = zeros(num_points, num_points);

    % 双层循环计算容量
    for i = 1:num_points
        for j = 1:num_points
            a = a_vals(i);
            b = b_vals(j);
            H = create_H_matrix(a, b, D, f);
            Capacity(i,j) = los_mimo_capacity(H, P, N0, B);
        end
    end

    % 绘制3D曲面图
    figure;
    surf(b_vals, a_vals, Capacity);
    
    xlabel('接收端天线间距 b (m)');
    ylabel('发射端天线间距 a (m)');
    zlabel('容量 C (bps/Hz)');
    title(sprintf('2x2 LoS MIMO 容量三维曲面图\n f=%.1f GHz, D=%.1f km', f/1e9, D/1e3));
                
    colorbar;
    shading interp; % 平滑曲面色彩
    view(45,30); % 设置视角，自己拖动也可以
    grid on;
end

%% water-filling 优化不同a,b 下的信道容量
capacity_vs_ab_with_waterfill([0.5, 50], [0.5, 50], 10e3, 18e9, 1, 4e-21, 1000);
function [C_uniform, best_a_uniform, best_b_uniform, C_waterfill, WF_best_a, WF_best_b] = capacity_vs_ab_with_waterfill( ...
    a_range, b_range, D, f, P, N0, num_points)

    global c;    
    B = 120e6;      % 信道带宽Hz
    
    a_vals = linspace(a_range(1), a_range(2), num_points);
    b_vals = linspace(b_range(1), b_range(2), num_points);

    C_uniform = zeros(num_points, num_points);
    C_waterfill = zeros(num_points, num_points);
    sigma_values = zeros(num_points * num_points, 4); % 列为 [a, b, sigma1, sigma2]

    idx = 0; % 用于线性索引奇异值表格行
    
    for i = 1:num_points
        for j = 1:num_points
            a = a_vals(i);
            b = b_vals(j);

            H = create_H_matrix(a, b, D, f);

            % 计算奇异值
            S = svd(H);
            sigma1 = S(1);
            sigma2 = S(2);

            idx = idx + 1;
            sigma_values(idx, :) = [a, b, sigma1, sigma2];

            % 1) 计算均匀功率分配容量
            C_uniform(i,j) = los_mimo_capacity(H, P, N0, B);

            % 2) 计算水填充功率分配容量
            [C_wf, ~, ~] = water_filling_capacity_bisect(H, P, N0, B);
            C_waterfill(i,j) = C_wf;
        end
    end

    % % 保存奇异值数据到csv
    % T = array2table(sigma_values, 'VariableNames', {'a', 'b', 'sigma1', 'sigma2'});
    % writetable(T, 'sigma_values.csv');

    % 找最大容量和对应位置，画图部分同你代码...
    [max_C, linear_idx] = max(C_waterfill(:));
    [row_idx, col_idx] = ind2sub(size(C_waterfill), linear_idx);
    WF_best_a = a_vals(row_idx);
    WF_best_b = b_vals(col_idx);
    fprintf('最大水填充容量为 %.4f bits/s/Hz，发生在 a = %.4f m, b = %.4f m\n', max_C, WF_best_a, WF_best_b);

    [max_C_uniform, idx_uniform] = max(C_uniform(:));
    [r_u, c_u] = ind2sub(size(C_uniform), idx_uniform);
    best_a_uniform = a_vals(r_u);
    best_b_uniform = b_vals(c_u);
    fprintf('最大均匀容量为 %.4f bits/s/Hz，发生在 a = %.4f m, b = %.4f m\n', max_C_uniform, best_a_uniform, best_b_uniform);

    figure;
    subplot(1,2,1);
    imagesc(b_vals, a_vals, C_uniform);
    hold on; % 保留图像坐标轴，允许添加标记
    plot(best_b_uniform, best_a_uniform, 'r*', 'MarkerSize', 10); % 红色星号标记最大点
    text(best_b_uniform + 0.1, best_a_uniform - 0.7, ...
     sprintf('(%.1f, %.1f)', best_b_uniform, best_a_uniform), ...
     'Color', 'black', 'FontSize', 8);

    hold off;
    colorbar;
    xlabel('b (m)');
    ylabel('a (m)');
    title('均匀功率分配容量 (bits/s/Hz)');
    set(gca,'YDir','normal');

    subplot(1,2,2);
    imagesc(b_vals, a_vals, C_waterfill);
    hold on; % 保留图像坐标轴，允许添加标记
    plot(WF_best_b, WF_best_a, 'r*', 'MarkerSize', 10); % 红色星号标记最大点
    text(WF_best_b + 0.1, WF_best_a - 0.6, ...
     sprintf('(%.1f, %.1f)', WF_best_b, WF_best_a), ...
     'Color', 'black', 'FontSize', 8);

    hold off;
    colorbar;
    xlabel('b (m)');
    ylabel('a (m)');
    title('Water-filling功率分配容量 (bits/s/Hz)');
    set(gca, 'YDir', 'normal');

end


%% different freq v.s. capacity
% plot_capacity_heatmaps_all_freqs(); % 计算不同频率下，变化的a和b情况下的信道容量
function plot_capacity_heatmaps_all_freqs()
% plot_capacity_heatmaps_all_freqs
% 绘制多个频率下的2x2 LoS MIMO容量热图

    % 固定参数
    D = 10e3;          % 链路距离 10 km
    P = 1;              % 归一化功率
    N0 = 1e-9;          % 噪声功率
    num_points = 1000;   % a, b 采样点数
    B = 56e6;
    a_vals = linspace(0.5, 10, num_points); % 发射端天线间距
    b_vals = linspace(0.5, 10, num_points); % 接收端天线间距

    % 频率列表
    freq_list = [18e9, 80e9, 170e9]; % Hz
    freq_labels = {'18 GHz', '80 GHz', '170 GHz'};

    % 绘图
    
    figure;
    fig = gcf;
    tiledlayout(1, length(freq_list), 'Padding', 'compact');
    for idx = 1:length(freq_list)
        f = freq_list(idx);
        Capacity = zeros(num_points);

        for i = 1:num_points
            for j = 1:num_points
                a = a_vals(i);
                b = b_vals(j);
                H = create_H_matrix(a, b, D, f);
                Capacity(i,j) = los_mimo_capacity(H, P, N0, B);
            end
        end

        nexttile;
        imagesc(b_vals, a_vals, Capacity);
        xlabel('接收端天线间距 b (m)');
        ylabel('发射端天线间距 a (m)');
        title(['频率 = ' freq_labels{idx}]);
        colorbar;
        set(gca, 'YDir', 'normal');
        axis square;  % 保证正方形
    end

    sgtitle(sprintf('2x2 LoS MIMO 容量热图（D = %.1f km）', D / 1e3));
end

%% 不同距离D的情况下的信道容量分析
% plot_capacity_heatmaps_freq_distance(80e9, [1e3, 5e3]);
function plot_capacity_heatmaps_freq_distance(freq, D_range)
% 绘制固定频率下，不同链路距离 D ∈ [Dmin, Dmax] 下的容量热图
% 并自动标记最大容量点
% 参数:
%   freq    - 频率 (Hz)，如 80e9
%   D_range - 距离范围 [Dmin, Dmax]，单位 m，如 [1000, 20000]

    % 固定参数
    global c
    P = 1; 
    lambda = c / freq;     % 波长
    N0 = 1e-9;          
    num_points = 500; 
    B = 56e6;
    a_vals = linspace(0.5, 10, num_points); 
    b_vals = linspace(0.5, 10, num_points); 
    
    % 自动生成距离列表（默认取4个等间距点）
    Dmin = D_range(1);
    Dmax = D_range(2);
    num_D = 6;
    D_list = linspace(Dmin, Dmax, num_D);

    % 智能布局：自动行列
    n_cols = ceil(sqrt(num_D));
    n_rows = ceil(num_D / n_cols);

    % 图设置
    figure;
    fig = gcf;
    fig.Position = [100, 100, 300 * n_cols, 300 * n_rows];
    tiledlayout(n_rows, n_cols, 'Padding', 'compact', 'TileSpacing', 'compact');

    for idx = 1:num_D
        D = D_list(idx);
        Capacity = zeros(num_points, num_points);

        for i = 1:num_points
            for j = 1:num_points
                a = a_vals(i);
                b = b_vals(j);
                H = create_H_matrix(a, b, D, freq);
                Capacity(i,j) = los_mimo_capacity(H, P, N0, B);
            end
        end

        % 找最大容量点
        [C_max, idx_max] = max(Capacity(:));

        [i_opt, j_opt] = ind2sub(size(Capacity), idx_max);
        a_opt = a_vals(i_opt);
        b_opt = b_vals(j_opt);

        % 绘图
        nexttile;
        imagesc(b_vals, a_vals, Capacity);
        hold on;
        scatter(b_opt, a_opt, 50, 'ro', 'filled'); % 红色圆点
        % 标注位置和容量值，简洁写坐标和容量
        text(b_opt, a_opt, ...
            sprintf('%.2f, %.2f\nC=%.3f', b_opt, a_opt, C_max), ...
            'Color', 'red', 'FontWeight', 'bold', 'FontSize', 9, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

        ab_target = lambda * D /2;
        a_plot = linspace(min(a_vals), max(a_vals), 1000);
        b_plot = ab_target ./ a_plot;
        plot(b_plot, a_plot, 'black--', 'LineWidth', 1.5);

        hold off;

        xlabel('接收端天线间距 b (m)');
        ylabel('发射端天线间距 a (m)');
        title(sprintf('D = %.1f km', D / 1000));
        colorbar;
        set(gca, 'YDir', 'normal');
        axis square;
    end

    sgtitle(sprintf('频率 = %.1f GHz 不同链路距离下的容量热图（红点为容量最大点）', freq / 1e9));
end



%% Water-Filling 优化信道功率分配
H = create_H_matrix(50, 50, 15e3, 18e9);  % 更大奇异值差异
P = 1;
N0 = 1e-14;
B = 56e6;
[C, p, gain] = water_filling_capacity_bisect(H, P, N0, B);
disp('信道增益:'); disp(gain);
disp('功率分配:'); disp(p);
disp(['最大容量: ', num2str(C)]);

function [C, p_opt, sigma2] = water_filling_capacity_bisect(H, P, N0, B)
    [~, S, ~] = svd(H);
    sigma = diag(S);
    sigma2 = sigma.^2;
    n = length(sigma2);
    
    % 定义求和函数
    power_sum = @(mu) sum(max(0, mu - N0 ./ sigma2));
    
    % 二分法参数
    lower = min(N0 ./ sigma2);
    upper = P + max(N0 ./ sigma2);
    tol = 1e-9;
    max_iter = 1000;
    
    % 二分法求mu
    for iter = 1:max_iter
        mu = (lower + upper) / 2;
        val = power_sum(mu) - P;
        if abs(val) < tol
            break;
        elseif val > 0
            upper = mu;
        else
            lower = mu;
        end
    end
    
    % 计算功率分配
    p_opt = max(0, mu - N0 ./ sigma2);
    
    
    % 计算容量
    C = B * sum(log2(1 + p_opt .* sigma2 / (B *N0)));
end

%% 信道容量C和风吹位移的关系
[C_uniform, best_a_Uni, best_b_Uni, C_waterfill, best_a_WF, best_b_WF] = capacity_vs_ab_with_waterfill( ...
    [0.5, 10], [0.5, 10], 10e3, 18e9, 1, 4e-21, 1000);
capacity_vs_displacement(best_a_WF, best_b_WF, best_a_Uni, best_b_Uni);

function capacity_vs_displacement(best_a_WF, best_b_WF, best_a_Uni, best_b_Uni)
    % 参数设定
    D = 1.5e4;        % 传输距离 15km
    f = 18e9;         % 载波频率 18GHz
    P = 1;            % 总功率 (W)
    N0 = 4e-21;       % 噪声功率谱密度 (W/Hz)
    B = 120e6;        % 带宽 120MHz
    x_vals = linspace(-1, 1, 100); % 风吹导致的偏移，-1~1米

    C_uniform = zeros(size(x_vals));
    C_waterfill = zeros(size(x_vals));

    % 计算偏移对容量的影响
    for idx = 1:length(x_vals)
        x = x_vals(idx);
        
        % 实际天线间距加上位移 (Uni)
        a_Uni = best_a_Uni + x; 
        b_Uni = best_b_Uni;
        
        % 实际天线间距加上位移 (WF)
        a_WF = best_a_WF + x; 
        b_WF = best_b_WF;
        
        % 计算信道矩阵H
        H_Uni = create_H_matrix(a_Uni, b_Uni, D, f);
        H_WF = create_H_matrix(a_WF, b_WF, D, f);
        
        % 计算容量
        C_uniform(idx) = los_mimo_capacity(H_Uni, P, N0, B); % 传统均分功率
        [C_waterfill(idx), ~, ~] = water_filling_capacity_bisect(H_WF, P, N0, B); % Water-filling
    end

    % 画图对比容量变化
    figure;
    plot(x_vals, C_uniform, '-b', 'LineWidth', 2);
    hold on;
    plot(x_vals, C_waterfill, '-r', 'LineWidth', 2);
    grid on;
    xlabel('天线间距偏移 x (m)');
    ylabel('信道容量 (bits/s)');
    legend('均匀功率分配', 'Water-filling功率分配', 'Location', 'Best');
    title(sprintf('天线间距偏移对容量影响 (a_{WF} = %.1fm, b_{WF} = %.1fm)', best_a_WF, best_b_WF));
    hold off;
    
    % ==== 计算并打印最佳点处的奇异值 ====
    H_best = create_H_matrix(best_a_WF, best_b_WF, D, f);
    s_best = svd(H_best);
    if length(s_best) >= 2
        fprintf('最大容量点对应的奇异值：sigma1 = %.4e, sigma2 = %.4e\n', ...
                s_best(1), s_best(2));
    else
        fprintf('警告：信道矩阵奇异值不足2个\n');
    end
    
    % ==== 构建奇异值热图矩阵 ====
    num_points = 100;
    a_range = linspace(best_a_WF - 1, best_a_WF + 1, num_points);
    b_range = linspace(best_b_WF - 1, best_b_WF + 1, num_points);
    
    sigma1_map = zeros(num_points, num_points);
    sigma2_map = zeros(num_points, num_points);
    
    % 计算热图区域内的奇异值
    for i = 1:num_points
        for j = 1:num_points
            H = create_H_matrix(a_range(i), b_range(j), D, f);
            s = svd(H);
            sigma1_map(i, j) = s(1);
            if length(s) >= 2
                sigma2_map(i, j) = s(2);
            else
                sigma2_map(i, j) = 0; % 处理奇异值不足的情况
            end
        end
    end

    % 绘制奇异值热图
    figure;
    subplot(1, 2, 1);
    imagesc(b_range, a_range, sigma1_map);
    colorbar;
    xlabel('天线间距 b (m)');
    ylabel('天线间距 a (m)');
    title('奇异值 \sigma_1');
    axis xy; % 确保坐标轴方向正确
    hold on;
    plot(best_b_WF, best_a_WF, 'rx', 'MarkerSize', 10, 'LineWidth', 2); % 标记最佳点
    hold off;
    
    subplot(1, 2, 2);
    imagesc(b_range, a_range, sigma2_map);
    colorbar;
    xlabel('天线间距 b (m)');
    ylabel('天线间距 a (m)');
    title('奇异值 \sigma_2');
    axis xy; % 确保坐标轴方向正确
    hold on;
    plot(best_b_WF, best_a_WF, 'rx', 'MarkerSize', 10, 'LineWidth', 2); % 标记最佳点
    hold off;
    
    % 设置主图标题
    sgtitle(sprintf('奇异值分布 (中心点 a=%.1fm, b=%.1fm)', best_a_WF, best_b_WF));
end