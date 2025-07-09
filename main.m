%%%%%%%%% Information Theory Project -- Nokia: LoS of MIMO %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% main %%%%%%%%%%%%%%%%%%%%%%%%%
function main()
    addpath('utils');
    % 参数设置
    c = 3e8;                                 % 光速
    num_points = 2000;
    a_vals = linspace(0.5, 10, num_points);  % 天线间距从0.5到10m
    b_vals = a_vals;                         % 对称情况，b = a
    D_list = [10e3, 5e3, 2e3, 1e3, 500];            % 链路长度 m
    % D_list = [10e3, 5e3, 2000, 1e3, 500];            % 链路长度 m
    Oscillation = 1;                               % Oscillation (m)
    freq_list = [18e9, 38e9, 80e9, 115e9, 170e9];   % Hz
    f = freq_list(1);                  % 载波频率（18 GHz）
    P = 1;                             % 发射功率Watt（归一化）
    N0 = 4e-21;                        % 噪声功率W/Hz (常温下是-174dBm)
    B = 250e6;                         % 信道带宽Hz

    %%
  
    % K取值为奇数，比如前5个奇数
    K_values = [1, 3, 5, 7, 9, 11, 13, 15, 17];

    a_values_test = sqrt((c/freq_list(3) * D_list(3) / 2) .* K_values);
    
    fprintf('正交位置 a=b 时的天线间距（米）:\n');
    for i = 1:length(a_values_test)
        fprintf('K = %d, a = b = %.4f m\n', K_values(i), a_values_test(i));
    end

    %%%%
    % % 在正交位置测试
    % d_orth = sqrt(c/freq_list(3) * D_list(3) / 2);
    % H_orth = create_H_matrix(d_orth, d_orth, D_list(3), freq_list(3), 0);
    % 
    % % 1. 检查H*H'是否对角
    % HH = H_orth * H_orth';
    % disp('H*H^H = ');
    % disp(HH);
    % % 2. 检查奇异值
    % [U,S_test,V] = svd(H_orth);
    % disp('奇异值 = ');
    % disp(diag(S_test));
    % 
    % % 3. 计算理论特征值
    % % 对于正交信道，特征值应相等
    % expected_eig = mean(diag(HH));
    % disp('理论特征值 = ');
    % disp([expected_eig, expected_eig]);
    % C_wf = water_filling_capacity_bisect(H_orth, energy, P, N0, B);
    % C_trad = los_mimo_capacity(H_orth, energy, P, N0, B);
    % fprintf('注水法容量: %.4f\n', C_wf);
    % fprintf('传统法容量: %.4f\n', C_trad);%%%%
    %%
    
    
    % plot_capacity_vs_spacing_symmetric_geometry(a_vals, num_points, D_list, freq_list, N0, B, P); % 计算对称情况下信道容量随a的变化
    % plot_capacity_vs_symmetric_water_filling(a_vals, num_points, D_list, freq_list, N0, B, P);
    plot_capacity_vs_symmetric_water_filling_and_Traditional(a_vals, num_points, D_list, freq_list, N0, B, P);
    % % 为每个频率分别画水填功率分配图
    % for idx = 1:length(freq_list)        
    %     plot_power_distribution_of_WF_per_freq(a_vals, b_vals, num_points, D_list(idx), freq_list(idx), N0, B, P);
    % end
    % 
    % plot_Ka_vs_a(a_vals, num_points, D_list, freq_list, N0, B, P); % 画出不同频率下的条件数折线图2D
    % plot_LoS_MIMO_capacity(a_vals, b_vals, num_points, D_list, freq_list, N0, B, P) % 计算不同a和b的情况下的信道容量
    % plot_LoS_MIMO_capacity_WF(a_vals, b_vals, num_points, D_list, freq_list, N0, B, P) % 计算不同a和b的情况下的信道容量- WF
    % 
    % plot_los_capacity_Oscillations_heatmap(a_vals, num_points, D_list, N0, B, P, freq_list, Oscillation);
    % plot_los_capacity_Oscillations_heatmap_WF(a_vals, num_points, D_list, N0, B, P, freq_list, Oscillation);
    % plot_capacity_heatmap_freq_vs_spacing(a_vals, num_points, D_list, freq_list, N0, B, P)
    % 
    %%%%%%%%%%%%%%%% 下面的好像多余了 %%%%%%%%%%%%%%
    % plot_capacity_heatmaps_all_freqs(a_vals, b_vals, num_points, D_list, N0, B, P, freq_list);
    % plot_capacity_heatmaps_all_freqs_WF(a_vals, b_vals, num_points, D, N0, B, P, freq_list);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --Functions-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 计算对称情况下信道容量随a(b)的变化 - 传统均分功率分配方法
function plot_capacity_vs_spacing_symmetric_geometry(a_vals, num_points, D_list, freq_list, N0, B, P)
% plot_capacity_vs_spacing 绘制对称天线间距 a=b 情况下的 2x2 LoS MIMO 容量变化曲线
    fig = figure; % 新建图窗
    hold on; % 保持画布，绘制多条曲线

    for idx = 1:length(freq_list)
        f = freq_list(idx);
        capacity_results = zeros(num_points,1);
        D = D_list(idx);
        for i = 1:num_points
            a = a_vals(i);
            b = a;
            H = create_H_matrix(a, b, D, f, 0);
            capacity_results(i) = los_mimo_capacity(H, P, N0, B);
        end
        
        plot(a_vals, capacity_results, 'LineWidth', 1.5, 'DisplayName', sprintf('Freq = %.0f GHz', f/1e9));
    end

    hold off;
    xlabel('Antenna Spacing a = b (m)');
    ylabel('Capacity (bits/s)');
    title('2x2 LoS MIMO Capacity vs Symmetric Antenna Spacing for Multiple Frequencies');
    grid on;
    legend('show'); % 显示图例
    save_figure_custom(fig, 'Capacity_a_eq_b', '~', B, '~');
end

%% 计算对称情况下信道容量随a(b)的变化 - Water Filling功率分配方法
function plot_capacity_vs_symmetric_water_filling(a_vals, num_points, D_list, freq_list, N0, B, P)
    fig = figure; % 新建图窗
    hold on; % 保持画布，绘制多条曲线
    for idx = 1:length(freq_list)
        f = freq_list(idx);
        capacity_results = zeros(num_points,1);
        D = D_list(idx);
        for i = 1:num_points
            a = a_vals(i);
            b = a;
            H = create_H_matrix(a, b, D, f, 0);
            [capacity_results(i), ~, ~, ~, ~] = water_filling_capacity_bisect(H, P, N0, B);
        end
        plot(a_vals, capacity_results, 'LineWidth', 1.5, 'DisplayName', sprintf('Freq = %.0f GHz', f/1e9));
    end

    hold off;
    xlabel('Antenna Spacing a = b (m)');
    ylabel('Capacity (bits/s)');
    title('2x2 LoS MIMO Capacity vs Symmetric Antenna Spacing for Multiple Frequencies-Water Filling');
    grid on;
    legend('show'); % 显示图例
    save_figure_custom(fig, 'Capacity_a_eq_b_WF', '~', B, '~');

end

%% 把WF和传统方法合一起画一个图，单频率的信道容(a=b)
function plot_capacity_vs_symmetric_water_filling_and_Traditional(a_vals, num_points, D_list, freq_list, N0, B, P)
    fig = figure; % 新建图窗
    hold on; % 保持画布，绘制多条曲线

    f = freq_list(3);
    capacity_results_WF = zeros(num_points,1);
    capacity_results_Trad = zeros(num_points,1);
    D = D_list(3);
    for i = 1:num_points
        a = a_vals(i);
        b = a;
        H = create_H_matrix(a, b, D, f, 0);
        [capacity_results_WF(i), ~, ~, ~, ~] = water_filling_capacity_bisect(H, P, N0, B);
        capacity_results_Trad(i)= los_mimo_capacity(H, P, N0, B);
    end
    plot(a_vals, capacity_results_WF, 'LineWidth', 1.5, 'DisplayName', sprintf('WF: Freq = %.0f GHz', f/1e9));
    plot(a_vals, capacity_results_Trad, 'LineWidth', 1.5, 'LineStyle', '--','DisplayName', sprintf('Trad: Freq = %.0f GHz', f/1e9));

    hold off;
    xlabel('Antenna Spacing a = b (m)');
    ylabel('Capacity (bits/s)');
    title('2x2 LoS MIMO Capacity vs Symmetric Antenna Spacing for Multiple Frequencies-Water Filling and Tradional');
    grid on;
    legend('show'); % 显示图例
    save_figure_custom(fig, 'Capacity_a_eq_b_WFandTrad', '~', B, '~');

end

%% 条件数 (Ka) 展示
function plot_Ka_vs_a(a_vals, num_points, D_list, freq_list, N0, B, P)
    % ====== 开始画图（Ka vs antenna spacing）======

    for idx = 1:length(freq_list)
        fig_Ka = figure;   
        f = freq_list(idx);
        D = D_list(idx);
        for i = 1:num_points
            a = a_vals(i);
            b = a;
            H = create_H_matrix(a, b, D, f, 0);
            [~, ~, ~, ~, Ka(i, idx)] = water_filling_capacity_bisect(H, P, N0, B);
        end

        semilogy(a_vals, log2(Ka(:, idx)), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Freq = %.1f GHz', freq_list(idx)/1e9));
        % === 图形标签 ===
        xlabel('Antenna Spacing a = b (m)', 'FontSize', 12);    % 横坐标
        ylabel('Condition Number K_a (log2 scale)', 'FontSize', 12); % 纵坐标
        title('Condition Number vs. Symmetric Antenna Spacing for LoS MIMO', 'FontSize', 13); % 标题
        legend('show', 'Location', 'best');                      % 显示图例
        grid on;
        % ylim([0, 2]);  % 因为Ka≥1，所以下限设1
        save_figure_custom(fig_Ka, 'Num_Condition', f, B, D);
        
    end

end
%% 单频率单图绘制功率分配 (only WF has this)
function plot_power_distribution_of_WF_per_freq(a_vals, b_vals, num_points, D, f, N0, B, P)
    fig = figure;
    hold on;
    colors = lines(2);  % 2条线，2个功率分配分量
    
    % 准备存储
    p_opt_all = zeros(2, num_points);  % 2个信道功率分配
    sigma2_all = zeros(2, num_points); % 对应奇异值平方
    
    for i = 1:num_points
        a = a_vals(i);
        b = b_vals(i);
        H = create_H_matrix(a, b, D, f, 0);
        [~, p_opt, sigma2, ~] = water_filling_capacity_bisect(H, P, N0, B);
        p_opt_all(:,i) = p_opt(:);
        sigma2_all(:,i) = sigma2(:);
    end

    % 画功率分配折线图
    for k = 1:2
        plot(a_vals, p_opt_all(k,:)/P, '-o', 'Color', colors(k,:), 'DisplayName', sprintf('Stream %d', k));
    end

    hold off;
    xlabel('Antenna Spacing a = b (m)');
    ylabel('Power Allocation (Watt)');
    title(sprintf('Water Filling Power Distribution @ %.0f GHz', f/1e9));
    legend('show');
    grid on;
    ylim([0, 1]); % 
    save_figure_custom(fig, 'Distribution', f, B, D);
end


%% a和b变化，其他固定，得到容量C 热图, 传统均分功率方法
function plot_LoS_MIMO_capacity(a_vals, b_vals, num_points, D_list, freq_list, N0, B, P)
    freq_list = [18e9, 80e9, 170e9];
    D_list = [10e3, 2e3, 500];
    num_freqs = length(freq_list);
    Capacity_all = cell(1, num_freqs);
    freq_labels = {'18 GHz', '80 GHz', '170 GHz'};

    % 计算容量矩阵
    for idx = 1:num_freqs
        D = D_list(idx);
        f = freq_list(idx);
        Capacity = zeros(num_points, num_points);

        for i = 1:num_points
            for j = 1:num_points
                a = a_vals(i);
                b = b_vals(j);
                H = create_H_matrix(a, b, D, f, 0);
                Capacity(i,j) = los_mimo_capacity(H, P, N0, B);
            end
        end
        Capacity_all{idx} = Capacity;
    end

    % 计算行列数，接近正方形布局
    num_rows = ceil(sqrt(num_freqs));
    num_cols = ceil(num_freqs / num_rows);

    % 创建布局
    fig = figure('Name', '2x2 LoS MIMO 容量热图 多频率', 'NumberTitle', 'off');
    tlo = tiledlayout(num_rows, num_cols, 'Padding', 'compact', 'TileSpacing', 'compact');

    is3D = false;
    draw(a_vals, b_vals, freq_list, freq_labels, Capacity_all, is3D);

    % 切换按钮
    uicontrol('Style', 'pushbutton', ...
              'String', '2D/3D', ...
              'Position', [20 20 80 30], ...
              'Callback', @(src, event) toggleView());

    function toggleView()
        is3D = ~is3D;
        draw(a_vals, b_vals, freq_list, freq_labels, Capacity_all, is3D);
    end
    save_figure_custom(fig, 'Capacity_vs_a_b', '~', '~', '~');
end

%% a和b变化，其他固定，得到容量C 热图, Water Filling 方法
function plot_LoS_MIMO_capacity_WF(a_vals, b_vals, num_points, D_list, freq_list, N0, B, P)
    freq_list = [18e9, 80e9, 170e9];
    D_list = [10e3, 2e3, 500];
    num_freqs = length(freq_list);
    Capacity_all = cell(1, num_freqs);
    freq_labels = {'18 GHz', '80 GHz', '170 GHz'};

    % 计算容量矩阵
    for idx = 1:num_freqs
        D = D_list(idx);
        f = freq_list(idx);
        Capacity = zeros(num_points, num_points);

        for i = 1:num_points
            for j = 1:num_points
                a = a_vals(i);
                b = b_vals(j);
                H = create_H_matrix(a, b, D, f, 0);
                Capacity(i,j) = water_filling_capacity_bisect(H, P, N0, B);
            end
        end
        Capacity_all{idx} = Capacity;
    end

    % 计算行列数，接近正方形布局
    num_rows = ceil(sqrt(num_freqs));
    num_cols = ceil(num_freqs / num_rows);

    % 创建布局
    fig = figure('Name', '2x2 LoS MIMO Capacity Heatmap Multi-frequency - WF', 'NumberTitle', 'off');
    tlo = tiledlayout(num_rows, num_cols, 'Padding', 'compact', 'TileSpacing', 'compact');

    is3D = false;
    draw(a_vals, b_vals, freq_list, freq_labels, Capacity_all, is3D);

    % 切换按钮
    uicontrol('Style', 'pushbutton', ...
              'String', '2D/3D', ...
              'Position', [20 20 80 30], ...
              'Callback', @(src, event) toggleView());

    function toggleView()
        is3D = ~is3D;
        draw(a_vals, b_vals, freq_list, freq_labels, Capacity_all, is3D);
    end
    save_figure_custom(fig, 'WF_Capacity_vs_a_b', '~', '~', '~');
end



%% 天线间距波动x = (-1m, +1m) 对信道容量的影响, 此时a=b, 增加了参数 x=Oscillation
function plot_los_capacity_Oscillations_heatmap(a_vals, num_points, D_list, N0, B, P, freq_list, Oscillation)
    % Oscillation 范围：正负 x 米的扰动
    freq_list = [18e9, 80e9, 170e9];
    D_list = [10e3, 2e3, 500];
    num_freqs = length(freq_list);
    freq_labels = {'18 GHz', '80 GHz', '170 GHz'};
    x_vals = linspace(-Oscillation, Oscillation, num_points); % 垂直方向
    Capacity_all = cell(1, num_freqs);

    % 计算容量矩阵
    for idx = 1:num_freqs
        D = D_list(idx);
        f = freq_list(idx);
        C_map = zeros(num_points); % x (行) × a (列)

        for i = 1:num_points  % Oscillation
            x_shift = x_vals(i);
            for j = 1:num_points  % a
                a = a_vals(j);
                b = a;  % 对称
                H = create_H_matrix(a, b, D, f, x_shift);
                C_map(i, j) = los_mimo_capacity(H, P, N0, B);
            end
        end
        Capacity_all{idx} = C_map;
    end

    % 图像窗口和布局
    fig = figure('Name', 'LoS MIMO Capacity Heatmap: Antenna Shift vs. a', 'NumberTitle', 'off');
    num_rows = ceil(sqrt(num_freqs));
    num_cols = ceil(num_freqs / num_rows);
    tlo = tiledlayout(num_rows, num_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    is3D = false;  % 初始为2D视图
    draw();        % 画初始图

    % 切换按钮
    uicontrol('Style', 'pushbutton', ...
              'String', '2D/3D', ...
              'Position', [20 20 80 30], ...
              'Callback', @(src, event) toggle_view());

    sgtitle(sprintf('Channel Capacity Heatmap: Impact of a and Antenna Displacement Max Amplitude x = %.1f m - Conventional Method', Oscillation));

    % 回调函数：切换视图
    function toggle_view()
        is3D = ~is3D;
        draw();
    end

    % 绘图函数
    function draw()
        for idx = 1:num_freqs
            nexttile(idx);
            cla;
            [A_mesh, X_mesh] = meshgrid(a_vals, x_vals);
            C = Capacity_all{idx};

            if is3D
                surf(A_mesh, X_mesh, C, 'EdgeColor', 'none');
                view(45, 30);
            else
                imagesc(a_vals, x_vals, C);
                set(gca, 'YDir', 'normal');
            end

            xlabel('a = b (m)');
            ylabel('RX Antenna-1 Shift x (m)');
            title(freq_labels{idx});
            colorbar;
        end
    end
    save_figure_custom(fig, 'Capacity_vs_Oscillation', '~', '~', '~');
end


%% WF Oscillation
function plot_los_capacity_Oscillations_heatmap_WF(a_vals, num_points, D_list, N0, B, P, freq_list, Oscillation)
    % Oscillation 范围：正负 x 米的扰动
    freq_list = [18e9, 80e9, 170e9];
    D_list = [10e3, 2e3, 500];
    num_freqs = length(freq_list);
    freq_labels = {'18 GHz', '80 GHz', '170 GHz'};
    x_vals = linspace(-Oscillation, Oscillation, num_points); % 垂直方向
    Capacity_all = cell(1, num_freqs);

    % 计算容量矩阵
    for idx = 1:num_freqs
        D = D_list(idx);
        f = freq_list(idx);
        C_map = zeros(num_points); % x (行) × a (列)

        for i = 1:num_points  % Oscillation
            x_shift = x_vals(i);
            for j = 1:num_points  % a
                a = a_vals(j);
                b = a;  % 对称
                H = create_H_matrix(a, b, D, f, x_shift);
                C_map(i, j) = water_filling_capacity_bisect(H, P, N0, B);
            end
        end
        Capacity_all{idx} = C_map;
    end

    % 图像窗口和布局
    fig = figure('Name', 'LoS MIMO Capacity Heatmap - Antenna Displacement vs a - WF', 'NumberTitle', 'off');
    num_rows = ceil(sqrt(num_freqs));
    num_cols = ceil(num_freqs / num_rows);
    tlo = tiledlayout(num_rows, num_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    is3D = false;  % 初始为2D视图
    draw();        % 画初始图

    % 切换按钮
    uicontrol('Style', 'pushbutton', ...
              'String', '2D/3D', ...
              'Position', [20 20 80 30], ...
              'Callback', @(src, event) toggle_view());

    sgtitle(sprintf('Capacity (heatmap)：The impact of "a" and Oscillation x = %.1f m - WF', Oscillation));

    % 回调函数：切换视图
    function toggle_view()
        is3D = ~is3D;
        draw();
    end

    % 绘图函数
    function draw()
        for idx = 1:num_freqs
            nexttile(idx);
            cla;
            [A_mesh, X_mesh] = meshgrid(a_vals, x_vals);
            C = Capacity_all{idx};

            if is3D
                surf(A_mesh, X_mesh, C, 'EdgeColor', 'none');
                view(45, 30);
            else
                imagesc(a_vals, x_vals, C);
                set(gca, 'YDir', 'normal');
            end

            xlabel('a = b (m)');
            ylabel('RX1 moves x (m)');
            title(freq_labels{idx});
            colorbar;
        end
    end
    save_figure_custom(fig, 'WF_Capacity_vs_Oscillation', '~', '~', '~');
end

%% Capacity vs. freq
function plot_capacity_heatmap_freq_vs_spacing(a_vals, num_points, D_list, freq_list, N0, B, P)
% 绘制 LoS 2x2 MIMO 系统中频率 vs. 天线间距 下的信道容量热图
% a_vals: 一组天线间距值 (a = b)
% freq_list: 一组载波频率值
% 其他参数与之前相同

    num_f = length(freq_list);
    capacity_map = zeros(num_a, num_f);

    for i = 1:num_points
        a = a_vals(i);
        b = a; % 对称结构
        for j = 1:num_f
            f = freq_list(j);
            H = create_H_matrix(a, b, D, f, 0);
            [C, ~, ~, ~, ~] = water_filling_capacity_bisect(H, P, N0, B);
            capacity_map(i, j) = C;
        end
    end

    % 绘制热图
    fig = figure;
    imagesc(freq_list/1e9, a_vals, capacity_map); % 横轴频率 (GHz)，纵轴天线间距 (m)
    colorbar;
    xlabel('Carrier Frequency (GHz)', 'FontSize', 12);
    ylabel('Antenna Spacing a = b (m)', 'FontSize', 12);
    title('Capacity Heatmap vs. Frequency and Antenna Spacing', 'FontSize', 13);
    set(gca, 'YDir', 'normal'); % 保持 a 从小到大向上画
    save_figure_custom(fig, 'WF_Capacity_vs_freq_3D', '~', '~', '~');
end



%%%%%%%%%%%%%%%%%%%%%%%%%  下面的好像是重复了，白写了操 %%%%%%%%%%%%%%%%%%%%%
%% 传统方法的功率分配：不同频率下的信道容量展示
function plot_capacity_heatmaps_all_freqs(a_vals, b_vals, num_points, D_list, N0, B, P, freq_list)
    freq_labels = {'18 GHz', '38 GHz', '80 GHz', '115 GHz', '170 GHz'};
    fig = figure('Name', '2x2 LoS MIMO Capacity Heatmaps', 'NumberTitle', 'off');
    tlo = tiledlayout(1, length(freq_list), 'Padding', 'compact');

    % 状态变量：是否为3D视图
    is3D = false;

    % 数据缓存
    Capacity_all = cell(1, length(freq_list));

    % 先预计算所有频率下的容量
    for idx = 1:length(freq_list)
        f = freq_list(idx);
        D = D_list(idx);
        Capacity = zeros(num_points);

        for i = 1:num_points
            for j = 1:num_points
                a = a_vals(i);
                b = b_vals(j);
                H = create_H_matrix(a, b, D, f, 0);
                Capacity(i,j) = los_mimo_capacity(H, P, N0, B);
            end
        end

        Capacity_all{idx} = Capacity;
    end   

    % 计算行列数，接近正方形布局
    num_rows = ceil(sqrt(length(freq_list)));
    num_cols = ceil(length(freq_list) / num_rows);

    % 创建布局
    tlo = tiledlayout(num_rows, num_cols, 'Padding', 'compact', 'TileSpacing', 'compact');

    draw(a_vals, b_vals, freq_list, freq_labels, Capacity_all, is3D);

    % 添加切换按钮
    uicontrol('Style', 'pushbutton', ...
              'String', '切换 2D/3D', ...
              'Position', [10 10 100 30], ...
              'Callback', @(src, event) toggle_view());

    % 回调函数：切换状态并重绘
    function toggle_view()
        is3D = ~is3D;
        draw(a_vals, b_vals, freq_list, freq_labels, Capacity_all, is3D);
    end

    sgtitle(sprintf('2x2 LoS MIMO 容量图 - 传统均分功率方法'));
end


%% Water Filling Algorithm 功率分配：不同频率下的信道容量展示
function plot_capacity_heatmaps_all_freqs_WF(a_vals, b_vals, num_points, D_list, N0, B, P, freq_list)
    freq_labels = {'18 GHz', '80 GHz', '170 GHz'};
    fig = figure('Name', '2x2 LoS MIMO Capacity Heatmaps', 'NumberTitle', 'off');
    tlo = tiledlayout(1, length(freq_list), 'Padding', 'compact');

    % 状态变量：是否为3D视图
    is3D = false;

    % 数据缓存
    Capacity_all = cell(1, length(freq_list));

    % 先预计算所有频率下的容量
    for idx = 1:length(freq_list)
        f = freq_list(idx);
        D = D_list(idx);
        Capacity = zeros(num_points);

        for i = 1:num_points
            for j = 1:num_points
                a = a_vals(i);
                b = b_vals(j);
                H = create_H_matrix(a, b, D, f, 0);
                Capacity(i,j) = water_filling_capacity_bisect(H, P, N0, B);
            end
        end

        Capacity_all{idx} = Capacity;
    end    
    
    draw(a_vals, b_vals, freq_list, freq_labels, Capacity_all, is3D);

    % 添加切换按钮
    uicontrol('Style', 'pushbutton', ...
              'String', '切换 2D/3D', ...
              'Position', [10 10 100 30], ...
              'Callback', @(src, event) toggle_view());

    % 回调函数：切换状态并重绘
    function toggle_view()
        is3D = ~is3D;
        draw(a_vals, b_vals, freq_list, freq_labels, Capacity_all, is3D);
    end

    sgtitle(sprintf('2x2 LoS MIMO 容量图（D = %.1f km）- Water Filling', D / 1e3));
end