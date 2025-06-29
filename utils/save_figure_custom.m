%% %% 稳健的保存图像函数
%%  使用示例(灵活省略参数):
% save_figure_custom(fig, 'Capacity', 80e9, 120e6, 2000);     % 全部有
% save_figure_custom(fig, 'Capacity', 80e9, '~', 2000);       % 跳过带宽
% save_figure_custom(fig, 'Capacity', '~', '~', 5000);        % 只有D
% save_figure_custom(fig, 'Capacity', '~', '~', '~');         % 只有前缀
% 保存结果文件名示例：
% Capacity_Freq_80GHz_BW_120MHz_D_2km.png
% 
% Capacity_Freq_80GHz_D_2km.png
% 
% Capacity_D_5km.png
% 
% Capacity.png
%%
function save_figure_custom(fig, prefix, freq_Hz, B, D)
    name_parts = {prefix};

    if nargin >= 3 && ~isempty(freq_Hz) && ~isequal(freq_Hz, '~')
        name_parts{end+1} = sprintf('Freq_%.0fGHz', freq_Hz/1e9);
    end
    if nargin >= 4 && ~isempty(B) && ~isequal(B, '~')
        name_parts{end+1} = sprintf('BW_%.0fMHz', B/1e6);
    end
    if nargin >= 5 && ~isempty(D) && ~isequal(D, '~')
        name_parts{end+1} = sprintf('D_%.0fkm', D/1e3);
    end

    filename_base = strjoin(name_parts, '_');

    % ==== 核心修改：定位到项目主目录 ====
    current_file_path = mfilename('fullpath');
    [current_folder, ~, ~] = fileparts(current_file_path);
    project_root = fileparts(current_folder);
    folder = fullfile(project_root, 'output_data');
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    % ==== 最大化窗口 ====
    try
        set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        drawnow;
    catch
        warning('⚠️ 无法最大化窗口。');
    end

    % ==== 矢量 PDF + 位图 PNG ====
    filepath_pdf = fullfile(folder, [filename_base, '.pdf']);
    print(fig, filepath_pdf, '-dpdf', '-bestfit');

    filepath_png = fullfile(folder, [filename_base, '.png']);
    saveas(fig, filepath_png);

    fprintf('图像已保存为 PDF：%s\n', filepath_pdf);
    fprintf('图像已保存为 PNG：%s\n', filepath_png);
end
