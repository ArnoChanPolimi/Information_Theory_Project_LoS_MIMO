% %% %% 稳健的保存图像函数
% %%  使用示例(灵活省略参数):
% % save_figure_custom(fig, 'Capacity', 80e9, 120e6, 2000);     % 全部有
% % save_figure_custom(fig, 'Capacity', 80e9, '~', 2000);       % 跳过带宽
% % save_figure_custom(fig, 'Capacity', '~', '~', 5000);        % 只有D
% % save_figure_custom(fig, 'Capacity', '~', '~', '~');         % 只有前缀
% % 保存结果文件名示例：
% % Capacity_Freq_80GHz_BW_120MHz_D_2km.png
% % 
% % Capacity_Freq_80GHz_D_2km.png
% % 
% % Capacity_D_5km.png
% % 
% % Capacity.png
% %%
% function save_figure_custom(fig, prefix, freq_Hz, B, D)
%     name_parts = {prefix};
% 
%     if nargin >= 3 && ~isempty(freq_Hz) && ~isequal(freq_Hz, '~')
%         name_parts{end+1} = sprintf('Freq_%.0fGHz', freq_Hz/1e9);
%     end
%     if nargin >= 4 && ~isempty(B) && ~isequal(B, '~')
%         name_parts{end+1} = sprintf('BW_%.0fMHz', B/1e6);
%     end
%     if nargin >= 5 && ~isempty(D) && ~isequal(D, '~')
%         name_parts{end+1} = sprintf('D_%.0fkm', D/1e3);
%     end
% 
%     filename_base = strjoin(name_parts, '_');
% 
%     % ==== 核心修改：定位到项目主目录 ====
%     current_file_path = mfilename('fullpath');
%     [current_folder, ~, ~] = fileparts(current_file_path);
%     project_root = fileparts(current_folder);
%     folder = fullfile(project_root, 'output_data');
%     if ~exist(folder, 'dir')
%         mkdir(folder);
%     end
% 
%     % ==== 最大化窗口 ====
%     try
%         set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
%         drawnow;
%     catch
%         warning('⚠️ 无法最大化窗口。');
%     end
% 
%     % ==== 矢量 PDF + 位图 PNG ====
%     filepath_pdf = fullfile(folder, [filename_base, '.pdf']);
%     print(fig, filepath_pdf, '-dpdf', '-bestfit');
% 
%     filepath_png = fullfile(folder, [filename_base, '.png']);
%     saveas(fig, filepath_png);
% 
%     fprintf('图像已保存为 PDF：%s\n', filepath_pdf);
%     fprintf('图像已保存为 PNG：%s\n', filepath_png);
% end
function save_figure_custom(fig, prefix, freq_Hz, B, D, style)
% 保存图像为 PDF（矢量）+ PNG（位图），用于 Overleaf 或论文插图
%
% 参数说明：
%   fig      - 图句柄
%   prefix   - 文件名前缀，例如 'CapacityMap'
%   freq_Hz  - 频率（Hz），例如 10e9
%   B        - 带宽（Hz），例如 20e6
%   D        - 传输距离（m），例如 1000
%   style    - 'vector'（默认）或 'bitmap'

    if nargin < 6
        style = 'vector'; % 默认保存为矢量图
    end

    % ===== 构造文件名 =====
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

    % ===== 项目根目录下的 output_data 文件夹 =====
    current_file_path = mfilename('fullpath');
    [current_folder, ~, ~] = fileparts(current_file_path);
    project_root = fileparts(current_folder);
    folder = fullfile(project_root, 'output_data');
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    % ===== 最大化窗口并去除多余元素 =====
    try
        set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        set(fig, 'MenuBar', 'none');
        set(fig, 'ToolBar', 'none');
        drawnow;
    catch
        warning('⚠️ 无法最大化窗口或隐藏菜单。');
    end

    % ===== 保存 PDF 矢量图（推荐用于论文） =====
    filepath_pdf = fullfile(folder, [filename_base, '.pdf']);
    if strcmpi(style, 'vector')
        exportgraphics(fig, filepath_pdf, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'none');
    else
        % 若选为 bitmap，则使用 print
        print(fig, filepath_pdf, '-dpdf', '-bestfit');
    end

    % ===== 保存 PNG 位图 =====
    filepath_png = fullfile(folder, [filename_base, '.png']);
    exportgraphics(fig, filepath_png, ...
        'Resolution', 300, ...
        'BackgroundColor', 'white');

    % ===== 输出路径 =====
    fprintf('[✔] 图像已保存：\n📄 PDF: \033[1m%s\033[0m\n🖼️ PNG: \033[1m%s\033[0m\n', ...
        filepath_pdf, filepath_png);
end
