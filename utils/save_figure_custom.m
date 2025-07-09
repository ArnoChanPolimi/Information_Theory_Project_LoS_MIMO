function save_figure_custom(fig, prefix, freq_Hz, B, D, style)
% Save figure as PDF (vector) + PNG (bitmap), for use in Overleaf or papers
%
% Parameters:
%   fig      - figure handle
%   prefix   - filename prefix, e.g., 'CapacityMap'
%   freq_Hz  - frequency (Hz), e.g., 10e9
%   B        - bandwidth (Hz), e.g., 20e6
%   D        - transmission distance (m), e.g., 1000
%   style    - 'vector' (default) or 'bitmap'

    if nargin < 6
        style = 'vector'; % default to vector format
    end

    % ===== Construct filename =====
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

    % ===== Save to output_data folder under project root =====
    current_file_path = mfilename('fullpath');
    [current_folder, ~, ~] = fileparts(current_file_path);
    project_root = fileparts(current_folder);
    folder = fullfile(project_root, 'output_data');
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    % ===== Maximize window and remove GUI elements =====
    try
        set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        set(fig, 'MenuBar', 'none');
        set(fig, 'ToolBar', 'none');
        drawnow;
    catch
        warning('Warning! Failed to maximize window or hide menu.');
    end

    % ===== Save PDF vector image (recommended for papers) =====
    filepath_pdf = fullfile(folder, [filename_base, '.pdf']);
    if strcmpi(style, 'vector')
        exportgraphics(fig, filepath_pdf, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'none');
    else
        % If bitmap style is selected, use print
        print(fig, filepath_pdf, '-dpdf', '-bestfit');
    end

    % ===== Save PNG bitmap =====
    filepath_png = fullfile(folder, [filename_base, '.png']);
    exportgraphics(fig, filepath_png, ...
        'Resolution', 300, ...
        'BackgroundColor', 'white');

    % ===== Output saved paths =====
    fprintf('Figure saved:\n PDF: \033[1m%s\033[0m\n PNG: \033[1m%s\033[0m\n', ...
        filepath_pdf, filepath_png);
end
