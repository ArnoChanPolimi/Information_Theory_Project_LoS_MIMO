% 图像绘制函数（根据 is3D 切换）
function draw(a_vals, b_vals, freq_list, freq_labels, Capacity_all, is3D)
    for k = 1:length(freq_list)
        nexttile(k);
        cla; % 清除当前tile内容
        C = Capacity_all{k};
    
        if is3D
            [B_mesh, A_mesh] = meshgrid(b_vals, a_vals);
            surf(B_mesh, A_mesh, C, 'EdgeColor', 'none');
            view(45,30);
        else
            imagesc(b_vals, a_vals, C);
            set(gca, 'YDir', 'normal');
        end
    
        xlabel('b (m)');
        ylabel('a (m)');
        title(['Freq = ', freq_labels{k}]);
        colorbar;
        axis square;
    end
end