%%归一化后，信道矩阵的总能量（所有元素平方和）为1。
function [H_norm, energy] = normalize_channel(H_raw)
    % 计算信道矩阵的Frobenius范数平方（总能量）
    energy = trace(H_raw * H_raw');
    
    % 归一化
    H_norm = H_raw / sqrt(energy);
    
    % 验证归一化后的能量
    energy_norm = trace(H_norm * H_norm');
    % fprintf('归一化前信道能量: %.6e\n', energy);
    % fprintf('归一化后信道能量: %.6f\n', energy_norm);
end
