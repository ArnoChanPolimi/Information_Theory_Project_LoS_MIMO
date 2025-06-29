%%% 计算LoS信道容量C %%%
function C = los_mimo_capacity(H, P, N0, B)           
    
    % 计算容量
    Nt = 2; % 发射天线数
    SNR = P/(B*N0); % 信噪比
    
    % 计算容量 C = log2(det(I + (P/N0)*H*H^H))
    C = B * real(log2(det(eye(Nt) + (SNR/Nt)*(H*H'))));
end