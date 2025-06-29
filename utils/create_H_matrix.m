%%%% 信道转移矩阵 H %%%%%
function H = create_H_matrix(a, b, D, f, Oscillation)
% 计算路径距离 d_ij
    
    c = 3e8;
    d11 = D + Oscillation;
    d12 = sqrt((D + Oscillation)^2 + a^2);
    d21 = sqrt(D^2 + b^2);
    d22 = sqrt(D^2 + (a - b)^2);
    lambda = c / f; % 波长
    
    % 计算相位 phi_ij
    phi11 = 2*pi*d11/lambda;
    phi12 = 2*pi*d12/lambda;
    phi21 = 2*pi*d21/lambda;
    phi22 = 2*pi*d22/lambda;
    
    % 构造信道矩阵 H
    H = [ (lambda/d11)*exp(-1j * phi11),   (lambda/d12)*exp(-1j * phi12);
          (lambda/d21)*exp(-1j * phi21),   (lambda/d22)*exp(-1j * phi22) ];
end