%%% Compute LoS channel capacity C --- traditional equal power allocation %%%
function C = los_mimo_capacity(H, P, N0, B)
% fprintf('\n=============================================  START  ===============================================\n');
    
    cond_number = cond(H);
    % fprintf('Trad: Condition number of H: %.3e\n', cond_number);
    
    % H = [1 0; 0 1]; % for debugging

    % Compute capacity
    Nt = 2; % Number of transmit antennas
    SNR = P / (Nt * B * N0); % Signal-to-noise ratio
    

    [sigma, sigma2] = compute_singular_values(H);
    sum_sigma2 = sum(sigma2);
    fprintf('Trad: 奇异值平方之和 = %.6e\n', sum_sigma2);

    % Capacity: C = B * log2(det(I + (P / B N0) * H H^H))
    C = B * sum(log2(1 + SNR * sigma2));

    if abs(sigma2(1) - sigma2(2)) < 1e-13
        fprintf('\n=============================================  START  ===============================================\n');

        fprintf('Trad: H matrix = ');
        disp(H);
        fprintf('Trad: SNR = %.3e \n', SNR);
        fprintf('Singular values of H (traditional): sigma(1) = %.3e, sigma(2) = %.3e\n', sigma(1), sigma(2));
        fprintf('Squared singular values: sigma2(1) = %.3e, sigma2(2) = %.3e;\n', sigma2(1), sigma2(2));

        for i = 1:Nt
            fprintf('Trad: SNR * sigma2_sorted(%d) = %.3e\n',  i, SNR * sigma2(i));
        end

        fprintf('Traditional Capacity: C = %.3e bps\n', C);
    end
    

end
