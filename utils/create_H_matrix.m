%%%% Channel transfer matrix H %%%%
function [H_scaled] = create_H_matrix(a, b, D, f, Oscillation)
% Compute path distances d_ij

    c = 3e8;
    d11 = D + Oscillation;
    d12 = sqrt((D + Oscillation)^2 + a^2);
    d21 = sqrt(D^2 + b^2);
    d22 = sqrt(D^2 + (a - b)^2);
    lambda = c / f; % Wavelength

    % Compute phase shifts phi_ij
    phi11 = 2*pi*d11/lambda;
    phi12 = 2*pi*d12/lambda;
    phi21 = 2*pi*d21/lambda;
    phi22 = 2*pi*d22/lambda;

    % Construct the channel matrix H
    % H = [ (lambda/d11)*exp(-1j * phi11),   (lambda/d12)*exp(-1j * phi12);
    %       (lambda/d21)*exp(-1j * phi21),   (lambda/d22)*exp(-1j * phi22) ];
    H = [ (lambda/(4*pi*d11))*exp(-1j * phi11),   (lambda/(4*pi*d12))*exp(-1j * phi12);
          (lambda/(4*pi*d21))*exp(-1j * phi21),   (lambda/(4*pi*d22))*exp(-1j * phi22) ];
    % H = [ (1/d11^2)*exp(-1j * phi11),   (1/d12^2)*exp(-1j * phi12);
    %       (1/d21^2)*exp(-1j * phi21),   (1/d22^2)*exp(-1j * phi22) ];
    
    H_scaled = 5e1 * H; %5e-6 * H / norm(H, 'fro');  % 统一能量（最重要）; % Optional scaling (e.g., * 1e2)
    % H_scaled

    % H = [1 0; 0 1]; % for debugging, orthogonal case
    % H = [1 0; 1/sqrt(2) 1/sqrt(2)]; % for debugging, non-orthogonal case
    % fprintf('rank(H) = %d\n', rank(H));
    if rank(H) < 2
        fprintf('WARNING: rank(H) = %d < 2 → Channel is rank-deficient, potential capacity loss.\n', rank(H));
    end

end
