%%% Compute the singular values of the channel matrix H %%%
function [sigma, sigma2] = compute_singular_values(H)
    % compute_singular_values computes the singular values of a given matrix H
    %
    % Input:
    %   H - Channel matrix (complex-valued)
    %
    % Output:
    %   sigma  - Vector of singular values (in descending order)
    %   sigma2 - Squared singular values

    % Compute singular values
    [~, S, ~] = svd(H, 'econ');  % Debug: compute singular values of H
    sigma = diag(S);

    % Compute squared singular values
    sigma2 = sigma.^2;
end

