%%%%%%% This function is currently unused %%%%%%%%%%%%%%%%%%
%% After normalization, the total energy of the channel matrix (sum of squares of all elements) is 1.
function [H_norm, energy] = normalize_channel(H_raw)
    % Calculate the Frobenius norm squared (total energy) of the channel matrix
    energy = trace(H_raw * H_raw');
    
    % Normalize the channel matrix
    H_norm = H_raw / sqrt(energy);
    
    % Verify normalized energy
    energy_norm = trace(H_norm * H_norm');
    % fprintf('Channel energy before normalization: %.6e\n', energy);
    % fprintf('Channel energy after normalization: %.6f\n', energy_norm);
end
