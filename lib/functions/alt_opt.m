function [W_MMSE, HH_s, HH_real] = alt_opt(setup, channel, csi_mode, sigma_s, sigma_n)
% Unified alternating optimization for SIM and BD architectures, any R layers.
%
% Inputs:
%   setup     - struct with arch, n_layers, sys, W_T, W_T_1
%   channel   - struct with H (and H_est if estimated)
%   csi_mode  - 'perfectcsi' or 'estimated'
%   sigma_s   - signal std
%   sigma_n   - noise std
%
% Outputs:
%   W_MMSE  - MMSE filter
%   HH_s    - effective channel using H (or H_est)
%   HH_real - effective channel using true H (for received signal generation)

    arch     = setup.arch;
    R_layers = setup.n_layers;
    sys      = setup.sys;
    Nc       = sys.Nc;
    W_T      = setup.W_T.';
    W_0      = setup.W_T_1.';
    
    % --- Select working channel ---
    switch lower(csi_mode)
        case 'perfectcsi', H_s = channel.H;
        case 'estimated',  H_s = channel.H_est;
        otherwise, error('csi_mode not listed: %s', csi_mode);
    end
    
    % --- Initialize R_r per architecture ---
    R_cells = cell(R_layers, 1);
    for r = 1:R_layers
        R_cells{r} = sample_layer(arch, Nc);
    end
    
    % --- Build initial cascade G = W_0 * R_1 * W_T * R_2 * ... * R_R ---
    G_SIM = build_AO_cascade(W_0, R_cells, W_T);
    
    Cs = (sigma_s^2) * eye(sys.n_tx);
    Cw = (sigma_n^2) * eye(sys.n_rx);
    
    HH_s   = G_SIM * H_s;
    W_MMSE = Cs * HH_s' / (HH_s * Cs * HH_s' + Cw);
    
    % --- Alternating optimization loop ---
    for iAO = 1:sys.iterAO
        % Backward pass: B_R = H_s, B_r = W_T * R_{r+1} * B_{r+1}
        B_cells = cell(R_layers, 1);
        B_cells{R_layers} = H_s;
        for r = R_layers-1:-1:1
            B_cells{r} = W_T * R_cells{r+1} * B_cells{r+1};
        end
        
        % Forward pass: update each R_r and accumulate A
        A = W_0;
        for r = 1:R_layers
            R_cells{r} = update_layer(arch, W_MMSE, A, B_cells{r}, sys);
            if r < R_layers
                A = A * R_cells{r} * W_T;
            else
                G_SIM = A * R_cells{r};
            end
        end
        
        HH_s   = G_SIM * H_s;
        W_MMSE = Cs * HH_s' / (HH_s * Cs * HH_s' + Cw);
    end
    
    HH_real = G_SIM * channel.H;
end


function R = sample_layer(arch, Nc)
    switch arch
        case "SIM", R = diag(exp(1j * 2*pi * rand(Nc, 1)));
        case "BD",  [R, ~] = qr(randn(Nc) + 1i*randn(Nc));
        otherwise, error('Unknown arch: %s', arch);
    end
end

function G = build_AO_cascade(W_0, R_cells, W_T)
    R_layers = numel(R_cells);
    G = W_0 * R_cells{1};
    for r = 2:R_layers
        G = G * W_T * R_cells{r};
    end
end

function R = update_layer(arch, W_MMSE, A, B, sys)
% Updates layer R_r in the alternating optimization.
%
% Inputs:
%   arch    - "SIM" or "BD"
%   W_MMSE  - current MMSE filter
%   A       - effective channel forward from W_0 up to layer r-1
%   B       - effective channel backward from H_s up to layer r+1
%             (B already encodes channel uncertainty through H_s = H_est)
%   sys     - system struct (used by SIM helper)
%
% Output:
%   R       - updated layer matrix (diagonal phase for SIM, unitary for BD)
%
% Note on csi_mode passed to computeThetaDiagonal2 (SIM case):
%   We call it with 'perfectcsi' even when the outer alternatingOpt is in
%   'estimated' mode. This is intentional: B is already constructed from
%   H_s, which is the *estimated* channel in 'estimated' mode. The estimation
%   error is therefore embedded inside B (and A by extension). At the
%   single-layer update step, we treat A and B as known and compute the
%   optimal diagonal phase — the imperfect-CSI nature lives one level above,
%   not inside this update. This avoids re-implementing the helper and
%   reuses the existing computeThetaDiagonal2 cleanly.

    switch arch
        case "SIM"
            channel_mod = struct();
            channel_mod.H = zeros(sys.n_rx, sys.n_tx);
            channel_mod.G = A;
            channel_mod.F = B;
            % NOTE: 'perfectcsi' here is correct — see header note.
            % A and B already carry the estimation error via H_s.
            R = computeThetaDiagonal2(sys, channel_mod, W_MMSE, 'perfectcsi');
            
        case "BD"
            % BD update via Procrustes (SVD).
            Theta_o = pinv(B * W_MMSE * A);
            [U, ~, V] = svd(Theta_o);
            R = U * V';
            
        otherwise
            error('update_layer: unknown arch %s', arch);
    end
end