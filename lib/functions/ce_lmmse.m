function [ch_est, nmse] = ce_lmmse(setup, channel, pilots_mode, TxSymbolData, RecoveredPacket)
% Unified LMMSE channel estimator for SIM and BD cascades.
%
% Inputs:
%   setup           - struct with sys, estimation, idd, W_T, W_T_1, arch, n_layers, pwr
%   channel         - struct with H, Hkloss, Corr_T
%   pilots_mode     - 'initial' | 'improved' | 'benchmark'
%   TxSymbolData    - transmitted symbols (n_tx × N_sym)
%   RecoveredPacket - only used in 'improved' mode; pass [] otherwise
%
% Outputs:
%   H_est           - estimated channel (Nc × n_tx)
%   nmse            - normalized MSE against channel.H

    % --- Unpack (declares dependencies explicitly) ---
    sys = setup.sys;
    est = setup.estimation;
    idd = setup.idd3;
    W_T   = setup.W_T.';
    W_T_1 = setup.W_T_1.';
    arch  = setup.arch;
    R_cas = setup.n_layers;
    
    H  = channel.H;
    Nc = sys.Nc;
    K  = sys.n_tx;
    M  = sys.n_rx;
    fi = est.pilots_coarse / K;
    M_ORD = sys.M_ORD;
    
    % --- Power normalization (computed once) ---
    sigma2_nest = (10^(est.noise_dbm/10)) / (sys.PWR(setup.pwr) * idd.R);
    sigma2_s    = 1;
    sigma_n     = sqrt(sigma2_nest/2);
    
    % --- Build sensing matrix C per mode ---
    switch lower(pilots_mode)
        case 'initial'
            T0     = est.pilots_coarse;
            T0_eff = T0 / fi;
            pos    = 1 + sys.PacketDataLength/M_ORD;
            pilots = TxSymbolData(:, pos:pos+T0_eff-1) * sqrt(sigma2_s);
            
            % Pre-allocate W_G: each iteration adds M rows
            W_G = zeros(fi*M, Nc);
            for u = 1:fi
                rows = (u-1)*M + (1:M);
                W_G(rows, :) = build_cascade(arch, R_cas, Nc, W_T, W_T_1);
            end
            C = kron(pilots.', W_G);
            
        case 'improved'
            T0 = length(TxSymbolData)/2 + est.pilots_coarse;
            rx_symbols = qammod(RecoveredPacket', 2^M_ORD, ...
                                'InputType','bit', 'UnitAveragePower',true).' * sqrt(sigma2_s);
            pilots = rx_symbols(:, 1:T0) .* exp(1i*pi/4);
            
            % Pre-allocate C: M*T0 rows × Nc*K cols
            C = zeros(M*T0, Nc*K);
            estado_completo = rng();
assignin('base', 'full_state_new', estado_completo.State);
            for u = 1:T0
                aux = build_cascade(arch, R_cas, Nc, W_T, W_T_1);
                rows = (u-1)*M + (1:M);
                C(rows, :) = kron(pilots(:,u).', aux);
            end

        case 'benchmark'
            T0 = length(TxSymbolData)/2 + est.pilots_coarse;
            pilots = TxSymbolData(:, 1:T0) .* exp(1i*pi/4);
            
            C = zeros(M*T0, Nc*K);
            for u = 1:T0
                aux = build_cascade(arch, R_cas, Nc, W_T, W_T_1);
                rows = (u-1)*M + (1:M);
                C(rows, :) = kron(pilots(:,u).', aux);
            end
            
        otherwise
            error('CE_LMMSE:badMode', 'Unknown pilots_mode: %s', pilots_mode);
    end
    
    % --- Observation model ---
    n_obs = size(C, 1);
    noise = sigma_n * (randn(n_obs,1) + 1j*randn(n_obs,1));
    y_p   = C * H(:) + noise;
    
    % --- LMMSE estimation ---
    Rh_vec = kron(diag(channel.Hkloss), channel.Corr_T);
    A_mat  = sigma2_nest * eye(n_obs) + sigma2_s * (C * Rh_vec * C');
    B_mat  = Rh_vec * C' * sqrt(sigma2_s);
    h_est  = B_mat * (A_mat \ y_p);
    
    H_est = reshape(h_est, Nc, K);
    nmse  = norm(H_est - H, 'fro')^2 / norm(H, 'fro')^2;

    ch_est = channel;
    ch_est.H_est = H_est;
end


function aux = build_cascade(arch, n_layers, Nc, W_T, W_T_1)
% Generates one realization: W_T_1 * Theta_1 * W_T * Theta_2 * ... * Theta_R
    aux = W_T_1 * sample_layer(arch, Nc);
    for r = 2:n_layers
        aux = aux * W_T * sample_layer(arch, Nc);
    end
end


function Theta = sample_layer(arch, Nc)
% Random unitary layer per architecture.
    switch arch
        case "SIM"
            Theta = diag(exp(1j * 2*pi * rand(Nc,1)));
        case "BD"
            [Theta, ~] = qr(randn(Nc) + 1i*randn(Nc));
        otherwise
            error('sample_layer:badArch', 'Unknown arch: %s', arch);
    end
end