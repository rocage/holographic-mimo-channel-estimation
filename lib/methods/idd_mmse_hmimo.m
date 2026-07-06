function [bErrors, Rsum_aux, RecoveredPacket] = idd_mmse_hmimo(setup, ch, csi_mode, TxSymbolData, Packet)
% Iterative Detection and Decoding (IDD) with MMSE-SIC detector.
% Unified for SIM and BD architectures, any number of layers R >= 1.
%
% Uses:
%   - alt_opt   : alternating optimization to design the metasurface layers
%                 (dispatches SIM/BD internally via setup.arch)
%   - compute_wk: linear MMSE filter with soft interference cancellation
%   - decodeLogDomain_box_plus_quantizer2: LDPC sum-product decoder
%
% Inputs:
%   setup        - struct from init_system (sys, estimation, idd, W_T, W_T_1,
%                  arch ∈ {"SIM","BD"}, n_layers >= 1)
%   ch           - channel struct (H, and H_est if csi_mode = 'estimated')
%   csi_mode     - 'perfectcsi' or 'estimated'
%   TxSymbolData - transmitted symbols (n_tx x N_sym)
%   Packet       - original bits (for BER computation)
%
% Outputs:
%   bErrors         - number of bit errors on data payload
%   Rsum_aux        - sum-rate approximation from post-SIC SNIR
%   RecoveredPacket - decoded packet bits (fed back to the pipeline for refinement)
%
% Author: Roberto C. G. Porto — camara@ime.eb.br
% Last update: Jun/2026 (unified SIM/BD, any R)

% --- Configuration by CSI mode ---
switch lower(csi_mode)
    case 'perfectcsi'
        refine_bool = false;
        data_start  = 0;
    case 'estimated'
        refine_bool = true;
        data_start  = setup.estimation.pilots_coarse * setup.sys.M_ORD;
    otherwise
        error('idd_mmse: unknown csi_mode "%s"', csi_mode);
end

% --- Unpack (declares dependencies explicitly) ---
sys       = setup.sys;
est       = setup.estimation;
idd       = setup.idd;

n_tx      = sys.n_tx;
n_rx      = sys.n_rx;
M_ORD     = sys.M_ORD;
N         = idd.N;
IDD       = idd.IDD;
iter      = idd.iter;
threshold = idd.threshold;
A_LDPC    = idd.A;
S_matrix  = idd.S_matrix;
Ex        = mean(abs(S_matrix).^2);

% --- Power normalization ---
sigma_n = sqrt(sys.sigma_n2_mW / (sys.PWR * idd.R));
sigma_s = 1;

% --- Preallocation ---
La              = zeros(n_tx, N);
Lc              = zeros(n_tx, N);
Ld              = zeros(n_tx, N);
RecoveredData   = zeros(n_tx, sys.PacketDataLength);
RecoveredPacket = zeros(n_tx, N);
data_flag       = zeros(n_tx, 1);
Detectedsymbols = zeros(n_tx, N/M_ORD);
gamma_k         = zeros(n_tx, 1);

% --- Metasurface design via alternating optimization ---
% alt_opt dispatches SIM vs BD internally via setup.arch, and handles
% any R >= 1 (BD-RIS classic for R=1, SIM/BD-SIM for R > 1).
[~, HH_est, HH_real] = alt_opt(setup, ch, csi_mode, sigma_s, sigma_n);

% --- Received signal (uses HH_real: optimized cascade × true channel) ---
noise = (sigma_n/sqrt(2)) * (randn(n_rx, length(TxSymbolData)) + 1i*randn(n_rx, length(TxSymbolData)));
ReceivedSymbols = HH_real * TxSymbolData + noise;

% --- IDD loop ---
for IDD_counter = 1:IDD
    
    % Pilot refinement: fix LLRs at known pilot positions
    if refine_bool
        Lc = encodedPilotRefinement3(sys, est.pilots_coarse, Lc, Packet, 2000);
    end
    
    % --- Detection: soft MMSE-SIC per symbol ---
    for ii = 1:N/M_ORD
        Lc_Sic  = Lc(:, M_ORD*ii - (M_ORD-1) : M_ORD*ii);
        P_prior = num2mtx(sys, Lc_Sic);
        
        % Soft symbol mean and variance from LDPC extrinsic info
        x_soft = round(P_prior * S_matrix, 10);
        s_qam  = zeros(n_tx, M_ORD^2);
        for jj = 1:M_ORD^2
            s_qam(:, jj) = (abs(S_matrix(jj) - x_soft).^2) .* P_prior(:, jj);
        end
        sigma2_vec = sum(s_qam, 2);
        
        % Linear MMSE filter with soft interference cancellation
        Wk = compute_wk(HH_est, n_tx, n_rx, sigma2_vec, Ex, sigma_s, sigma_n);
        
        % Per-user detection
        for k = 1:n_tx
            wk = Wk(k, :)';
            
            % Cancel soft interference from other users
            CancelSoftInfo = ReceivedSymbols(:, ii) ...
                - sum(repmat(x_soft.', size(HH_est, 1), 1) .* HH_est, 2) ...
                + x_soft(k) * HH_est(:, k);
            
            Detectedsymbols(k, ii) = wk' * CancelSoftInfo;
            
            % Effective SNR (Miuk) and residual noise variance (Itak^2)
            Miuk = wk' * HH_est(:, k);
            Itak = sqrt(sigma_s^2 * (Miuk - Miuk^2));
            
            % Symbol-level likelihood P(y|x)
            PL = 1/pi/Itak^2 * exp(-1/Itak^2 * abs(Detectedsymbols(k, ii) - Miuk .* S_matrix).^2);
            
            % Extrinsic LLR per coded bit
            for jj = 1:M_ORD
                [pos1, pos0] = llr_pl_pr(sys, jj);
                aux = real(P_prior(k, pos1) * PL(pos1)) / real(P_prior(k, pos0) * PL(pos0));
                Ld(k, M_ORD*(ii-1) + jj) = log(aux) - Lc_Sic(k, jj);
            end
            
            % SNIR estimate (DOI: 10.1109/26.774855)
            gamma_k(k) = 1 * inv(inv(Miuk) - 1);
        end
    end
    
    % LLR clipping (prevents saturation in LDPC decoder)
    Ld(Ld >  threshold) =  threshold;
    Ld(Ld < -threshold) = -threshold;
    
    % --- Decoding: LDPC sum-product per user ---
    for i = 1:n_tx
        if data_flag(i) == 0
            [La(i, :), Ld(i, :)] = decodeLogDomain_box_plus_quantizer2(A_LDPC, iter, -Ld(i, :));
            Lc(i, :) = La(i, :) - Ld(i, :);
            
            % Hard decisions
            RecoveredData(i, :)   = double(La(i, N/2+1:end) > 0);
            RecoveredPacket(i, :) = double(La(i, :) > 0);
            
            % Early stop for this user if packet matches
            if sum(mod(RecoveredData(i, :) + Packet(i, :), 2)) == 0
                data_flag(i, :) = 1;
            end
        end
    end
    
    % Global early stop when all users converged
    if sum(data_flag - ones(n_tx, 1)) == 0
        break;
    end
end

% --- Metrics ---
bit_diff = RecoveredData(:, 1+data_start:end) - Packet(:, 1+data_start:end);
bErrors  = sum(abs(bit_diff(:)));
Rsum_aux = sum(log2(gamma_k + 1));

end