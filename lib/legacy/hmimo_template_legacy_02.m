% =========================================================================
% FILE: hmimo_template.m
% DESCRIPTION: Main simulation script for Improved HMIMO channel estimation
% AUTHOR: Roberto C. G. Porto
% EMAIL: camara@ime.eb.br
% GITHUB: https://github.com/rocage/
% DATE: Jun/2026
% =========================================================================

%% 1. Cleanup and Setup
clear; clc; close all;
format short G;

%% 2. Scenario Definition
Pw_dBm_vec       = 3:1:5;       % current TX power dBm
NumberOfPackets  = [10 20 30];  % packets per Pw value
n_layers         = 3;
arch             = "SIM";       % "SIM" or "BD"

if numel(NumberOfPackets) ~= numel(Pw_dBm_vec)
    error('NumberOfPackets and Pw_dBm_vec must have the same length');
end

%% 3. Stages
stage_names = ["Coarse", "Refined", "Perfect CSI"];
nStages = numel(stage_names);

%% 4. Allocation
ber_metric  = zeros(numel(Pw_dBm_vec), nStages);
rsum_metric = zeros(size(ber_metric));
nmse_metric = zeros(size(ber_metric));

total_timer = tic;

%% 5. Main Loop
for idx = 1:numel(Pw_dBm_vec)
    
    % --- Re-init setup with current Pw (ensures sigma_n etc. are updated) ---
    setup = init_system_01(Pw_dBm_vec(idx), n_layers, arch);
    
    % --- temporary, just for sanity (legacy) test ---
    % --- IDD with 3 iterations ---
    setup.idd3       = setup.idd;
    setup.idd3.IDD   = 3;
    setup.idd1       = setup.idd;
    setup.idd1.IDD   = 1;
    [W_L,W_0,~] = generate_WSIM(setup.sys.n_layers,setup.sys.lambda,setup.sys.Nc,setup.sys.n_rx);
    
    % --- Per-packet accumulators ---
    nPackets = NumberOfPackets(idx);
    ber_accum  = zeros(nPackets, nStages);
    rsum_accum = zeros(size(ber_accum));
    nmse_accum = zeros(size(ber_accum));
    
    loop_timer = tic;
    fprintf('\n[Progress] Pw = %d dBm... ', Pw_dBm_vec(idx));
    
    for p = 1:nPackets
        ch = generate_rice(setup);
        [TxSymbolData, Packet] = generate_signal(setup);
        
        out = run_pipeline(setup, ch, TxSymbolData, Packet);
        % out = run_pipeline(setup, ch, TxSymbolData, Packet,W_L,W_0,idx);
        
        ber_accum(p,:)  = [out.coarse.bErrors / setup.sys.DataLength, ...
                           out.refined.bErrors / setup.sys.DataLength, ...
                           NaN];
        rsum_accum(p,:) = [out.coarse.rsum, out.refined.rsum, NaN];
        nmse_accum(p,:) = [out.coarse.nmse, out.refined.nmse, out.perfect.nmse];
    end
    
    ber_metric(idx,:)  = mean(ber_accum);
    rsum_metric(idx,:) = mean(rsum_accum);
    nmse_metric(idx,:) = mean(nmse_accum);
    
    fprintf('done in %.2f s', toc(loop_timer));
end

fprintf('\n\n=== SIMULATION COMPLETED in %.2f seconds ===\n', toc(total_timer));

%% 6. Plots
figure;
semilogy(Pw_dBm_vec, ber_metric, '-o');
legend(stage_names, 'Location', 'best');
xlabel('Pw (dBm)'); ylabel('BER');
grid on; title('BER vs Pw');

figure;
plot(Pw_dBm_vec, rsum_metric, '-o');
legend(stage_names, 'Location', 'best');
xlabel('Pw (dBm)'); ylabel('Sum-rate (bit/s/Hz)');
grid on; title('Sum-rate vs Pw');

figure;
semilogy(Pw_dBm_vec, nmse_metric, '-o');
legend(stage_names, 'Location', 'best');
xlabel('Pw (dBm)'); ylabel('NMSE');
grid on; title('NMSE vs Pw');


% lib/run_pipeline.m
% function out = run_pipeline(setup, ch, TxSymbolData, Packet,W_L,W_0,idx)
function out = run_pipeline(setup, ch, TxSymbolData, Packet)
    out = struct();
    sys = setup.sys;
    estimation = setup.estimation;
    idd3 = setup.idd3;
    idx = 1;
    W_L = setup.W_T;
    W_0 = setup.W_T_1;
    [ch_est,out.coarse.nmse] = CE_SIM_RIS_LMMSE_EP_3layers(sys, ch, estimation, 'initial', TxSymbolData, 0, W_L, W_0, idd3, idx);
    [out.coarse.bErrors, out.coarse.rsum , RecoveredPacket] = IDD_MMSE_SIM_EP(sys, ch_est, estimation, idd3, 'estimated',W_L, W_0,TxSymbolData, Packet, idx);
    
    [ch_est,out.refined.nmse] = CE_SIM_RIS_LMMSE_EP_3layers(sys, ch, estimation, 'improved', TxSymbolData, RecoveredPacket, W_L, W_0, idd3, idx);
    [out.refined.bErrors, out.refined.rsum, ~] = IDD_MMSE_SIM_EP(sys, ch_est, estimation, idd3, 'estimated',W_L, W_0,TxSymbolData, Packet, idx);
    
    [~,out.perfect.nmse] = CE_SIM_RIS_LMMSE_EP_3layers(sys, ch, estimation, 'benchmark', TxSymbolData, 0, W_L, W_0, idd3, idx);
end








