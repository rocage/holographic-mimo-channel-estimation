% =========================================================================
% FILE: hmimo_template.m
% DESCRIPTION: Main simulation script for HMIMO channel estimation
% (Refactoring)
% AUTHOR: Roberto C. G. Porto
% EMAIL: camara@ime.eb.br
% GITHUB: https://github.com/rocage/holographic-mimo-channel-estimation
% DATE: Jun/2026 (big refactoring)
% =========================================================================

%% 1. Cleanup and Setup
clear; clc; close all;
format short G;

%% 2. Scenario Definition
Pw_dBm_vec       = 3:1:5;       % current TX power dBm
NumberOfPackets  = [10 20 30]/10;  % packets per Pw value

% Pw_dBm_vec = 6:1:10;
% NumberOfPackets = [150 150 200 200 250]/10;
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
    [W_L,W_0,~] = generate_WSIM(setup.n_layers,setup.sys.lambda,setup.sys.Nc,setup.sys.n_rx);
    
    % --- Per-packet accumulators ---
    nPackets = NumberOfPackets(idx);
    ber_accum  = zeros(nPackets, nStages);
    rsum_accum = zeros(size(ber_accum));
    nmse_accum = zeros(size(ber_accum));
    
    loop_timer = tic;
    fprintf('\n[Progress] Pw = %d dBm... ', Pw_dBm_vec(idx));
    
    parfor p = 1:nPackets
        ch = generate_rice(setup);
        [TxSymbolData, Packet] = generate_signal(setup);
        
        out = run_pipeline_hmimo(setup, ch, TxSymbolData, Packet);
        
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

    % --- save partial results ---
    save('results/hmimo_template.mat');
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
plot(Pw_dBm_vec, pow2db(nmse_metric), '-o');
legend(stage_names, 'Location', 'best');
xlabel('Pw (dBm)'); ylabel('NMSE');
grid on; title('NMSE vs Pw');







