% =========================================================================
% FILE: main.m
% DESCRIPTION: Main simulation script for Improved HMIMO channel estimation
% AUTHOR: Roberto C. G. Porto
% EMAIL: camara@ime.eb.br
% GITHUB: https://github.com/rocage/
% DATE: Jun/2026
% =========================================================================

% 1. Cleanup and Setup
clear; clc; close all;
format short G; % Clean numeric format for console output

% 2. Scenario Definition (Default parameters)
Pw_dBm = 3:1:5;  % current TX power dBm 
NumberOfPackets = [10 20 30];      
n_layers = 3;
setup = system_init_concise(Pw_dBm, n_layers);

disp(setup); 
if numel(NumberOfPackets) ~= numel(Pw_dBm), error('Error Packets and Pw_user_dBm vectors must have the same size'); end

% 3. Stages Evaluated in Simulation
stage_names = ["Coarse", "Refined", "Perfect CSI"];
nStages = numel(stage_names);

% 4. Data Allocation and auxiliary variables
ber_metric = zeros(length(Pw_dBm), nStages);
rsum_metric = zeros(size(ber_metric));
nmse_metric = zeros(size(ber_metric));
total_timer = tic;

% 5. Generate stacked Intelligent Metasurfaces
[W_L,W_0,setup.enviroment.Corr_T] = generate_WSIM(setup.sys.n_layers,setup.sys.lambda,setup.sys.Nc,setup.sys.n_rx);

%% 6. Main Loop
for idx = 1:length(Pw_dBm)

    % current parameter
    current_var = Pw_dBm(idx);   
    % update the SNR vector for the next iteration
    setup.TxPower_dBm = current_var;

    % variable Initialization
    ber_accum = zeros(numel(NumberOfPackets), nStages);
    rsum_accum = zeros(size(ber_accum));
    nmse_accum = zeros(size(ber_accum));
    % simulation tracking begin
    loop_timer = tic;
    fprintf('\n[Progress] Simulating: %d... ', current_var);
    % algorithm loop
    parfor p = 1:NumberOfPackets(idx)

        % A. Channel and Signal Generation
        ch = generate_rice(setup);
        % [TxSymbolData, Packet] = generate_bd_ris(setup.sys, setup.estimation, setup.idd3);
        [TxSymbolData, Packet] = generate_signal(setup);

        out = run_pipeline(setup, ch, TxSymbolData, Packet,W_L,W_0,idx);

        % slicing for parfor
        ber_accum(p,:)  = [out.coarse.bErrors/setup.sys.DataLength;  out.refined.bErrors/setup.sys.DataLength;  NaN];
        rsum_accum(p,:) = [out.coarse.rsum; out.refined.rsum; out.perfect.rsum];
        nmse_accum(p,:) = [out.coarse.nmse; out.refined.nmse; NaN];
    end

    % average over packets for each metric
    ber_metric(idx,:) = mean(ber_accum);
    rsum_metric(idx,:) = mean(rsum_accum);
    nmse_metric(idx,:) = mean(nmse_accum);
    
    % simulation tracking end
    elapsed = toc(loop_timer);
    fprintf('Done in %.2f s.', elapsed);
end

total_time = toc(total_timer);
fprintf('\n\n=== SIMULATION COMPLETED in %.2f seconds ===\n', total_time);


semilogy(Pw_dBm, ber_metric)
legend(stage_names)
grid on

figure;
plot(Pw_dBm, rsum_metric)
legend(stage_names)
grid on

figure;
plot(Pw_dBm, nmse_metric)
legend(stage_names)
grid on


function [setup] = system_init_concise(Pw_dBm, n_layers)
    [setup.sys,setup.enviroment,setup.estimation,~,setup.idd3] = initializeSystem_Extensions_1024_VTC_3_layers(Pw_dBm);
    setup.sys.n_layers = n_layers;
    setup.sys.DataLength = (setup.sys.PacketDataLength-setup.estimation.bd*setup.sys.M_ORD)*setup.sys.n_tx;
end

% lib/run_pipeline.m
function out = run_pipeline(setup, ch, TxSymbolData, Packet,W_L,W_0,idx)
    out = struct();
    sys = setup.sys;
    estimation = setup.estimation;
    idd3 = setup.idd3;
    [ch_est,out.coarse.nmse] = CE_SIM_RIS_LMMSE_EP_3layers(sys, ch, estimation, 'initial', TxSymbolData, 0, W_L, W_0, idd3, idx);
    [out.coarse.bErrors, out.coarse.rsum , RecoveredPacket] = IDD_MMSE_SIM_EP(sys, ch_est, estimation, idd3, 'estimated',W_L, W_0,TxSymbolData, Packet, idx);
    
    [ch_est,out.refined.nmse] = CE_SIM_RIS_LMMSE_EP_3layers(sys, ch, estimation, 'improved', TxSymbolData, RecoveredPacket, W_L, W_0, idd3, idx);
    [out.refined.bErrors, out.refined.rsum, ~] = IDD_MMSE_SIM_EP(sys, ch_est, estimation, idd3, 'estimated',W_L, W_0,TxSymbolData, Packet, idx);
    
    [~,out.perfect.nmse] = CE_SIM_RIS_LMMSE_EP_3layers(sys, ch, estimation, 'benchmark', TxSymbolData, 0, W_L, W_0, idd3, idx);
end



% lib/run_pipeline.m
% function out = run_pipeline_template(setup, channel, TxData, Packet)
% out = struct();
% [out.coarse,  channel_est, RecoveredPacket] = run_coarse_stage(setup, channel, TxData, Packet);
% [out.refined, channel_est, RecoveredPacket] = run_refined_stage(setup, channel, TxData, Packet, RecoveredPacket);
% out.perfect                                = run_perfect_stage(setup, channel, TxData);
% end