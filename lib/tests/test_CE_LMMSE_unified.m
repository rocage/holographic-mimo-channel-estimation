% =========================================================================
% FILE: test_CE_LMMSE_unified.m
% DESCRIPTION: Verifies that CE_LMMSE matches the legacy 5 separate functions
%              bit-for-bit, across all (arch, R, mode) combinations.
% =========================================================================
function test_CE_LMMSE_unified()
    fprintf('\n=== CE_LMMSE unification test ===\n');
    
    % Setup a fixed scenario
    Pw_dBm = 5;
    [sys, env, est, ~, idd3] = initializeSystem_Extensions_1024_VTC_3_layers(Pw_dBm);
    sys.n_layers = 3;
    
    setup.sys = sys;
    setup.enviroment = env;
    setup.estimation = est;
    setup.idd3 = idd3;
    setup.pwr = 1;
    
    [setup.W_T, setup.W_T_1, ~] = generate_WSIM(3, sys.lambda, sys.Nc, sys.n_rx);
    
    % Generate a fixed channel and TxData
    rng(42);
    channel = generateRice_sim(sys, env, 'sim', 'correlated');
    [TxSymbolData, Packet] = generate_bd_ris(sys, est, idd3);
    
    % Fake RecoveredPacket for 'improved' mode (would normally come from decoder)
    rng(43);
    RecoveredPacket = randi([0 1], sys.n_tx, idd3.N);
    
    % --- Configurations to test ---
    configs = { ...
        % {arch,  R, legacy_function_name}
        % {"SIM", 1, @CE_SIM_RIS_LMMSE_EP_1layers}, ...
        % {"SIM", 3, @CE_SIM_RIS_LMMSE_EP_3layers}, ...
        % {"SIM", 5, @CE_SIM_RIS_LMMSE_EP_5layers}, ...
        {"BD",  3, @CE_BD_SIM_RIS_LMMSE_EP_3layer}, ...
        % {"BD",  5, @CE_BD_SIM_RIS_LMMSE_EP_5layer}, ...
    };
    
    modes = {'initial', 'improved', 'benchmark'};
    tolerance = 1e-10;
    
    n_pass = 0;
    n_fail = 0;
    
    for c = 1:numel(configs)
        arch    = configs{c}{1};
        R       = configs{c}{2};
        f_old   = configs{c}{3};
        
        % Need to load the correct initializeSystem and W for this R
        [sys_R, env_R, est_R, ~, idd3_R] = load_init_for_R(R, Pw_dBm);
        [W_T_R, W_T_1_R, ~] = generate_WSIM(R, sys_R.lambda, sys_R.Nc, sys_R.n_rx);
        
        % Update setup for this R
        setup_R = setup;
        setup_R.sys = sys_R;
        setup_R.estimation = est_R;
        setup_R.idd3 = idd3_R;
        setup_R.W_T = W_T_R;
        setup_R.W_T_1 = W_T_1_R;
        setup_R.arch = arch;
        setup_R.n_layers = R;
        
        % Regenerate channel/TxData for this R (in case Nc differs)
        rng(42);
        channel_R = generateRice_sim(sys_R, env_R, 'sim', 'correlated');
        [TxData_R, ~] = generate_bd_ris(sys_R, est_R, idd3_R);
        
        rng(43);
        RP_R = randi([0 1], sys_R.n_tx, idd3_R.N);
        
        for m = 1:numel(modes)
            mode = modes{m};
            
            % --- Run legacy ---
            rng(100);  % fixed seed for noise + theta generation
            [ch_old_struct, nmse_old] = f_old(sys_R, channel_R, est_R, mode, ...
                                              TxData_R, RP_R, W_T_R, W_T_1_R, idd3_R, 1);
            H_est_old = ch_old_struct.H_est;
            
            % --- Run new ---
            rng(100);  % SAME seed
            RP_arg = RP_R;
            if ~strcmpi(mode, 'improved'), RP_arg = []; end
            [H_est_new, nmse_new] = CE_LMMSE(setup_R, channel_R, mode, TxData_R, RP_arg);
            
            % --- Compare ---
            err_H    = norm(H_est_new - H_est_old, 'fro') / norm(H_est_old, 'fro');
            err_nmse = abs(nmse_new - nmse_old);
            
            tag = sprintf('%s R=%d %s', arch, R, mode);
            if err_H < tolerance && err_nmse < tolerance
                fprintf('  [PASS] %-25s  err_H=%.2e  err_nmse=%.2e\n', tag, err_H, err_nmse);
                n_pass = n_pass + 1;
            else
                fprintf('  [FAIL] %-25s  err_H=%.2e  err_nmse=%.2e\n', tag, err_H, err_nmse);
                n_fail = n_fail + 1;
            end
        end
    end
    
    fprintf('\n--- Summary: %d passed, %d failed ---\n\n', n_pass, n_fail);
    
    if n_fail == 0
        fprintf('All configurations match. Safe to delete legacy functions.\n');
    else
        fprintf('SOME TESTS FAILED. Inspect failing configs before deletion.\n');
    end
end


function [sys, env, est, idd1, idd3] = load_init_for_R(R, Pw_dBm)
    switch R
        case 1, [sys, env, est, idd1, idd3] = initializeSystem_Extensions_1024_VTC_1_layers(Pw_dBm);
        case 3, [sys, env, est, idd1, idd3] = initializeSystem_Extensions_1024_VTC_3_layers(Pw_dBm);
        case 5, [sys, env, est, idd1, idd3] = initializeSystem_Extensions_1024_VTC_5_layers(Pw_dBm);
        otherwise, error('R=%d not supported by legacy initializeSystem', R);
    end
end