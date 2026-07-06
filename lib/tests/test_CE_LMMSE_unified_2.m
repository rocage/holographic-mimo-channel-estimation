% =========================================================================
% FILE: test_CE_LMMSE_statistical.m
% DESCRIPTION: Validates that CE_LMMSE is statistically equivalent to legacy
%              implementations across N independent trials. Does NOT require
%              bit-identical rand consumption.
% =========================================================================

clear; clc;

% --- Configuração do teste ---
N_trials = 100;
Pw_dBm   = 5;
modes    = {'initial', 'improved', 'benchmark'};

% Configs a testar
configs = { ...
    % {arch,  R, legacy_function_handle, init_function_handle}
    % {"SIM", 1, @CE_SIM_RIS_LMMSE_EP_1layers,  @initializeSystem_Extensions_1024_VTC_1_layers}, ...
    % {"SIM", 3, @CE_SIM_RIS_LMMSE_EP_3layers,  @initializeSystem_Extensions_1024_VTC_3_layers}, ...
    % {"SIM", 5, @CE_SIM_RIS_LMMSE_EP_5layers,  @initializeSystem_Extensions_1024_VTC_5_layers}, ...
    {"BD",  3, @CE_BD_SIM_RIS_LMMSE_EP_3layer, @initializeSystem_Extensions_1024_VTC_3_layers}, ...
    {"BD",  5, @CE_BD_SIM_RIS_LMMSE_EP_5layer, @initializeSystem_Extensions_1024_VTC_5_layers}, ...
};

fprintf('\n=== CE_LMMSE statistical equivalence test ===\n');
fprintf('N_trials = %d per (arch, R, mode) combination\n\n', N_trials);

for c = 1:numel(configs)
    arch    = configs{c}{1};
    R       = configs{c}{2};
    f_old   = configs{c}{3};
    f_init  = configs{c}{4};
    
    % --- Init para esta config ---
    [sys, env, est, ~, idd3] = f_init(Pw_dBm);
    sys.n_layers = R;
    
    setup.sys        = sys;
    setup.enviroment = env;
    setup.estimation = est;
    setup.idd3       = idd3;
    setup.arch       = arch;
    setup.n_layers   = R;
    setup.pwr        = 1;
    
    [W_T, W_T_1, ~] = generate_WSIM(R, sys.lambda, sys.Nc, sys.n_rx);
    setup.W_T   = W_T;
    setup.W_T_1 = W_T_1;
    
    fprintf('--- %s R=%d ---\n', arch, R);
    
    for m = 1:numel(modes)
        mode = modes{m};
        
        nmse_old_vec = zeros(N_trials, 1);
        nmse_new_vec = zeros(N_trials, 1);
        
        parfor  trial = 1:N_trials
            % Mesmo canal e dados em cada trial (compara em condições iguais)
            rng(trial);
            channel = generateRice_sim(sys, env, 'sim', 'correlated');
            [TxData, ~] = generate_bd_ris(sys, est, idd3);
            
            % RecoveredPacket pra modo 'improved' (canal já decodificado)
            if strcmpi(mode, 'improved')
                rng(trial + 10000);   % seed diferente para evitar correlação
                RP = randi([0 1], sys.n_tx, idd3.N);
            else
                RP = [];
            end
            
            % --- Versão legacy ---
            rng(trial + 20000);  % seed específica para esta trial
            [ch_old_struct, nmse_old_vec(trial)] = f_old( ...
                sys, channel, est, mode, TxData, RP, W_T, W_T_1, idd3, 1);
            
            % --- Versão nova ---
            rng(trial + 30000);  % seed diferente — não esperamos bit-identidade
            [~, nmse_new_vec(trial)] = CE_LMMSE(setup, channel, mode, TxData, RP);
        end
        
        % --- Estatísticas ---
        mean_old = mean(nmse_old_vec);
        mean_new = mean(nmse_new_vec);
        std_old  = std(nmse_old_vec);
        std_new  = std(nmse_new_vec);
        
        % Erro relativo entre médias
        rel_diff = abs(mean_old - mean_new) / mean_old;
        
        % t-test: H0 = médias iguais
        [~, p_val] = ttest2(nmse_old_vec, nmse_new_vec);
        
        % Veredito: passou se médias dentro de erro estatístico
        % Threshold conservador: diferença das médias < std/sqrt(N) * 3
        threshold = 3 * std_old / sqrt(N_trials);
        passed = (abs(mean_old - mean_new) < threshold) || (p_val > 0.05);
        
        if passed
            status = '[PASS]';
        else
            status = '[FAIL]';
        end
        
        fprintf('  %s %-11s  old=%.4e±%.4e  new=%.4e±%.4e  rel_diff=%.2f%%  p=%.3f\n', ...
                status, mode, mean_old, std_old, mean_new, std_new, ...
                100*rel_diff, p_val);
    end
    fprintf('\n');
end

fprintf('=== Done ===\n');