function setup = init_system_01(Pw_dBm, n_layers, arch)
% Centralized system initialization for HMIMO channel estimation experiments.
%
% Inputs:
%   Pw_dBm    - transmit power per user (dBm), scalar or vector
%   n_layers  - number of metasurface layers R (default: 1)
%   arch      - "SIM" or "BD" (default: "SIM")
%
% Output:
%   setup     - struct with fields:
%               .sys, .enviroment, .estimation, .idd1, .idd3
%               .arch, .n_layers, .pwr
%               .W_T, .W_T_1
%
% Based on initializeSystem_Extensions_1024_VTC_*_layers (legacy).
% Author: Roberto C. G. Porto — camara@ime.eb.br

    % --- Defaults ---
    if nargin < 2 || isempty(n_layers), n_layers = 1; end
    if nargin < 3 || isempty(arch),     arch     = "SIM"; end
    
    sigma_n2_dbm = -101;
    
    % --- System parameters ---
    QAM_M = 4;
    setup.sys.n_tx        = 4;
    setup.sys.n_rx        = 8;
    setup.sys.Nc          = 64;
    setup.sys.M_ORD       = log2(QAM_M);
    setup.sys.iterAO      = 5;

    % --- Power normalization ---
    setup.sys.sigma_n2_mW = 10^(sigma_n2_dbm/10);
    setup.sys.PWR         = 10^(Pw_dBm/10);

    % --- Active RIS parameters (legacy) ---
    setup.sys.pwr_div   = 1;
    setup.sys.ris_noise = 10^(-100/10);
    
    % --- SIM physical parameters ---
    c  = 3e8;
    f0 = 6e9;
    setup.sys.lambda = c / f0;
    
    % --- Fading and user positions ---
    large_fading_1 = 2.2;
    large_fading_2 = 2.8;
    large_fading_3 = 2.1;
    setup.enviroment.PLs   = @(dis) (10^(-3.73)) ./ (dis.^large_fading_1);
    setup.enviroment.PLw   = @(dis) (10^(-4.12)) ./ (dis.^large_fading_2);
    setup.enviroment.PLsim = @(dis) (10^(-5.20)) ./ (dis.^large_fading_3);
    setup.enviroment.kappa = 2;
    setup.enviroment.pos_ap            = [0; 0; 0];
    setup.enviroment.pos_ris           = [250; 10; 0];
    setup.enviroment.pos_center_users  = [250; 0; 0];
    setup.enviroment.radius_users      = 10;
    
    % --- Channel estimation parameters ---
    setup.estimation.noise_dbm = sigma_n2_dbm;
    setup.estimation.pilots_coarse         = ceil(setup.sys.Nc / setup.sys.n_rx) * setup.sys.n_tx;
    
    % --- LDPC parameters ---
    % NOTE: file is ldpc_mtx_512 despite legacy naming "_1024_"
    load('ldpc_mtx_512.mat', 'A_LDPC', 'G_LDPC');
    [M, N_ldpc] = size(A_LDPC);
    PacketDataLength = N_ldpc - M;
    
    setup.idd.N         = N_ldpc;
    setup.idd.A         = A_LDPC;
    setup.idd.G         = G_LDPC;
    setup.idd.R         = PacketDataLength / N_ldpc;
    setup.idd.iter      = 10;
    setup.idd.IDD       = 3;
    setup.idd.threshold = 2000;
    
    DecimalPacketCodeWord = (0:QAM_M-1)';
    setup.idd.S_matrix = qammod(DecimalPacketCodeWord', QAM_M, 'gray', ...
                                 'UnitAveragePower', true).';
    
    setup.sys.PacketDataLength = PacketDataLength;
    
    % --- Architecture and power index ---
    setup.arch = arch;
    setup.n_layers    = n_layers;
    setup.pwr  = 1;   % index into sys.PWR (1 since we pass Pw_dBm scalar per call)
    
    % --- Metasurface transmission matrices ---
    [setup.W_T, setup.W_T_1, setup.enviroment.Corr_T] = ...
        generate_WSIM(n_layers, setup.sys.lambda, setup.sys.Nc, setup.sys.n_rx);
    
    % --- Derived quantities ---
    setup.sys.DataLength = (setup.sys.PacketDataLength - ...
        setup.estimation.pilots_coarse * setup.sys.M_ORD) * setup.sys.n_tx;

    % --- Derived quantities ---
    % NOTE: sigma_n depends on PWR via this formula. Re-call init_system
    % whenever you change PWR (the main loop already does this).
    setup.sys.sigma_n = sqrt(setup.sys.sigma_n2_mW / (setup.sys.PWR * setup.idd.R));
end