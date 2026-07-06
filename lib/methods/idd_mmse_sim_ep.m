function [bErrors, Rsum_aux, RecoveredPacket] = idd_mmse_sim_ep(setup, ch, csi_mode,TxSymbolData, Packet)
% (sys, ch, estimation, idd, csi_mode,W_L, W_0,TxSymbolData, Packet, pwr)

switch lower(csi_mode)
    case 'perfectcsi'
        refine_bool = false;
        data_start = 0;
    case 'estimated'
        refine_bool = true;
        data_start = setup.estimation.pilots_coarse*setup.sys.M_ORD;
    otherwise
        error('Method not listed: %s', csi_mode);
end

%Init parameters for the method -------------------------
PacketDataLength = setup.sys.PacketDataLength;
n_tx = setup.sys.n_tx;
n_rx = setup.sys.n_rx;
M_ORD = setup.sys.M_ORD;
sigma_n2_mW = setup.sys.sigma_n2_mW;
PWR = setup.sys.PWR;
IDD = setup.idd.IDD;
iter = setup.idd.iter;
% G_LDPC = setup.idd.G;
A_LDPC = setup.idd.A;
S_matrix = setup.idd.S_matrix;
R = setup.idd.R;
N = setup.idd.N;
Ex = mean(abs(S_matrix).^2);
threshold = setup.idd.threshold;
% Power Normalization -----------------------------------
sigma_n = sqrt(sigma_n2_mW/(PWR*R));
sigma_s = 1;
% IID scheme -----------------------------------
La = zeros(n_tx,N);
Lc = zeros(n_tx,N);
Ld = zeros(n_tx,N);
RecoveredData = zeros(n_tx,PacketDataLength);
RecoveredPacket = zeros(n_tx,N);
data_flag = zeros(n_tx,1);
Detectedsymbols = zeros(n_tx,N/M_ORD);
gamma_k=zeros(n_tx,1);

W_L = setup.W_T;
W_0 = setup.W_T_1;

W_L = W_L.';
W_0 = W_0.';

% Alternate Optimization ------------------------------------------------------------------
        switch lower(setup.n_layers)
            case 1
                [~, HH_est, HH_real] = alternatingOpt_1(setup.sys,ch,csi_mode,sigma_s,sigma_n,W_L, W_0);
            case 3
                [~, HH_est, HH_real] = alternatingOpt_3(setup.sys,ch,csi_mode,sigma_s,sigma_n,W_L, W_0);
            case 5
                [~, HH_est, HH_real] = alternatingOpt_5(setup.sys,ch,csi_mode,sigma_s,sigma_n,W_L, W_0);
            otherwise
                error('The number of layers selected in sys.n_layers is not valid, check the initfunction');
        end

% Channel ----------------------------------------------------------------------------------
noise = (sigma_n/sqrt(2))*(randn(n_rx,length(TxSymbolData))+1i*randn(n_rx,length(TxSymbolData)));
ReceivedSymbols = HH_real*TxSymbolData + noise; % TENHO
% QUE AJEITAR AQUI

for IDD_counter = 1:IDD
    if refine_bool
        Lc = encodedPilotRefinement3(setup.sys, setup.estimation.pilots_coarse, Lc, Packet,2000);
    end
    for ii = 1:N/M_ORD
        % evaluate the a-priori probability:
        Lc_Sic = Lc(:,M_ORD*ii-(M_ORD-1):M_ORD*ii);
        P_prior = num2mtx(setup.sys,Lc_Sic);

        % symbol mean
        x_soft = round(P_prior*S_matrix,10);

        % symbol variance
        s_qam = zeros (n_tx,M_ORD^2);
        for jj = 1:M_ORD^2
            s_qam(:,jj) = (abs(S_matrix(jj)-x_soft).^2).*P_prior(:,jj);
        end
        sigma2_vec = sum(s_qam,2);

        % Filter design will depend on the estimation mode
        Wk = compute_wk(HH_est,n_tx,n_rx,sigma2_vec,Ex,sigma_s,sigma_n);

        for k = 1:n_tx
            % select the linear filter for the user 'k'
            wk = Wk(k,:)';

            % cancel the soft estimates
            CancelSoftInfo = ReceivedSymbols(:,ii)-sum(repmat(x_soft.',size(HH_est,1),1).*HH_est,2) + x_soft(k)*HH_est(:,k) ;

            % detection estimate
            Detectedsymbols(k,ii) = wk'*CancelSoftInfo;

            % compute Miuk and Itak
            Miuk = wk'*HH_est(:,k);
            Itak = sqrt(sigma_s^2*(Miuk-Miuk^2));

            % evaluate the likelihood function P(x_est|x)
            PL = 1/pi/Itak^2*exp(-1/Itak^2*abs(Detectedsymbols(k,ii)-Miuk.*S_matrix).^2);

            % LLR value, for user 'k'
            for jj = 1:M_ORD
                [pos1,pos0] = llr_pl_pr(setup.sys,jj);
                aux = real(P_prior(k,pos1)*PL(pos1))/real(P_prior(k,pos0)*PL(pos0));
                Ld(k,M_ORD*(ii-1)+jj) = log(aux)- Lc_Sic(k,jj);
            end

            % SNIR computation - formula from DOI: 10.1109/26.774855
            gamma_k(k) = 1*inv(inv(Miuk)-1);
        end
    end

    Ld(Ld > threshold) =   threshold;
    Ld(Ld < -threshold) = -threshold;

    for i=1:n_tx
        if data_flag(i) == 0
            [La(i,:),Ld(i,:)] = decodeLogDomain_box_plus_quantizer2(A_LDPC,iter,-Ld(i,:));
            Lc(i,:) = La(i,:)-Ld(i,:);

            % Hard decision
            RecoveredData(i,:) = double(La(i,N/2+1:end)>0);
            RecoveredPacket(i,:) = double(La(i,:)>0);
            if sum(mod(RecoveredData(i,:)+Packet(i,:),2)) == 0
                data_flag(i,:) = 1;
            end
        end
    end

    if sum(data_flag-ones(n_tx,1)) == 0
        break;
    end

end

bit_diff = RecoveredData(:,1+data_start:end)-Packet(:,1+data_start:end);
bErrors = sum(abs(bit_diff(:)));

%  Sum-Rate  ---------------------------------------------------------------
Rsum_aux = sum(log2(gamma_k+1));

end

function [W_MMSE, HH_s, HH_real] = alternatingOpt_1(sys,channel,csi_mode,sigma_s,sigma_n,W_T, W_0)

%Init channel (if it is estimated or not)----------------
switch lower(csi_mode)
    case 'perfectcsi'
        H_s = channel.H;
    case 'estimated'
        H_s = channel.H_est;
    otherwise
        error('Method not listed: %s', csi_mode);
end

%Init channel for the method -------------------------
R1=diag(exp(1j*2*pi*rand(sys.Nc,1)));

G_SIM = (W_0) * R1; 

Cs = (sigma_s^2)*eye(sys.n_tx);
Cw = (sigma_n^2)*eye(sys.n_rx);

HH_s = G_SIM*H_s;

% computation of the linear filter -------------------------
W_MMSE = Cs*HH_s'/(HH_s*Cs*HH_s'+Cw);


channel_mod = channel;
sys_mod = sys;
% sys_mod.n_ = sys.Nc;
% alternating optmization -------------------------
for i = 1:sys.iterAO
    % calculation of the theta - IMPLEMENTED JUST FOR DIAGONAL
    A = W_0;
    B1 = H_s;

    channel_mod.H = zeros(sys.n_rx,sys.n_tx);
    channel_mod.G = A;
    channel_mod.F = B1;
    [R1] = computeThetaDiagonal2(sys_mod,channel_mod,W_MMSE,'perfectcsi');

    G_SIM = A*R1;

    HH_s = G_SIM*B1;

    % calculation of the linear filter
    W_MMSE = Cs*HH_s'/(HH_s*Cs*HH_s'+Cw);
end

    HH_real = G_SIM*channel.H;
end

function [W_MMSE, HH_s, HH_real] = alternatingOpt_3(sys,channel,csi_mode,sigma_s,sigma_n,W_T, W_0)

%Init channel (if it is estimated or not)----------------
switch lower(csi_mode)
    case 'perfectcsi'
        H_s = channel.H;
    case 'estimated'
        H_s = channel.H_est;
    otherwise
        error('Method not listed: %s', csi_mode);
end

%Init channel for the method -------------------------
R1=diag(exp(1j*2*pi*rand(sys.Nc,1)));
R2=diag(exp(1j*2*pi*rand(sys.Nc,1)));
R3=diag(exp(1j*2*pi*rand(sys.Nc,1)));

G_SIM = (W_0) * R1 * W_T * R2 * W_T * R3;

Cs = (sigma_s^2)*eye(sys.n_tx);
Cw = (sigma_n^2)*eye(sys.n_rx);

HH_s = G_SIM*H_s;

% computation of the linear filter -------------------------
W_MMSE = Cs*HH_s'/(HH_s*Cs*HH_s'+Cw);


channel_mod = channel;
sys_mod = sys;
% sys_mod.n_ = sys.Nc;
% alternating optmization -------------------------
for i = 1:sys.iterAO
    % calculation of the theta - IMPLEMENTED JUST FOR DIAGONAL
    A = W_0;
    B3 = H_s;
    B2 = W_T*R3*B3;
    B1 = W_T*R2*B2;

    channel_mod.H = zeros(sys.n_rx,sys.n_tx);
    channel_mod.G = A;
    channel_mod.F = B1;
    [R1] = computeThetaDiagonal2(sys_mod,channel_mod,W_MMSE,'perfectcsi');

    A = A*R1*W_T;
    channel_mod.G = A;
    channel_mod.F = B2;
    [R2] = computeThetaDiagonal2(sys_mod,channel_mod,W_MMSE,'perfectcsi');

    A = A*R2*W_T;
    channel_mod.G = A;
    channel_mod.F = B3;
    [R3] = computeThetaDiagonal2(sys_mod,channel_mod,W_MMSE,'perfectcsi');

    % channel update
    G_SIM = A*R3;

    HH_s = G_SIM*B3;

    % calculation of the linear filter
    W_MMSE = Cs*HH_s'/(HH_s*Cs*HH_s'+Cw);
end

    HH_real = G_SIM*channel.H;
end

function [W_MMSE, HH_s, HH_real] = alternatingOpt_5(sys,channel,csi_mode,sigma_s,sigma_n,W_T, W_0)

%Init channel (if it is estimated or not)----------------
switch lower(csi_mode)
    case 'perfectcsi'
        H_s = channel.H;
    case 'estimated'
        H_s = channel.H_est;
    otherwise
        error('Method not listed: %s', csi_mode);
end

%Init channel for the method -------------------------
R1=diag(exp(1j*2*pi*rand(sys.Nc,1)));
R2=diag(exp(1j*2*pi*rand(sys.Nc,1)));
R3=diag(exp(1j*2*pi*rand(sys.Nc,1)));
R4=diag(exp(1j*2*pi*rand(sys.Nc,1)));
R5=diag(exp(1j*2*pi*rand(sys.Nc,1)));

G_SIM = (W_0) * R1 * W_T * R2 * W_T * R3* W_T * R4* W_T * R5;

Cs = (sigma_s^2)*eye(sys.n_tx);
Cw = (sigma_n^2)*eye(sys.n_rx);

HH_s = G_SIM*H_s;

% computation of the linear filter -------------------------
W_MMSE = Cs*HH_s'/(HH_s*Cs*HH_s'+Cw);


channel_mod = channel;
sys_mod = sys;
% sys_mod.n_ = sys.Nc;
% alternating optmization -------------------------
for i = 1:sys.iterAO
    % calculation of the theta - IMPLEMENTED JUST FOR DIAGONAL
    A = W_0;
    B5 = H_s;
    B4 = W_T*R5*B5;
    B3 = W_T*R4*B4;
    B2 = W_T*R3*B3;
    B1 = W_T*R2*B2;

    channel_mod.H = zeros(sys.n_rx,sys.n_tx);
    channel_mod.G = A;
    channel_mod.F = B1;
    [R1] = computeThetaDiagonal2(sys_mod,channel_mod,W_MMSE,'perfectcsi');

    A = A*R1*W_T;
    channel_mod.G = A;
    channel_mod.F = B2;
    [R2] = computeThetaDiagonal2(sys_mod,channel_mod,W_MMSE,'perfectcsi');

    A = A*R2*W_T;
    channel_mod.G = A;
    channel_mod.F = B3;
    [R3] = computeThetaDiagonal2(sys_mod,channel_mod,W_MMSE,'perfectcsi');

    A = A*R3*W_T;
    channel_mod.G = A;
    channel_mod.F = B4;
    [R4] = computeThetaDiagonal2(sys_mod,channel_mod,W_MMSE,'perfectcsi');

    A = A*R4*W_T;
    channel_mod.G = A;
    channel_mod.F = B5;
    [R5] = computeThetaDiagonal2(sys_mod,channel_mod,W_MMSE,'perfectcsi');

    % channel update
    G_SIM = A*R5;

    HH_s = G_SIM*B5;

    % calculation of the linear filter
    W_MMSE = Cs*HH_s'/(HH_s*Cs*HH_s'+Cw);
end

    HH_real = G_SIM*channel.H;
end