clear; clc;
rng(100);

% Parâmetros mínimos
Nc = 64;  % ajusta se for outro valor no seu setup
R = 3;

W_T = randn(Nc) + 1i*randn(Nc);
W_T_1 = randn(Nc) + 1i*randn(Nc);

% --- VERSÃO ORIGINAL (copiada literalmente do CE_BD_SIM_RIS_LMMSE_EP_3layer) ---
rng(42);  % seed fixa
[Theta1_old,~] = qr(randn(Nc) + 1i*randn(Nc));
[Theta2_old,~] = qr(randn(Nc) + 1i*randn(Nc));
[Theta3_old,~] = qr(randn(Nc) + 1i*randn(Nc));
aux_old = W_T_1*Theta1_old*W_T*Theta2_old*W_T*Theta3_old;

% --- VERSÃO UNIFICADA (build_cascade) ---
rng(42);  % MESMA seed
aux_new = build_cascade("BD", R, Nc, W_T, W_T_1);

% Comparar
fprintf('Diff aux: %.3e\n', norm(aux_old - aux_new, 'fro'));
fprintf('Diff norm: %.3e\n', abs(norm(aux_old,'fro') - norm(aux_new,'fro')));

function aux = build_cascade(arch, n_layers, Nc, W_T, W_T_1)
    aux = W_T_1 * sample_layer(arch, Nc);
    for r = 2:n_layers
        aux = aux * W_T * sample_layer(arch, Nc);
    end
end

function Theta = sample_layer(arch, Nc)
    switch arch
        case "SIM"
            Theta = diag(exp(1j * 2*pi * rand(Nc,1)));
        case "BD"
            [Theta, ~] = qr(randn(Nc) + 1i*randn(Nc));
    end
end