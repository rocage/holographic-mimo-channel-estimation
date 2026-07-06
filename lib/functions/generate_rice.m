
function [ch] = generate_rice(setup)

        % A. Init 'sys' parameters 
        n_tx = setup.sys.n_tx;
        Nc = setup.sys.Nc;
        pos_ap = setup.enviroment.pos_ap;
        lambda = setup.sys.lambda;

        % B. Init 'enviroment' parameters
        pos_center_users = setup.enviroment.pos_center_users;
        radius_users = setup.enviroment.radius_users;
        kappa = setup.enviroment.kappa;
        Corr_T = setup.enviroment.Corr_T;
        HPL = setup.enviroment.PLsim;

        % C. Geral parameters
        d = lambda/2;
        k_wave = 2*pi/lambda;
        Nx = sqrt(Nc);
        Ny = sqrt(Nc);
        user_position = zeros(n_tx,3);

        
        % D. Sort positions of users and calculate the distances 
        for i = 1:n_tx
            user_position(i,:) = pos_center_users' +[radius_users*cos(2*pi*rand()) radius_users*sin(2*pi*rand()) 0];
        end

        [dist_ap_k,angle_ap_k] = rangeangle(pos_ap,user_position');
        dist_ap_k = dist_ap_k.';
        phi_users = deg2rad(angle_ap_k(1,:));

        
        % E. Sort nlos channel coefficients using rayleight Distribution 
        H_nlos = ( randn(n_tx,Nc)+1j*randn(n_tx,Nc))/sqrt(2);
        H_nlos = (Corr_T)^(1/2)*H_nlos';

        % E. Sort los channel coefficients using users position 
        for k = 1:n_tx
            phi = phi_users(k);

            % E.1 fase depende apenas do seno do azimute
            iy = 0:Ny-1;
            psi_y = k_wave * d * sin(phi);
            a_y = exp(-1j * iy * psi_y).'; % Vetor Coluna (Ny x 1)

            % E.2 Vetor Vertical (X)
            % Como theta = 0, não há defasagem vertical relativa entre linhas
            a_x = ones(Nx, 1);

            % E.3 Produto de Kronecker
            H_los(:, k) = kron(a_x, a_y);
        end
   
        % F. compute rice coefficients
        H_rice = sqrt(kappa/(kappa+1))*H_los + sqrt(1/(kappa+1))*H_nlos;
        H_rice_loss = ((H_rice.').*sqrt(HPL(dist_ap_k))).';

        % G. Create the output structure
        ch.H = H_rice_loss;
        ch.Hkloss=HPL(dist_ap_k);
        ch.Corr_T = Corr_T;
        
end
