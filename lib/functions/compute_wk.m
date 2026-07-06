function [Wk] = compute_wk(HH,n_tx,n_rx,sigma2_vec,Ex,sigma_s,sigma_n)

                Wk = zeros(n_tx,n_rx);

                for k = 1:n_tx
                    % Delta_k - covariance matrix
                    aux = sigma2_vec;
                    aux(k) = 1;
                    Delta_x = diag(aux/Ex);

                    % matrix evaluation (LRA-MMSE)
                    G_MMSE = sigma_n^2/sigma_s^2*eye(size(HH,1)) + HH*Delta_x*HH';
                    
                    % slower implementation --------------------------
                    %wk = pinv(G_MMSE)*HH(:,k);
                    % more fast implematation https://www.mathworks.com/help/matlab/ref/pinv.html
                    wk = lsqminnorm(G_MMSE,HH(:,k));
                    
                    Wk(k,:) = wk';
                end