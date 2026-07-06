function [W_T,W_T_1,Corr_T] = generate_WSIM(n_layers,lambda,N,M)

n_max = sqrt(N);
d_layer = 5*lambda/n_layers;
dx = lambda/2;
dy = lambda/2;
d_element_SIM = lambda/2;


W_T = zeros(N, N);
Corr_T = zeros(N, N);
W_T_1 = zeros(N, M);


for pos1 = 1 : N
    pos_z = ceil(pos1/n_max);
    pos_x = mod(pos1-1,n_max)+1;
    for pos2 = 1 : N
        n_z = ceil(pos2/n_max);
        n_x = mod(pos2-1,n_max)+1;
        d_temp  = sqrt( (pos_x-n_x)^2 +  (pos_z-n_z) ^2 )*d_element_SIM;
        d_temp2 = sqrt(d_layer^2 + d_temp^2);
        cosX = d_layer/d_temp2;

        W_T(pos2,pos1) = dx*dy*(cosX/d_temp2*(1/2/pi/d_temp2-1i/lambda))*exp(1i*2*pi*d_temp2/lambda);
        Corr_T(pos2,pos1) = sinc(2*d_temp/lambda);
    end

    for pos2 = 1:M
        d_transmit = sqrt(d_layer^2 + ...
            ( (pos_x-(1+n_max)/2)*d_element_SIM - (pos2-(1+M)/2)*lambda/2)^2 + ...
            ( (pos_z-(1+n_max)/2)*d_element_SIM )^2 );
        W_T_1(pos1,pos2) = lambda^2/4*(d_layer/d_transmit/d_transmit*(1/2/pi/d_transmit-1i/lambda))*exp(1i*2*pi*d_transmit/lambda);
    end
end
end


