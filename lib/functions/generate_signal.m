%OBS: Normalization (the noise/signal by R is normalized at the main code)

function [TxSymbolData, Packet] = generate_signal(setup)
    % Variables initialization
    sigma_s = 1;
    G_LDPC = setup.idd3.G;
    N = setup.idd3.N;
    PacketDataLength = setup.sys.PacketDataLength;
    n_tx = setup.sys.n_tx;
    M_ORD = setup.sys.M_ORD;
    
    % Pilots
    num_pilots = setup.estimation.pilots_coarse;
    num_bit_pilots = num_pilots*M_ORD;
    
    if num_pilots
        pilots = hadamard(num_pilots);
        pilots = pilots(1:n_tx,:).*exp(1i*pi/4)*1.34;
        data_pilots  = qamdemod(pilots.',2^M_ORD,'OutputType','bit','UnitAveragePower',true).';
    else
        data_pilots = [];
    end
    
    % LDPC coding and modulation ---------------------------------------------------------------
    PacketBitData = randi([0 1], n_tx ,PacketDataLength-num_bit_pilots);
    Packet = [data_pilots PacketBitData];
    TxDatabits = zeros(n_tx,N);
    for i = 1:n_tx
        TxDatabits(i,:)= mod(Packet(i,:)*G_LDPC,2);
    end
    TxSymbolData = qammod(TxDatabits',2^M_ORD,'InputType','bit','UnitAveragePower',true).'*(sigma_s);

end