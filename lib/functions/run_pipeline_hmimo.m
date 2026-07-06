function out = run_pipeline_hmimo(setup, ch, TxSymbolData, Packet)
    out = struct();
    % --- stage 1 ---
    [ch_est,out.coarse.nmse] = ce_lmmse(setup, ch, 'initial', TxSymbolData, 0);
    [out.coarse.bErrors, out.coarse.rsum, RecoveredPacket] = idd_mmse_hmimo(setup, ch_est, 'estimated',TxSymbolData, Packet);
    % --- stage 2 ---
    [ch_est,out.refined.nmse] = ce_lmmse(setup, ch, 'improved', TxSymbolData, RecoveredPacket);
    [out.refined.bErrors, out.refined.rsum, ~] = idd_mmse_hmimo(setup, ch_est, 'estimated',TxSymbolData, Packet);
    % --- stage 3 ---
    [~,out.perfect.nmse] = ce_lmmse(setup, ch, 'benchmark', TxSymbolData, 0);
end