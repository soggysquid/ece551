% Get results from simulation
if runsim
    px_latency = min(find(px_valid.data));
    xq = xq_re.data + 1j*xq_im.data;
    xq = xq(find(xq_valid.data));
    xq = reshape(xq,[Nfft,nsect,ncols]);
    start=min(find(xwin_valid.data));
    Px = px.data(find(px_valid.data));
    Xf = Xk_re.data(find(Xk_valid.data)) + 1j*Xk_im.data(find(Xk_valid.data));
    if width_out < width_perio - bp_perio
        Px = Px*2^(width_perio-bp_perio-width_out);
    end
    blkexp = blk_exp.data(min(find(blk_exp_status.data)));
    Qfft = fft(x(:,1:nsect,1:ncols)-xq)/2^Nmax;
    Qpsd = Qfft.*conj(Qfft);
    Qpsd = circshift(Qpsd, Nfft/2-1);
    if hwver==1
        Px_u = Px/2^(2*(Nmax-N));
        % Px = Px*2^(2*(-L));
        Px = Px/(2^(2*Nmax));
    elseif hwver == 2
       Px_u = 2^(2*blkexp)*Px; 
       Px = 2^(2*blkexp)*Px/2^(2*N);  
    else
        % Px_u = Px*2^(2*blkexp-L-(L-N));
        % Px_u = 2^(2*blkexp)*Px*2^(N-L);
        Px_u = 2^(2*blkexp)*Px;
        % Px = 2^(2*blkexp)*Px*2^(-L-N);  
        % Px = 2^(2*blkexp)*Px*2^(-L-N)*2^(L-N); 
        Px = 2^(2*blkexp)*Px/2^(2*N);  
    end

end
% Px = reshape(Px, [2^N, 1, ceil(length(Px)/Nfft)]); % changed this feb4
% Px = reshape(Px, [2^N, 1]); 
% Px = Px(:,1);                                                                                                               

if ~exist('ncols','var')
    ncols = 1;
end
