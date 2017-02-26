% Get results from simulation
if runsim
    px_latency = min(find(px_valid.data));
    % fft_latency = min(find(Xk_valid.data));
    % latency = fft_latency;
    % latency = min(find(bart_valid.data));
    % Px_b = px_bart.data(latency:latency+Nfft*ncols-1);
    % Px_b = reshape(Px_b, [Nfft, 1, ncols])/nsect;
    % Px_b = reshape(Px_b, [Nfft, 1, ncols]);
    xq = xq_re.data + 1j*xq_im.data;
    % start=min(find(xq_valid.data));
    % xq=xq(start+1:start+xlength);
    xq = xq(find(xq_valid.data));
    xq = reshape(xq,[Nfft,nsect,ncols]);
    start=min(find(xwin_valid.data));
    % xwin = xwin_re.data + 1j*xwin_im.data;
    % xwin = xwin(start+1:start+xlength);
    % xwin = xwin(find(xwin_valid.data));
    % xwin = xwin(2:end-1);
    % xwin = xwin(1:end-1);  
    % xwin = reshape(xwin,[Nfft,nsect,ncols]);
    % Px = px.data(px_latency:px_latency+xlength-1);
    Px = px.data(find(px_valid.data));
    % Px = reshape(Px, [Nfft, nsect, ncols]);
    Xf = Xk_re.data(find(Xk_valid.data)) + 1j*Xk_im.data(find(Xk_valid.data));
    if width_out < width_perio - bp_perio
        Px = Px*2^(width_perio-bp_perio-width_out);
    end
 % if scaling == 3
    blkexp = blk_exp.data(min(find(blk_exp_status.data)));
    Qfft = fft(x(:,1:nsect,1:ncols)-xq)/2^Nmax;
    Qpsd = Qfft.*conj(Qfft);
    % xwin = xq .* repmat(win, [1 nsect ncols]);
    Qpsd = circshift(Qpsd, Nfft/2-1);
    if hwver==1
        Px_u = Px/2^(2*(13-N));
        Px = Px*2^(2*(-L));
    else
        Px_u = Px*2^(2*blkexp-(L-N));
        Px = Px*2^(2*(blkexp-L));
    end

end
Px = reshape(Px, [2^N, 1, ceil(length(Px)/Nfft)]); % changed this feb4
Px = Px(:,1);                                                                                                               
% Px = Px
% Px = 2^bartmax*Px;
% Px = Px/2^(Nmax-N);
% Xf = Xk_re.data(fft_latency:fft_latency+xlength-1) ...
%   + 1j*Xk_im.data(fft_latency:fft_latency+xlength-1);
% Px = Px/2^(2*(N-blkexp))*1/2^(L-N);


if ~exist('ncols','var')
    ncols = 1;
end
% Xf = reshape(Xf, [Nfft, nsect, ncols]);