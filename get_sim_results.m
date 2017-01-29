% Get results from simulation
px_latency = min(find(px_valid.data));
fft_latency = min(find(Xk_valid.data));
latency = fft_latency;
% latency = min(find(bart_valid.data));
% Px_b = px_bart.data(latency:latency+Nfft*ncols-1);
% Px_b = reshape(Px_b, [Nfft, 1, ncols])/nsect;
% Px_b = reshape(Px_b, [Nfft, 1, ncols]);
xq = xq_re.data + 1j*xq_im.data;
start=min(find(xq_valid.data));
xq=xq(start+1:start+xlength);
xq = reshape(xq,[Nfft,nsect,ncols]);
start=min(find(xwin_valid.data));
xwin = xwin_re.data + 1j*xwin_im.data;
% xwin = xwin(start+1:start+xlength);
xwin = xwin(find(xwin_valid.data));
xwin = xwin(2:end);  
xwin = reshape(xwin,[Nfft,nsect,ncols]);
% Px = px.data(px_latency:px_latency+xlength-1);
Px = px.data(find(px_valid.data));
% Px = reshape(Px, [Nfft, nsect, ncols]);
Px = reshape(Px, [Nfft, 1, ncols]);
% Xf = Xk_re.data(fft_latency:fft_latency+xlength-1) ...
%   + 1j*Xk_im.data(fft_latency:fft_latency+xlength-1);
Xf = Xk_re.data(find(Xk_valid.data)) + 1j*Xk_im.data(find(Xk_valid.data));
% if scaling == 3
blkexp = blk_exp.data(min(find(blk_exp_status.data)));
Px = Px/2^(2*(N-blkexp))*1/2^(L-N);
    % Xf = Xf.*2^blkexp;
% else
%     blkexp = 0;
% end
if ~exist('ncols','var')
    ncols = 1;
end
% Xf = reshape(Xf, [Nfft, nsect, ncols]);