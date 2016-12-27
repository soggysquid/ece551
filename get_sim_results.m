% Get results from simulation
px_latency = min(find(px_valid.data));
fft_latency = min(find(Xk_valid.data));
start=min(find(xwin_valid.data));
xwin = xwin_re.data + 1j*xwin_im.data;
xwin = xwin(start+1:start+xlength);
xwin = reshape(xwin,[Nfft,nsect,ncols]);
Px = px.data(px_latency:px_latency+xlength-1);
Px = reshape(Px, [Nfft, nsect, ncols]);
Xf = Xk_re.data(fft_latency:fft_latency+xlength-1) ...
   + 1j*Xk_im.data(fft_latency:fft_latency+xlength-1);
if scaling == 3
    blkexp = blk_exp.data(max(find(blk_exp.data)));
    Px = Px/2^(2*(N-blkexp));
    % Xf = Xf.*2^blkexp;
else
    blkexp = 0;
end
if ~exist('ncols','var')
    ncols = 1;
end
Xf = reshape(Xf, [Nfft, nsect, ncols]);