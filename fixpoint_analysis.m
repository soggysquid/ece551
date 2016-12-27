% Analysis of hardware implemented PSD
% load('parameters');
if ~exist('doplot')
    doplot = 1;
end
%% Get FFT output
fft_latency = min(find(Xk_valid.data));
Xf = Xk_re.data(fft_latency:fft_latency+xlength-1) ...
   + 1j*Xk_im.data(fft_latency:fft_latency+xlength-1);
if ~exist('ncols')
    ncols = 1;
end
Xf = reshape(Xf, [Nfft, nsect, ncols]);
% Get periodgram output
perio_latency = min(find(px_valid.data));
Px = px.data(perio_latency:perio_latency+xlength-1);
Px = reshape(Px, [Nfft, nsect, ncols]);
% Get Bartlett output
latency = min(find(bart_valid.data));
Px_b = px_bart.data(latency:latency+Nfft*ncols-1);
% Px_b = reshape(Px_b, [Nfft, 1, ncols])/nsect;
Px_b = reshape(Px_b, [Nfft, 1, ncols]);
xq = xq_re.data + 1j*xq_im.data;
start=min(find(xq_valid.data));
xq=xq(start+1:start+xlength);
start=min(find(xwin_valid.data));
xwin = xwin_re.data + 1j*xwin_im.data;
xwin = xwin(start+1:start+xlength);
% xwin = xwin_re.data(find(xwin_valid.data)) + 1j*xwin_im.data(find(xwin_valid.data));
% xwin = xwin(1:Nfft*nsect*ncols);
% xq = xwin;
xq = reshape(xq,[Nfft,nsect,ncols]);
xwin = reshape(xwin,[Nfft,nsect,ncols]);
if scaling == 3
    blkexp = blk_exp.data(max(find(blk_exp.data)));
end

%% Analysis
% sf = 2^(N-1);
sf = 2^Nmax;
% sf=1;
sf2 = 2^N;
% sf = 1;
Qfft = fft(x(:,1:nsect,1:ncols)-xq)/sf;
Qpsd = Qfft.*conj(Qfft);
% xwin = xq .* repmat(win, [1 nsect ncols]);
Xf2 = fft(xwin, Nfft);
Xf3 = fft(x, Nfft);
if source ~= 3
    Xf2 = Xf2/sf2;
    Xf3 = Xf3/sf2;
    % Xf = Xf/sf;
    % Px = Px/sf^2;
    % Px_b = Px_b/sf^2*1/2^bartmax;
else
    Xf = Xf*sf;
    Px = Px*sf^2;
    Px_b = Px_b*sf^2;
    % Px_b = Px_b/2^bartmax;
end
if scaling == 3
    Xf = Xf*2^blkexp;
    Xf = Xf/sf2;
    Px = Px*2^(blkexp*2);
    Px = Px/sf2^2;
end
% Px_b = Px_b/nsect;
Px2 = Xf2 .* conj(Xf2) * Ewin;
Px3 = Xf3 .* conj(Xf3) * Ewin;
Px2_b = sum(Px2,2)/nsect;
Px3_b = sum(Px3,2)/nsect;
Qpsd = circshift(Qpsd, Nfft/2-1);
Xf = circshift(Xf, Nfft/2-1);
Xf2 = circshift(Xf2, Nfft/2-1);
Px = circshift(Px, Nfft/2-1);
Px_b = circshift(Px_b, Nfft/2-1);
Px2 = circshift(Px2, Nfft/2-1);
Px2_b = circshift(Px2_b, Nfft/2-1);
Px3 = circshift(Px3, Nfft/2-1);
Px3_b = circshift(Px3_b, Nfft/2-1);
findex = [-2^(N-1)+1:2^(N-1)]'/(Ts*2^N);
if doplot
   figure(1)
   subplot(2,1,1)
   plot_em(Px2, findex, 'Periodogram', fc, width_out, Px);
   subplot(2,1,2);
   plot_em(Px2_b, findex, 'Bartlett', fc, width_out, Px_b);
   figure(2)
   plot_em(Qpsd, findex, 'Quantization Noise PSD', fc, 64);
   % figure(3)
   % subplot(2,1,1)
   % plot_em(Px3, findex, 'No ADC Periodogram', fc);
   % subplot(2,1,2)
   % plot_em(Px3_b, findex, 'No ADC Bartlett', fc);
end
if source == 3
    pk_index = find(Px>min(mean(Px)));
    pk2_index = find(Px2>min(mean(Px2)));
    n_index = find(Px<min(mean(Px)));
    n2_index = find(Px2<min(mean(Px)));    
    Px_noise = Px(n_index);
    Px2_noise = Px2(n2_index);
    bias_fix = 10*log10(mean(Px_noise))
    bias_flt = 10*log10(mean(Px2_noise))
    var_fix = 10*log10(var(Px_noise))
    var_flt = 10*log10(var(Px2_noise))
end