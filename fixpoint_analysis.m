% Analysis of hardware implemented PSD
% load('parameters');
if ~exist('doplot')
    doplot = 1;
end
% Note to self: to merge two figures, fig1 and fig2
% L = findobj(1,'type','line');
% copyobj(L,findobj(2,'type','axes'));
%% Get FFT output
get_sim_results

%% Analysis

xwin2 = xq.*win;
Xf2 = fft(xwin2, Nfft);
Px2_u = (Xf2 .* conj(Xf2)) * Ewin;
Px2 = (Xf2/Nfft .* conj(Xf2/Nfft)) * Ewin;
nsect = 2^(L-N);
Px2_b = sum(Px2,2)/nsect;
Px = circshift(Px, Nfft/2-1);
Px2 = circshift(Px2, Nfft/2-1);
Px2_b = circshift(Px2_b, Nfft/2-1);
% findex = [-2^(N-1)+1:2^(N-1)]'/(Ts*2^N);
if doplot
   figure(1)
   % subplot(2,1,1)
   plot_em(Px2_b, 'Periodogram', width_out, Px);
   % subplot(2,1,2);
   % plot_em(Px2_b, findex, 'Bartlett', fc, width_out, Px_b);
   % figure(2)
   % plot_em(Qpsd, findex, 'Quantization Noise PSD', fc, 64);
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