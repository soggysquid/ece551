function [ Xper, W ] = func_make_bias_plot( x, w, k, Nfft,doplot,Ewin)
w=[zeros(Nfft/2*(Nfft-1),1);w;zeros(Nfft/2*(Nfft-1)-1,1)]; %zero-pad
Xper = zeros(Nfft^2,1);
Xper(1:Nfft:length(Xper)) = (fft(x).*conj(fft(x)))/Nfft^2;
Xper = Xper/Ewin;
% w_fft = 0.5*fft(w,Nfft*Nfft)/Nfft;
w_fft = fft(w,Nfft*Nfft)/Nfft;
w_fft = w_fft.*w_fft*0.5/Ewin;
w_fft1 = circshift(w_fft, Nfft*Nfft/2);
w_fft2 = circshift(w_fft, Nfft*(Nfft/2+1));
W=w_fft1+w_fft2;
if doplot
    scatter([0:length(Xper)-1]/Nfft,10*log10(Xper));
    hold('on')
    plot([0:length(Xper)-1]/Nfft,10*log10(abs(W)))
    xlim([0 Nfft]);
    ylim([-100 5]);
    ylabel('dB')
    xlabel('bin')
end


