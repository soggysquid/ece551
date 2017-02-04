function [ Xper, W ] = func_make_bias_plot( x, w, k, Nfft, type, Ewin)
doplot = 1;
w=[zeros(Nfft/2*(Nfft-1),1);w;zeros(Nfft/2*(Nfft-1)-1,1)]; %zero-pad
Xper = zeros(Nfft^2,1);
Xper(1:Nfft:length(Xper)) = (fft(x).*conj(fft(x)))/Nfft^2;
Xper = Xper/Ewin;
Xmax = 10*log10(max(Xper))
% w_fft = 0.5*fft(w,Nfft*Nfft)/Nfft;
w_fft = fft(w,Nfft*Nfft)/Nfft;
if type
    w_fft = w_fft.*w_fft*0.5/Ewin;
    w_fft1 = circshift(w_fft, Nfft*Nfft/2);
    w_fft2 = circshift(w_fft, Nfft*(Nfft/2+1));
    W=w_fft1+w_fft2;
else
    w_fft = w_fft.*w_fft/Ewin;
    W = circshift(w_fft, Nfft*Nfft/2);
end

unity = ones(length(Xper),1);
if doplot
    scatter([0:length(Xper)-1]/Nfft,10*log10(Xper));
    hold('on')
    plot([0:length(Xper)-1]/Nfft,10*log10(abs(W)))
    plot([0:length(Xper)-1]/Nfft,Xmax*unity,'--')
    x1 = 3;
    y1 = Xmax+2;
    txt = sprintf('%0.2f dB',Xmax);
    text(x1,y1,txt)
    xlim([0 Nfft]);
    ylim([-100 5]);
    ylabel('dB')
    xlabel('bin')
end


