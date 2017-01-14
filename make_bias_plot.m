clear;
source = 6; 
N=5;
L=N;
Nfft = 2^N;
sigma = 0;
xmin=0;
deltaf=0;
inv_delta=Nfft;
dir='.';
[x,findex,fs] = create_test_vector(source, Nfft, L, sigma, xmin, deltaf, dir);
% f1 = 1/8 + 0.5/Nfft;
% x = A1*(cos(2*pi*n*f1)+1j*sin(2*pi*n*f1));
Xper = zeros(Nfft*inv_delta,1);
Xper(1:Nfft:length(Xper)) = (fft(x).*conj(fft(x)))/Nfft^2;
figure(1);
% scatter([1:length(Xper)],10*log10(Xper));
scatter([0:length(Xper)-1]/Nfft,10*log10(Xper));
% scatter([0:length(Xper)-1]/Nfft,Xper);
%%
k=[-inv_delta/2:inv_delta/2];
n=[-Nfft/2:1/inv_delta:Nfft/2-1/inv_delta];
k=k(1:end)';
wb=[ones(inv_delta+1,1)].*(inv_delta-abs(k))/inv_delta;
wb=[zeros(inv_delta/2*(Nfft-1),1);wb;zeros(inv_delta/2*(Nfft-1)-1,1)];
wb_fft = 0.5*fft(wb,Nfft*inv_delta)/Nfft;
wb_fft1 = circshift(wb_fft, Nfft*inv_delta/2);
wb_fft2 = circshift(wb_fft, Nfft*(inv_delta/2+1));

Wb=wb_fft1+wb_fft2;
hold('on')
plot([0:length(Xper)-1]/Nfft,10*log10(abs(Wb)))
xlim([0 Nfft]);



 