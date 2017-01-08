clear;
source = 6; 
N=5;
L=N;
Nfft = 2^N;
sigma = 0;
xmin=0;
deltaf=0;
dir='.';
[x,findex,fs] = create_test_vector(source, Nfft, L, sigma, xmin, deltaf, dir);
k=[-Nfft/2:Nfft/2];
k=k(1:end-1)';
wb=[ones(Nfft,1)].*(Nfft-abs(k))/Nfft;
wr=[ones(Nfft,1)];
wh=hann(Nfft+1);
wh=wh(1:end-1);
xh=wh.*x;
figure(1)
[Xper,Wb] = func_make_bias_plot(x,wr,k,Nfft,1,1);
figure(2)
Ewin = sum(wh.^2)/Nfft;
[Xperh,Wh] = func_make_bias_plot(xh,wh,k,Nfft,1,Ewin);
