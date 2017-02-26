clear;
source = 6; 
N=5;
L=N;
Nfft = 2^N;
sigma = 0;
xmin=0;
deltaf=0;
dir='.';
A=1;
[x,findex,fs] = create_test_vector(source, Nfft, L, A, sigma, xmin, deltaf, dir);
source = 4;
[x1,findex,fs] = create_test_vector(source, Nfft, L, A, sigma, xmin, deltaf, dir);
k=[-Nfft/2:Nfft/2];
k=k(1:end-1)';
wb=[ones(Nfft,1)].*(Nfft-abs(k))/Nfft;
wr=[ones(Nfft,1)];
wh=hann(Nfft+1);
wh=wh(1:end-1);
xh=wh.*x;
xh1=wh.*x1;
figure(1)
Ewin = 1
[Xper,Wb] = func_make_bias_plot(x,wr,k,Nfft,1,Ewin);
figure(2)
[Xper,Wb] = func_make_bias_plot(x1,wr,k,Nfft,0,Ewin);
figure(3)
Ewin = sum(wh.^2)/Nfft;
[Xperh,Wh] = func_make_bias_plot(xh,wh,k,Nfft,1,Ewin);
figure(4)
[Xperh,Wh] = func_make_bias_plot(xh1,wh,k,Nfft,0,Ewin);

