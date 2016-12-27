% Test 1: Look at latency and error in periodogram for different lengths N
clear
M=1; % Number simulations to do
b = 32;  % output width
sample_width=16;
dir = date;
source = 12; % uniform random number (complex)
% source = 10; % gaussian noise zero-mean
w = 0;  % rectangular window
alpha = 0.0;
% alpha = 0.1;
doplot=0;
% figure(1);
deltaf = 5/2; % ignore this, it's junk, todo remove deltaf from call to setup
delta1_sq = 2^(-2*sample_width)/12;
Nmat = [8,9,10,11,12,13];
scaling = 3;
% Px_mse = zeros(length(Nmat));
% Pb_mse = Px_mse;
% Px3_mse = Px_mse;
% Pb3_mse = Px_mse;
ratio1 = zeros(length(Nmat),1);
ratio2 = ratio1;
blkexpvect=ratio1;
lat = ratio1;
K = ratio1;
Vxdelta = ratio1;
for i = 1:M
    ax=1;
for N = Nmat
    L = N;
    setup(b, source, N, L, w, alpha, deltaf, scaling);
    % setup(b, source, N, L, w, alpha);
    load(parfile);
    x = x + rand(size(x)) + 1j*rand(size(x));
    % matfile = ['s', int2str(source), 'w', int2str(w)...
    %     'N', int2str(N), 'L', int2str(L), 'b', int2str(b),...
    %     'alpha', int2str(alpha), 'scale', int2str(scaling)];
    % save([dir, '/', matfile, '_no2']);
    if scaling == 3
        sim('bart_block_floating_point')    
    else
        sim('bart_block');     
    end
    
    % fixpoint_analysis;
    get_sim_results
    sfactor = 1/2^(N-blkexp);
    blkexpvect(ax) = blkexp;
    Xf = Xf(:,1,2);
    Xf = Xf*sfactor;
    Px = Px(:,1,2);                 
    Xf2 = fft(xwin(:,1,2),Nfft);
    Xf2 = Xf2/Nfft;
    % Xf2 = Xf2(:,1,2);
    Px2 = Xf2.*conj(Xf2);
    xwin = xwin(:,1,2);
    K(ax) = 1/Nfft*sum(xwin.*conj(xwin));
    err1 =  sqrt(sum((Xf-Xf2).*conj(Xf-Xf2)))./sqrt(K(ax));
    % K(ax) = sum(xwin.*conj(xwin))/Nfft;
    % err1 =  sqrt(sum((Xf-Xf2).*conj(Xf-Xf2)))./sqrt(K(ax)*Nfft);
    % err2 = sum(abs(Px-Px2))./sqrt(K(ax));
    err2 = sum(abs(Px-Px2))/var(xwin);
    ratio1(ax) = ratio1(ax) + err1;
    ratio2(ax) = ratio2(ax) + err2;
    lat(ax) = px_latency-2^N;
    ax = ax + 1;
end
end
matfile = ['acc_','s', int2str(source), 'w', int2str(w)...
    'b', int2str(b),...
    'scale', int2str(scaling)];
delta1_sq = 2^(-2*sample_width)/12;
lubound = sqrt((Nmat-2.5)*delta1_sq);
% My upper bound seems to be off by factor of two
uubound = sqrt(2.^(Nmat+3)*delta1_sq)./sqrt(K'./2.^Nmat)*2;  
save([dir, '/', matfile],'ratio1','ratio2','Nmat','K','lubound','uubound','lat','M','scaling','blkexpvect');

%% Plot

figure(1)
ratio1 = ratio1/M;
scatter(Nmat,log2(ratio1),'*');
hold('on')
if scaling == 3
    plot(Nmat,log2(uubound));
    legend('error', 'upper-upper bound')
else
    plot(Nmat,(log2(lubound)));
    legend('error', 'lower-upper bound');
end
clear title;
title('Rms(err)/rms(result) FFT ');
xlabel('log2(Nfft)')
set(gca, 'XTick', Nmat)
ylabel('rms(err)/rms(result)')
print([dir, '/test1_log2mse.png'], '-dpng');
figure(2)
scatter(Nmat,ratio2,'*');
title('Rms(err)/rms(result) PSD ');
figure(3)
scatter(Nmat,lat./(2.^Nmat)');
title('Ratio of latency and FFT points');
xlabel('log2(Nf ft)');
set(gca, 'XTick', Nmat)
ylabel('cycles/Nfft');
print([dir, '/test1_lat.png'], '-dpng');
