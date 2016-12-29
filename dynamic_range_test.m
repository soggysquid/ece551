
% Test 1: Look at latency and error in periodogram for different lengths N
clear
b = 32;  % output width
sample_width=16;
dir = date;
source = 5; % uniform random number (complex)
% source = 10; % gaussian noise zero-mean
w = 0;  % rectangular window
alpha = 0.0;
% alpha = 0.1;
doplot=0;
% figure(1);
deltaf = 5/2; % ignore this, it's junk, todo remove deltaf from call to setup
Nmat = [8,9,10,11,12,13];
scaling = 3;
xminvect=zeros(length(Nmat),1);
i=1;
for N=Nmat
    cont = 1;
    expn = sample_width;
    expn = 15;
    Nfft = 2^N
    while cont & expn > 1
        L = N;
        setup(0, b, N, L, w, scaling);
        xmin = 2^-expn;
        create_test_vector(source,Nfft,L,alpha,xmin,deltaf,dir);
        if scaling == 3
            sim('bart_block_floating_point')  
        else
            sim('bart_block')  
        end
        get_sim_results;
        if length(find(Px))>1
            cont=0;
        else
            expn = expn-1
        end
    end
    xminvect(i)=expn;
    i=i+1;
end
            
