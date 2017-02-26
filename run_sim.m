% sample-width: fft sample width (adc width is fixed)
% width_out: width of output
% source: data set to use
% source:
%   1: iq .dat
%   2: Stephen's csv 
%   3: slot test
%   4: one tone in bin
%   5: two tones in bin
%   6: one tone between bins
%   7: two tones out of bin
%   8: impulse
%   9: all ones*A1
%   10: all zeros
%   11: a bunch of tones at harmonics of the sampling frequency
%   12: complex uniform random numbers between 0 and 1
% N: log2(length of fft)
% L: log2(length of data)
% w: window to use
%   0: rectangular
%   1: kaiser
%   2: hanning
%   3: blackman
% alpha: noise variance
% hwver:
%   0: block foating point
%   1: no scaling
%   2: block floating point w spectral smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set width_out to 0 to use full precision
width_out = 0;
source = 12;
A2 = 2^-13;
A2 = 0;
alpha = 0;
A = 1 - A2 - 2^-16 - 3*alpha;
% A = 1-2^-16;
N = 10;  % Length of FFT, must be < L
Nfft = 2^N;
L = 10;  % Length of sample
w = 0;
% alpha = 2^-15;
hwver = 1;
dir = 0;     % If you want to simulate
% If you're loading files...
runsim = 1;
if 1     % load instead of running results
    runsim = 0
    dir = '25-Feb-2017/errAnalysis_s12';
    i=1          
    avg=L-N;
    matfile = ['s', int2str(source), 'w', int2str(w)...
               'N', int2str(N), 'L', int2str(L), 'b',...
               int2str(width_out),'a', num2str(alpha),...
               'avg', int2str(avg), 'hw', int2str(hwver),...
               'iter', int2str(i), '.mat']
    load([dir,'/',matfile]);
%%
else
    rng(133);
    if hwver == 2
        w=2;
    end
    Error=0;
    setup(source,width_out,N,L,w,hwver,A,alpha,A2)
    load(vectfile);
    if ~Error
        if hwver == 0
            sim('bart_block_scaled.slx')
        elseif hwver == 1
            sim('bart_block.slx')
        else
            sim('bart_block_spectsmoothed_scaled.slx')
        end
    end
end
fixpoint_analysis