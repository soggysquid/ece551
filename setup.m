% Script to setup all the parameters and teuntitledst vectors that we're going to
% use for modelling and analysis
% clearvars -except latency sample_width twiddle_width;
function setup(source, width_out, N, L, w, hwver, A, alpha, A2)
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

switch nargin
    case 0
        % set width_out to 0 to use full precision
        width_out = 0;
        source = 4;
        Delta = 2/2^16;
        % alpha = sqrt(2^10*Delta^2); % This is a great value, don't erase
        alpha = 0.0167;
        A = 0;
        A2 = 1e-4;
        alpha = 2^-9;
        N = 13;  % Length of FFT, must be < L
        L = 13;  % Length of sample
        w = 2;
        % alpha = 2^-15;
        % alpha = 0.0;
        hwver = 1;
end
Nmax = 13;
deltaf = 1600/2^Nmax;      % Make deltaf 60 bins away
dosim=0;
if N < 8
    print 'Error, N < 8'
end
if ~exist('xmin')
    xmin = 2^-8;
    xmin=0;
end
if ~exist('alpha')
    alpha = 0;
end
if ~exist('A2')
    A2 = 1e-4;
end
sources = {' iq .dat ', ' usrp ', ...
    ' slot ', ' 1 tone in bin ', ' 2 tones in bins ', ' 1 inter-bin tone ',...
    ' 2 inter-bin tones ', 'impulse', 'all ones', 'all zeros'}; 
wins = {'rect', 'kaiser', 'hanning', 'blackman'};
% Nmax = 13; % log2 of max points of run-time configurable FFT
bartmax = Nmax-8;  % max number segments we can average 
float = 0;
adc_width = 16;
sample_width = adc_width;
dir = date;
if ~exist(dir, 'dir')
    mkdir(dir)
end
parfile = [dir, '/parameters'];  % name of .mat file to save parameters to
assignin('base', 'parfile', parfile)
imfile = 'data_imaginary4_11052015.xls';
refile = 'data_real4_11052015.xls';
datfile = 'LB3bp70f241sf64.dat';
% quantization options for simulating adc
rmode = 'trunc';
lmode = 'sat'; 
Nfft = 2^N;
latency = 2^(L+1) + 200;
simtime = latency + 2^L;
% alpha = 0.0;
usemex = 0;
% bp_reinterp = bp_fft+Nmax;
% width_perio = width_fft*2-1; % full precision  
% factors to put in ROM for Welch algorithm
% PWelch size FFT size is half Nmax because we have two of them
% nt = [0:Nmax/2-1]
% twiddles = exp(-j*pi*n/Nmax);
% realtwiddle = real(twiddles);
% imagtwiddle = imag(twiddles);

if hwver == 1  % no scaling
    width_fft = sample_width+Nmax;  % this is the unsigned width
    bp_fft = sample_width-1;
    % bp_fft = bp_fft+Nmax;
    bp_perio = bp_fft*2-3;  
    % width_perio = width_fft*2+1;
    width_perio = width_fft*2;   % addition doesn't add a bit because if 
                                 % one component is 2 the other must be 0
    % width_fpadd = width_perio+bartmax;
    width_fpadd = width_perio+bartmax;
    bp_reinterp = bp_fft;
% elseif hwver == 2
%     width_fft = sample_width+Nmax+1;
%     bp_fft = width_fft-2;
%     width_perio = 47; % Filter is restricted to 47 bit width
%     bp_perio = bp_fft*2-(width_fft*2-47);
%     width_fpadd = width_perio+bartmax;
%     % bp_reinterp = bp_fft;
else
    width_fft = sample_width;
    % bp_fft = width_fft-2;
    bp_fft = width_fft-1;
    bp_perio = bp_fft*2-3;  % Shift because of window energy
    width_perio = (width_fft-1)*2+1;
    width_fpadd = width_perio+bartmax;
end
bp_fpadd = bp_perio;
if ~width_out
    % width_out=width_fpadd;
    % width_out = width_fpadd;
    width_out = width_perio;
    % bp_out = bp_fpadd;
end
bp_out = max(width_out-(width_perio-bp_perio),0);
% bp_out = max(width_out-(width_perio-bp_perio),0);
% bp_out = bp_fpadd+bartmax;
% bp_bart = bp_fpadd+bartmax;   % for full precision bart output
% if width_out
%     % bp_out = bp_bart-(width_fpadd-width_out);
%     bp_out = bp_fpadd-(width_fpadd-width_out)+bartmax;
%     % bp_out = bp_fpadd-(width_fpadd-width_out);
% else
%     width_out = width_fpadd;
%     bp_out = bp_bart;
% end
%% Set up windowing
if hwver == 2
    w = 2
end
if w == 1
    win = kaiser(Nfft+1,9);
    Ewin = 1/(sum(kaiser(Nfft+1,9).^2)/Nfft);
elseif w == 2 
    win = hanning(Nfft+1);
    Ewin = 1/(sum(hanning(Nfft+1).^2)/Nfft);
elseif w == 3
    win = blackman(Nfft+1);
    Ewin = 1/(sum(blackman(Nfft+1).^2)/Nfft);
else
    win = ones(Nfft+1,1);
    Ewin = 1;
end
win = win(1:end-1);
% scale_sch = bi2de(repmat([1,0],1,N/2), 'left-msb');
scale_sch = bi2de(repmat([0,0],1,ceil(N/2)), 'left-msb');
if source ~= 0
    create_test_vector(source,Nfft,L,A,alpha,xmin,deltaf,dir,A2);
end
save(parfile);
end

