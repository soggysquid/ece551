% Script to setup all the parameters and test vectors that we're going to
% use for modelling and analysis
% clearvars -except latency sample_width twiddle_width;
function setup(width_out, source, N, L, w, alpha, deltaf, scaling)
% sample-width: fft sample width (adc width is fixed)
% width_out: width of output
% source: data set to use
% source:
%   0: real .dat
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
%   2: hamming
%   3: blackman
% alpha: noise variance
% scaling: 
%   1: full precision
%   2: scale periodogram
%   3: block floating point scaling

switch nargin
    case 0
        % set width_out to 0 to use full precision
        width_out = 32;
        source = 12;
        N = 8;
        L = 8;
        w = 0;
        % alpha = 2^-15;
        alpha = 0.0;
        deltaf = 5/2^N;
        scaling = 1;
end
if ~exist('alpha')
    alpha = 0;
end
sources = {' iq .dat ', ' usrp ', ...
    ' slot ', ' 1 tone in bin ', ' 2 tones in bins ', ' 1 inter-bin tone ',...
    ' 2 inter-bin tones ', 'impulse', 'all ones', 'all zeros'}; 
wins = {'rect', 'kaiser', 'hamming', 'blackman'};
Nmax = 13; % log2 of max points of run-time configurable FFT
bartmax = Nmax-3;  % max number segments we can average 
latency = 25000;
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
% alpha = 0.0;
usemex = 0;
% Figure out max bit growth for full-precision, unscaled
 % output width of abs value (unsigned) of fft (full-precision unscaled), 
width_fft = sample_width+Nmax;  % this is the unsigned width
bp_fft = sample_width-1;
bp_reinterp = bp_fft+Nmax;
% width_perio = width_fft*2-1; % full precision  
if scaling == 1
    width_perio = width_fft*2;
    bp_perio = bp_reinterp*2;  % full precision
elseif scaling == 3
    width_fft = sample_width;
    bp_fft = sample_width-1;
    width_perio = width_fft;
    bp_perio = bp_fft;  
    width_out = 0;
else
    width_perio = width_fft;
    bp_perio = bp_reinterp;  
end
width_fpadd = width_perio+bartmax;
bp_fpadd = bp_perio;
bp_bart = bp_fpadd+bartmax;   % for full precision bart output
if width_out
    % bp_out = bp_bart-(width_fpadd-width_out);
    bp_out = bp_fpadd-(width_fpadd-width_out)+bartmax;
    % bp_out = bp_fpadd-(width_fpadd-width_out);
else
    width_out = width_fpadd;
    bp_out = bp_bart;
end

%% Create a test vector
if source == 0 || source == 1
    fs = 48e6;
    fc = fs/2;
    data = openBin(datfile,inf,'int16'); 
    % .dat files are screwed up right now
    xr = data(2:4:end);
    xi = data(4:4:end);
    if source == 0
        x = xr/2^13;
    else
        xr = xr/2^15;
        xi = xi/2^15;
        x = xr + 1j*xi;
    end
    L = floor(log2(length(x)));
    P = 2^L; 
    x = x(1:P,:);
    nsect = length(x)/Nfft;
    % r = mod(length(x),Nfft);
    % P = length(x)-r;
    ncols = 1;
    x = reshape(x,[Nfft,nsect,ncols]);
elseif source == 2  % Use data from csv file
    fs = 48e6;
    xr = csvread(refile);
    xi = csvread(imfile);
    x = xr' + 1j*xi';
    L = floor(log2(length(x)));
    % P = floor(length(x)/2^L)*2^L;
    P = 2^L;
    x = x(1:P,:);
    s = size(x);
    ncols = s(2);
    % x = [x;zeros(Nfft*(ceil((latency+P)/Nfft))-length(x),ncols)];  
    % x = [x;zeros(ceil((latency+rec_length)/2^L)*2^L-rec_length, ncols)];
    fc = fs*[1:ncols]-fs/2;
    nsect = length(x)/Nfft;    
    x = reshape(x,[Nfft,nsect,ncols]);
elseif source == 3   % Slot test
    % Figure out why there is a DC value (has to do with the way I'm
    % rounding)
    P = 2^L;
    fs = Nfft*10e3;
    fc = 0;
    k0 = 1/64*Nfft;
    x = slot(Nfft, k0);
%    X = ones(1,Nfft);
%    X(1:k0) = 0;
%    X(Nfft-(k0+1):Nfft) = 0;
%    x = ifft(X);
    nsect = P/Nfft;
    x = repmat(x, [1,nsect]); 
elseif source >= 4 & source < 8  % tones
    % 4: one tone in bin
    % 5: two tones in bin
    % 6: one tone between bins
    % Oct 30 change to be continuous
    P = 2^L;
    nsect = P/Nfft;
    % n = 0:Nfft;
    n = 0:P-1;
    fs = 10e3*Nfft;
    fc = 0;
    if source == 4 | source == 5
        A1 = 1;
        A1 = 1-2^-15;
        A2 = max(2^-15,alpha);
        % A2 = 2^-15;
        f1 = 1/8;
        f2 = f1+deltaf;     
        x = A1*cos(2*pi*n*f1) + 1j*A1*sin(2*pi*n*f1);
        if source == 5
            x = x + A2*(cos(2*pi*n*f2)+1j*sin(2*pi*n*f2));
        end
        x = x/max(max(abs(x)));
    elseif source == 6
        A1 = 1;
        f1 = 1/8 + 0.5/Nfft;
        x = A1*(cos(2*pi*n*f1)+1j*sin(2*pi*n*f1));
        % x = A1*cos(2*pi*n*f1);
    elseif source == 7
        A1 = 1;
        A2 = max(2^-8,alpha);
        f1 = 1/8 + 0.5/Nfft;
        f2 = f1+deltaf;
        x = A1*(cos(2*pi*n*f1) + 1j*sin(2*pi*n*f1))+A2*(cos(2*pi*n*f2) + 1j*sin(2*pi*n*f2));
    end
    x = x + alpha*randn(size(x));
    % x = repmat(x(1:end-1)',[1,nsect]);
    x = reshape(x,[Nfft,nsect]);
elseif source == 8
    fs = Nmax*10e3;
    fc = 0;
    n=[1:Nfft];
    P = 2^L;
    x = zeros(1,Nfft);
    x(1) = 1;
    nsect = P/Nfft;
    ncols = 1;
    x = repmat(x', [1,nsect]); 
elseif source == 9 | source == 10 | source == 12
    fs = Nmax*10e3;
    A1 = 0.5;
    fc = 0;
    n=[1:Nfft];
    P = 2^L;
    if source == 9
        x = ones(1,Nfft)*A1 + 1j*ones(1,Nfft)*A1;
    elseif source == 12
        x = rand(1,Nfft) + 1j*rand(1,Nfft);
    else
        x = zeros(1,Nfft);
    end
    ncols = 2;
    nsect = P/Nfft;
    % x = x + 1j*x;
    x = repmat(x', [1,nsect,ncols]);     
elseif source == 11
    fs = Nmax*10e3;
    fc = 0;
    k = [1:10];
    A = 1/10;
    freqs = 1./(2.^k);
    P=2^L;
    n = [0:P-1];
    x = zeros(size(n));
    nsect = P/Nfft;
    for f=freqs
        x = x + A*(cos(2*pi*n*f)+1j*sin(2*pi*n*f));
    end
    x = reshape(x,[Nfft,nsect]);
end
% Quantize source if it is not coming form adc
%if source > 2 & usemex
   % x = 1/(2^adc_width-1) * round(x * (2^adc_width-1));
%    xr = fxquant(real(x), adc_width, rmode, lmode);
%    xi = fxquant(imag(x), adc_width, rmode, lmode);
%    x = xr + 1j*xi;
% end
xorig = x;
x = x + alpha*randn(size(x)) + 1j*alpha*randn(size(x));
%% Set up windowing
if w == 1
    win = kaiser(Nfft+1,9);
    Ewin = sum(kaiser(Nmax+1,9).^2)/Nmax;
elseif w == 2
    win = hamming(Nfft+1);
    Ewin = sum(hamming(Nmax+1).^2)/Nmax;
elseif w == 3
    win = blackman(Nfft+1);
    Ewin = sum(blackman(Nmax+1).^2)/Nmax;
else
    win = ones(Nfft+1,1);
    Ewin = 1;
end
win = win(1:end-1);
%% Set up models
if usemex    % We're doing a c-model simulation
    scale_sch = repmat([3],1,floor(N/2));
    if mod(N,2)
        scale_sch = [scale_sch,1];
    end
    generics = struct('C_NFFT_MAX', 12);
    generics.C_ARCH = 3;
    generics.C_HAS_NFFT = 1;
    if float
        generics.C_USE_FLT_PT = 1;
        generics.C_INPUT_WIDTH = 32; 
    else
        generics.C_USE_FLT_PT = 0;
        generics.C_INPUT_WIDTH = sample_width;
    end
    generics.C_TWIDDLE_WIDTH = twiddle_width;
    generics.C_HAS_SCALING = has_scaling;
    generics.C_HAS_BFP = 0;
    generics.C_HAS_ROUNDING = 1;
else
    % scale_sch = bi2de(repmat([1,0],1,N/2), 'left-msb');
    scale_sch = bi2de(repmat([0,0],1,ceil(N/2)), 'left-msb');
end
findex = [1:2^N]'*(fs)/(2^N+1);
xlength = length(x(:));
delay = 0;
nfft_ts = timeseries(N*ones(xlength,1), [1:xlength]);
nfft_valid = zeros(xlength,1);
nfft_valid(1:2) = 1;
nfft_valid_ts = timeseries(nfft_valid, [1:xlength]);
% Prepend a frame of zeros because first sample of first frame is dropped
% by fft core which might munge the first frame
xr_ts = timeseries([real(x(1)); real(x(:))], [1:xlength+1]); 
xi_ts = timeseries([imag(x(1)); imag(x(:))], [1:xlength+1]);
xv_ts = timeseries(ones(xlength+2,1),[1:xlength+2]);
Ts = 1/fs;
save(parfile);
end

