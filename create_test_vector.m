function [x,findex,fs] = create_test_vector(source, Nfft, L, alpha, xmin, deltaf, dir)
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
% Nfft: length of fft
% L: log2(length of data)
% alpha: should be sigma, guassian noise source
% xmin: min amplitude of second tone if there is one
% deltaf: of second tone if there is one
% dir: where to save the vect file

% if source == 0 || source == 1
%     fs = 48e6;
%     fc = fs/2;
%     data = openBin(datfile,inf,'int16'); 
%     % .dat files are screwed up right now
%     xr = data(2:4:end);
%     xi = data(4:4:end);
%     if source == 0
%         x = xr/2^13;
%     else
%         xr = xr/2^15;
%         xi = xi/2^15;
%         x = xr + 1j*xi;
%     end
%     L = floor(log2(length(x)));
%     P = 2^L; 
%     x = x(1:P,:);
%     nsect = length(x)/Nfft;
%     % r = mod(length(x),Nfft);
%     % P = length(x)-r;
%     ncols = 1;
%     x = reshape(x,[Nfft,nsect,ncols]);
% elseif source == 2  % Use data from csv file
%     fs = 48e6;
%     xr = csvread(refile);
%     xi = csvread(imfile);
%     x = xr' + 1j*xi';
%     L = floor(log2(length(x)));
%     % P = floor(length(x)/2^L)*2^L;
%     P = 2^L;
%     x = x(1:P,:);
%     s = size(x);
%     ncols = s(2);
%     % x = [x;zeros(Nfft*(ceil((latency+P)/Nfft))-length(x),ncols)];  
%     % x = [x;zeros(ceil((latency+rec_length)/2^L)*2^L-rec_length, ncols)];
%     fc = fs*[1:ncols]-fs/2;
%     nsect = length(x)/Nfft;    
%     x = reshape(x,[Nfft,nsect,ncols]);
% elseif source == 3   % Slot test
if ~exist('Nmax')
    Nmax = 13;
end
vectfile = [dir, '/testvector'];
assignin('base', 'vectfile', vectfile)
fs = 10e3*Nmax;
if source == 3
    % Figure out why there is a DC value (has to do with the way I'm
    % rounding)
    P = 2^L;
    % fs = Nfft*10e3;
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
    ncols = 1;
    % n = 0:Nfft;
    n = 0:P-1;
    % fs = 10e3*Nfft;
    fc = 0;
    if source == 4 | source == 5
        A1 = 1;
        A1 = 1-xmin;
        A2 = max(xmin,alpha);
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
    % fs = Nmax*10e3;
    fc = 0;
    n=[1:Nfft];
    P = 2^L;
    x = zeros(1,Nfft);
    x(1) = 1;
    nsect = P/Nfft;
    ncols = 1;
    x = repmat(x', [1,nsect]); 
elseif source == 9 | source == 10 | source == 12
    % fs = Nmax*10e3;
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
    ncols = 1;
    nsect = P/Nfft;
    % x = x + 1j*x;
    x = repmat(x', [1,nsect,ncols]);     
elseif source == 11
    % fs = Nmax*10e3;
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
findex = [1:Nfft]'*(fs)/(Nfft+1);
xorig = x;
x = x + alpha*randn(size(x)) + 1j*alpha*randn(size(x));
xlength = length(x(:));
nfft_ts = timeseries(log2(Nfft)*ones(xlength,1), [1:xlength]);
nfft_valid = zeros(xlength,1);
nfft_valid(1:2) = 1;
nfft_valid_ts = timeseries(nfft_valid, [1:xlength]);
% Prepend a frame of zeros because first sample of first frame is dropped
% by fft core which might munge the first frame
xr_ts = timeseries([real(x(1)); real(x(:))], [1:xlength+1]); 
xi_ts = timeseries([imag(x(1)); imag(x(:))], [1:xlength+1]);
xv_ts = timeseries(ones(xlength+1,1),[1:xlength+1]);
% xv_ts = timeseries([ones(3,1);zeros(xlength-3,1)],[1:xlength]);
Ts = 1/fs;
save(vectfile)