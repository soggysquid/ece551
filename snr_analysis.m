rng=123;
b = 16;         % ADC width, signed
Delta = 4/2^b;
alpha = sqrt(Delta^2/12);   % Set alpha to variance of quantization noise
width_out = 0;   % width_out 0 is full-precision
N = 10;
L = 10;
w = 0;
source = 4;
scaling = 0;
setup(source, width_out, N, L, w, scaling, alpha)
load(parfile)
load(vectfile)
matfile = ['s', int2str(source), 'w', int2str(w)...
        'N', int2str(N), 'L', int2str(L), 'widthout', int2str(width_out),...
        'alpha', num2str(alpha), 'scale', int2str(scaling)];
sim('bart_block');
get_sim_results;
save([dir, '/', matfile]);