% For a given confugration (scaled or unscaled), test with different params 
runSim=0; % Set to 0 if simulations have already been run and you need to load files
windows = [0,2,3];  % rect, hann, Black
widths = [0]; % precision of output
sample_width = 16;         % ADC width, signed
fftlength = [8,10,12];
averaging = [0];
sourceList = [4, 6];
Delta = 2/2^16;
alphaList = [sqrt(2^10*Delta)];
scaling = 1;

thisDir = [date, '/full_precision_perio'];
thisDir = ['04-Feb-2017/batch_tests'
if ~exist(dir, 'dir')
    mkdir(dir)
end
for w=windows 
    for width_out = widths
        for M = fftlength
            for avg = averaging
                for alpha = alphaList
                    for source = sourceList
                        for hwver = hwList
                            L = M-avg;
                            matfile = ['s', int2str(source), 'w', int2str(w)...
                                       'N', int2str(N), 'L', int2str(L), 'b', int2str(width_out),...
                                       'a', num2str(alpha), 'hw', int2str(scaling),'.mat'];
                            if runSim
                                rng=123;
                                setup(source, width_out, M, L, w, scaling, alpha)
                                load(parfile)
                                load(vectfile)
                                sim('bart_block');
                                get_sim_results;
                                save([thisDir, '/', matfile]);
                            else
                                load([thisDir, '/', matfile]);
                            end
                        if source ==4
                            if w == 0
                                bw = 0.89*2*pi/N;
                            elseif w == 2
                                deltaw = 1.44*2*pi/N;
                            else
                                deltaw = 1.68*2*pi/N;
                            end
                            [maxPx, imax] = max(Px);
                            fc = findex(imax);
                            mask = find(findex<(fc-bw/2) | findex>(fc+bw/2));
                            varPxx = var(Px(mask));
                            meanPxx = mean(Px(mask));                            
                    end
                end
            end
        end
    end
end