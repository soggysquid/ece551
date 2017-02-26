% Run a batch of tests and save the results of each test
% No
% clear
dir = date;
% dir = '21-Feb-2017';
seed=123;
runsim=1;
doErrAnalysis=1;
snrAnalysis=0;
pgAnalysis=0;
%%%%%%%%%%%%%%%%%
% Test Settings 
%%%%%%%%%%%%%%%%%
sample_width=16;
A = 1.0; % not used in source 12
alpha = 0.0; % also not used in source 12
if doErrAnalysis
    numsims=10;
    % widthOutList = [0];
    sourceList = [4];
    winList = [0];
    % mList = [8,9,10,11,12,13];
    % avgList = [0];
    % avg=0;
    % alphaList = [0,2^-18,2^-14,2^-10,2^-6,2^-4,2^-2];
    alphaList = [0.001,0.01,0.1];
    % hwList = [0,1];
    % testDir = [date, '/miscAnalysis2'];
    % testDir = [date, '/errAnalysis_s4'];
    if ~exist(testDir, 'dir')
        mkdir(testDir)
    end
    save([testDir, '/', 'settings.mat'], 'numsims', 'widthOutList', 'winList', ...
        'mList', 'avgList', 'hwList', 'alphaList', 'sourceList')
elseif pgAnalysis
    numsims=1;
    widthOutList = [32];
    sourceList = [6];
    winList = [0,2];
    mList = [8,10,12];
    avgList = [0];
    hwList = [0,1];
    testDir = [date, '/pgAnalysis'];
    if ~exist(testDir, 'dir')
        mkdir(testDir)
    end
    save([testDir, '/', 'settings.mat'], 'numsims', 'widthOutList', 'winList', ...
        'mList', 'avgList', 'hwList', 'alpha', 'sourceList')
end    
if ~exist(testDir, 'dir')
    mkdir(testDir)
end
mseIndex = 1;
errdB = zeros(200,1);
maxErr = zeros(200,1);
j=0;
for w=winList 
    for thisWidth = widthOutList
        for N = mList
            for avg = avgList
                L = N+avg;
                for alpha = alphaList
                    A = 1-2^-16-alpha*4;
                    for source = sourceList
                        for hwver = hwList
                            rng(seed);
                            for i = [1:numsims]
                                width_out = thisWidth;
                                Error=0;
                                setup(source, width_out, N, L, w, hwver, A, alpha);
                                
                                matfile = ['s', int2str(source), 'w', int2str(w)...
                                       'N', int2str(N), 'L', int2str(L), 'b', int2str(width_out),...
					  'a', num2str(alpha), 'avg', int2str(avg), 'hw', int2str(hwver), 'iter', int2str(i), '.mat'];
                                load(vectfile);
                                if runsim & ~Error
                                    j = j+1
                                    if hwver == 0
                                        sim('bart_block_scaled.slx')
                                    elseif hwver == 1
                                        sim('bart_block.slx')
                                    else
                                        sim('bart_block_spectsmoothed_scaled.slx')
                                    end
                                    get_sim_results;
                                    save([testDir, '/', matfile], 'Px', 'Px_u', 'xq', 'Xf', 'source', 'win', ...
                                        'Ewin', 'N', 'L', 'width_out', 'alpha', 'hwver', 'blkexp', 'px_latency'); 
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
