% Run a batch of tests and save the results of each test
clear
dir = date;
numSims=1; % Number simulations to do
seed=123;
%%%%%%%%%%%%%%%%%
% Test Settings %
%%%%%%%%%%%%%%%%%
widthOutList = [0,32];   
sample_width=16;
sourceList = [4,6];
winList = [0,2]; 
A = 0.9
alpha = .0167;
mList = [8,10,12];
avgList = [0,1,2]
hwList = [0,1,2];
testDir = [date, '/batch_tests'];
if ~exist(dir, 'dir')
    mkdir(dir)
end
for w=winList 
    for width_out = widthOutList
        for M = mlist
            for avg = avgList
                P = M+avg;
                for source = sourceList
                    for hw = hwList
                        matfile = ['s', int2str(source), 'w', int2str(w)...
                                   'M', int2str(M), 'P', int2str(P), 'b', int2str(width_out),...
                                   'a', num2str(alpha), 'hw', int2str(hw),'.mat'];
                        rng(seed);
                        setup(source, width_out, M, P, w, scaling, A, alpha)
                        if hwver == 0
                            sim('bart_block_scaled.slx')
                        elseif hwver == 1
                            sim('bart_block.slx')
                        else
                            sim('bart_block_spectsmoothed_scaled.slx')
                        end
                        get_sim_results;
                        save([thisDir, '/', matfile]);                        
                    end
                end
            end
        end
    end
end