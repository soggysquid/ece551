clear;
dir = date;
testDir = [date, '/errAnalysis_s4'];
load([testDir, '/', 'settings']);
impIndex = 0;
mIndex = 0;
% parseval = zeros(length(mList),5);
err = zeros(length(mList),length(hwList)*length(widthOutList));
scaling = zeros(length(mList),numsims);
w=winList(1);
s=4;
imp=strings(1,length(hwList)*length(widthOutList));
blkexpList = zeros(length(mList),4);
for hw=hwList
    for width = widthOutList
        impIndex=impIndex+1;
        mIndex = 0;
        for M = mList
            mIndex = mIndex+1;
            N=M;
            L=M;
            fftl = 2^N;
            thisErr = 0;
            thatErr = 0;
            for i = 1:numsims
                width_out = width;
                matfile = ['s', int2str(s), 'w', int2str(w)...
                           'N', int2str(N), 'L', int2str(L), 'b', int2str(width_out),...
                           'a', num2str(alpha), 'hw', int2str(hw), 'iter', int2str(i), '.mat'];
                load(matfile);
                Xf2 = fft(xq, fftl);
                Px2 = Xf2/fftl .* conj(Xf2/fftl);
                thatErr = thatErr + sum(Px2-Px)*fftl./sum(xq.*conj(xq));
                scaling(mIndex,i) = blkexp;
                thisErr = thisErr + sum(Px2-Px)*fftl./sum(xq.*conj(xq));    
            end
            err(mIndex,impIndex) = thisErr/numsims; 
            blkexpList(mIndex,impIndex) = blkexp;
        end
        imp(impIndex) = strcat(int2str(hw), int2str(width));
    end
end


                