% clear;
dir = date;
% dir = '25-Feb-2017';
testDir = [dir, '/errAnalysis_s4'];
% testDir = [dir, '/errAnalysis_s12_hw1_long'];
load([testDir, '/', 'settings']);
impNum = length(hwList)*length(widthOutList);
impIndex = 0;
if exist('alphaList')
    err = zeros(length(mList),impNum,length(alphaList));
    xferr = zeros(length(mList),impNum,length(alphaList));
    varErr = zeros(length(mList),impNum,length(alphaList));
    maxErr = zeros(length(mList),impNum,length(alphaList));
    xterm = zeros(length(mList),impNum,length(alphaList));
    % ubound = zeros(length(mList),impNum,length(alphaList),length(numsims));
    ubound = zeros(length(mList),impNum,length(alphaList));
    scaling = zeros(length(mList),impNum,length(alphaList));
    expErr = zeros(length(mList),impNum,length(alphaList));
else
    err = zeros(length(mList),impNum);
    XfErr = zeros(length(mList),impNum);
    varErr = zeros(length(mList),impNum);
    expErr = zeros(length(mList),impNum);
    maxErr = zeros(length(mList),impNum);
    scaling = zeros(length(mList),impNum);
    alphaList=alpha;
end
lat = zeros(length(mList),impNum);
w=winList(1);
s=12;
s=4;
% s=sourceList(1);
imp=strings(1,impNum);
blkexpList = zeros(length(mList),impNum);
hwIndex = 0;
for hw=hwList
    for width = widthOutList
        impIndex=impIndex+1;
        alphaIndex = 0;
        for alpha1 = alphaList
            alphaIndex=alphaIndex+1;
            mIndex = 0;
            for avg = avgList
                for M = mList
                    mIndex = mIndex+1;
                    N=M;
                    L=M;
                    Nfft = 2^M;
                    fftl = 2^N;
                    thisXfErr = 0;
                    thisErr = 0;
                    thisExpErr=0;
                    thisVarErr=0;
                    thisMaxErr = 0;   
                    thisUbound = 0;
                    thisXterm = 0;
                    for i = 1:numsims
                        width_out = width;
                        matfile = ['s', int2str(s), 'w', int2str(w)...
                                   'N', int2str(N), 'L', int2str(L), 'b', int2str(width_out),...
                                   'a', num2str(alpha1), 'avg', int2str(avg), 'hw', int2str(hw), 'iter', int2str(i), '.mat'];
                        load([testDir, '/', matfile]);
                        Xf2 = fft(xq, fftl);
%                         lat(mIndex,impIndex) = px_latency;
%                         if hw==1
%                             Px_u = Px_u/2^(2*(13-N));
%                         end
                        
                        Px2_u = Xf2 .* conj(Xf2);
                        Px2 = Px2_u/Nfft^2;
                        % thatErr = thatErr + sum(Px2_u-Px_u)./(Nfft*sum(xq.*conj(xq)));
                        scaling(mIndex,impIndex,alphaIndex) = blkexp;
                        % thisErr = thisErr + sum(Px_u)./(Nfft*sum(xq.*conj(xq)));
                        PSD = sum(Px2_u);
                        % thisErr = thisErr + abs(sum(Px_u)-PSD)/PSD; 
                        % thisErr = thisErr + sum(abs(Px_u-Px2_u))/(sum(xq.*conj(xq))*Nfft);
                        if hw==0
                            Xf = Xf*2^blkexp;
                        end
                        Xfe = real(Xf)-real(Xf2);
                        thisXfErr = thisXfErr + sum(Xfe.*conj(Xfe))/(Nfft*sum(xq.*conj(xq)));
                        Xfe = imag(Xf)-imag(Xf2);
                        thisXfErr = thisXfErr + sum(Xfe.*conj(Xfe))/(Nfft*sum(xq.*conj(xq)));
                        % thisXferr = thisXferr + Nfft*var(Xfe.*conj(Xfe))/(sum(xq.*conj(xq)));
                        thisErr = thisErr + sum(abs(Px2_u-Px_u))/(sum(xq.*conj(xq))*Nfft);
                        thisExpErr = thisExpErr + sum(abs(Xf.*conj(Xf)-Xf2.*conj(Xf2)))...
                            /(sum(xq.*conj(xq))*Nfft);
                        thisVarErr = thisVarErr + (sum((Px-Px2).^2))/PSD;
                        Px = max(Px,2^(-width_out+1));
                        mErr = 10*log10(max(max(Px./Px2),max(Px2./Px)));
                        thisMaxErr = max(mErr,thisMaxErr);
                        K = (1/Nfft)*sum(xq.*conj(xq));
                        % ubound(mIndex,impIndex,alphaIndex,i) = 2^(M+3)*(2^(-2*16)/12)/K
                        thisUbound = thisUbound + 2^(M+3)*(2^(-2*16)/12)/K;
    %                     thisXterm = thisXterm + ...
    %                         (sum(real(Xf).*(real(Xf2)-real(Xf)))...
    %                         + (sum(imag(Xf).*(imag(Xf2)-imag(Xf)))))/(sum(xq.*conj(xq))*Nfft);
                        thisXterm = thisXterm + sum(real(Xf).*real(Xf2)+imag(Xf).*imag(Xf2))/ ...
                                    (sum(xq.*conj(xq))*Nfft);    
                    end
                    xferr(mIndex,impIndex,alphaIndex) = thisXfErr/numsims;
                    xterm(mIndex,impIndex,alphaIndex) = thisXterm/numsims;
                    ubound(mIndex,impIndex,alphaIndex) = thisUbound/numsims;
                    err(mIndex,impIndex,alphaIndex) = thisErr/numsims;
                    experr(mIndex,impIndex,alphaIndex) = thisExpErr/numsims;
                    varErr(mIndex,impIndex,alphaIndex) = thisVarErr/numsims;
                    maxErr(mIndex,impIndex,alphaIndex) = thisMaxErr;
                    blkexpList(mIndex,impIndex,alphaIndex) = blkexp;
                end
            end
        end
        imp(impIndex) = strcat(int2str(hw), int2str(width));
    end
end


                