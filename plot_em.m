function plot_em( X1,findex, pl_title, fc, b, X2 )
%UNTITLED2 Summary of this function goes here
%   b is width of data
epsilon = 20*log10(2^-b);
s = size(X1);
bands = 2;
if length(s) == 3
    bands = s(3);
else
    bands = 1;
    X1 = reshape(X1,[s,1]);
    if exist('X2')
        X2 = reshape(X2,[s,1]);
    end
end
%findex = repmat(findex(:),[1,bands]) + ...
%    repmat(fc, [length(findex),1]);
gca();
for b=1:bands
    if exist('X2')
        % change this to if fc=0, use findex(:,1) otherwise findex(:,b) 
        plot(findex(:,1)/1e6, max(20*log10(X1(:,:,b)), epsilon), 'r', ...
            findex(:,1)/1e6, max(20*log10(X2(:,:,b)), epsilon), 'b');
    else
        plot(findex(:,1)/1e6, max(20*log10(X1(:,:,b)), epsilon), 'b');
    end
    hold('on')
%    if b==1
%        xlim([findex(1,1) findex(bands,end)]);
%    end
end
xlabel('Frequency MHz')
ylabel('dB')
if exist('X2')
    legend('SW', 'HW')
end
title(pl_title);
hold('off')
end