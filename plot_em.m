function plot_em( X1,pl_title, b, X2 )
%UNTITLED2 Summary of this function goes here
%   b is width of data
epsilon1 = 10*log10(2^-57);
epsilon2 = 10*log10(2^(-b+1));
epsilon3 = 10*log10(2^-255);

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
%         plot([1:length(X1)],max(10*log10(X1(:,:,b)), epsilon1), 'r', ...
%             [1:length(X2)],max(10*log10(X2(:,:,b)), epsilon2), 'b');
        plot([1:length(X1)],10*log10(X1(:,:,b)*1e3), 'r', ...
            [1:length(X2)],10*log10(X2(:,:,b)*1e3), 'b');
    else
%         plot(max(10*log10(X1(:,:,b)), epsilon3), 'b');
        plot(10*log10(X1(:,:,b)), 'b');
    end
    hold('on')
%    if b==1
%        xlim([findex(1,1) findex(bands,end)]);
%    end
end
xlabel('Bins')
ylabel('dBm')
xlim([0 length(X1)]);
if exist('X2')
    legend('SW', 'HW')
end
title(pl_title);
hold('off')
end