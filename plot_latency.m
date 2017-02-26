clear;
nList = [8:13]';
load('latency.mat');
% plot(nList, lat./2.^repmat(nList,[1,2]), 'o');
bar(nList, lat./2.^repmat(nList,[1,2]));
set(gca,'xtick',nList)
xlabel('FFT Points');
ylabel('Latency per Sample');
legend('BFP Scaling', 'No Scaling')
ylim([2 2.6])
% xticklabels(int2str(2.^nList));
xticklabels(num2cell(2.^nList));


