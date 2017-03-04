clear;
close;
figure(2)
dir = date;
testDir = [dir, '/errAnalysis_s12'];
load([testDir, '/', 'errAnalysis.mat']')
nList = [8:13]';
plot(nList, 10*log10(ubound(:,1)),'color','b');
set(gca,'xtick',nList)
hold('on');
plot(nList, 10*log10(xferr(:,1)),'s','color','b');
plot(nList, 10*log10(xferr(:,2)),'d','color','b');
xticklabels(num2cell(2.^nList));
% plot(nList, 10*log10(2*ubound(:,1,3)),'color','r');
% plot(nList, 10*log10(xferr(:,1,3)),'s','color','r');
% plot(nList, 10*log10(xferr(:,2,3)),'d','color','r');

legend('ubound','BFP-scaling', 'no scaling')
xlabel('FFT Points')
ylabel('dB')

figure(3)
plot(nList, 10*log10(err(:,1)),'s','color','b');
set(gca,'xtick',nList)
hold('on')
plot(nList, 10*log10(err(:,2)),'d','color','b');
% plot(nList, 10*log10(err(:,1,3)),'s','color','r');
% plot(nList, 10*log10(err(:,2,3)),'d','color','r');
xticklabels(num2cell(2.^nList));
legend('BFP-scaling', 'No scaling')
xlabel('FFT Points')
ylabel('dB')

ylim([-55 -20])
