figure(2)
load('s4error_feb27')
nList = [8:13]';
plot(nList, 10*log10(ubound(:,1)),'color','b');
set(gca,'xtick',nList)
hold('on');
plot(nList, 10*log10(xferr(:,1)),'s','color','b');
plot(nList, 10*log10(xferr(:,3)),'d','color','b');
xticklabels(num2cell(2.^nList));
% plot(nList, 10*log10(2*ubound(:,1,3)),'color','r');
% plot(nList, 10*log10(xferr(:,1,3)),'s','color','r');
% plot(nList, 10*log10(xferr(:,2,3)),'d','color','r');

legend('ubound','BFP-scaling', 'no scaling')
xlabel('FFT Points')
ylabel('dB')

figure(3)
plot(nList, 10*log10(err(:,1)/2^N),'s','color','b');
set(gca,'xtick',nList)
hold('on')
plot(nList, 10*log10(err(:,3)/2^N),'d','color','b');
% plot(nList, 10*log10(err(:,1,3)),'s','color','r');
% plot(nList, 10*log10(err(:,2,3)),'d','color','r');
xticklabels(num2cell(2.^nList));
legend('BFP-scaling', 'No scaling')
xlabel('FFT Points')
ylabel('dB')

ylim([-80 -45])
