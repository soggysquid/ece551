figure(2)
load('s4error_2')
nList = [8:13]';
plot(nList, 10*log10(2*ubound(:,1,1)),'color','b');
set(gca,'xtick',nList)
hold('on');
plot(nList, 10*log10(xferr(:,1,1)),'s','color','b');
plot(nList, 10*log10(xferr(:,2,1)),'d','color','b');
xticklabels(num2cell(2.^nList));
plot(nList, 10*log10(2*ubound(:,1,3)),'color','r');
plot(nList, 10*log10(xferr(:,1,3)),'s','color','r');
plot(nList, 10*log10(xferr(:,2,3)),'d','color','r');

legend('ubound, \sigma=0.001','\sigma=0.001, scaling', '\sigma=0.001, no scaling',...
    'ubound, \sigma=0.1', '\sigma=0.1, scaling', '\sigma=0.1, no scaling')
xlabel('FFT Points')
ylabel('dB')

figure(3)
plot(nList, 10*log10(err(:,1,1)),'s','color','b');
set(gca,'xtick',nList)
hold('on')
plot(nList, 10*log10(err(:,2,1)),'d','color','b');
plot(nList, 10*log10(err(:,1,3)),'s','color','r');
plot(nList, 10*log10(err(:,2,3)),'d','color','r');
xticklabels(num2cell(2.^nList));
legend('\sigma=0.001, scaling', '\sigma=0.001, no scaling',...
    'ubound, \sigma=0.1', '\sigma=0.1, scaling', '\sigma=0.1, no scaling')
xlabel('FFT Points')
ylabel('dB')