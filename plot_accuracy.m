clear;
close;
figure(2)
dir = date;
% dir = '04-Mar-2017';
testDir = [dir, '/errAnalysis_s4_hann'];
bartlett = 0;
load([testDir, '/', 'errAnalysis.mat']')

if bartlett
    nList = [1:6]';
    figure(3)
    % plot(nList, 10*log10(sqrt(ubound(:,1))),'color','b');
    hold('on')
    % plot(nList, 10*log10(sqrt(lbound(:,1))),'color','r');
    plot(nList, 10*log10(err(:,1)),'s','color','b');
    set(gca,'xtick',nList)
    % plot(nList, 10*log10(err(:,2)),'d','color','b');
    xticklabels(nList);
    % legend('ubound', 'lbound', 'BFP-scaling', 'No scaling')
    xlabel('Amount of Averaging')
    ylabel('dB')
    ylim([-55 -20])
else    
    nList = [8:13]';
    plot(nList, 10*log10(ubound(:,1)),'color','b');
    set(gca,'xtick',nList)
    hold('on');
    plot(nList, 10*log10(lbound(:,1)),'color','r');
    plot(nList, 10*log10(xferr(:,1)),'s','color','b');
    % plot(nList, 10*log10(xferr(:,2)),'d','color','b');
    xticklabels(num2cell(2.^nList));
    % plot(nList, 10*log10(2*ubound(:,1,3)),'color','r');
    % plot(nList, 10*log10(xferr(:,1,3)),'s','color','r');
    % plot(nList, 10*log10(xferr(:,2,3)),'d','color','r');
    legend('ubound','lbound','BFP-scaling', 'no scaling')
    xlabel('FFT Points')
    ylabel('dB')

    figure(3)
    plot(nList, 10*log10(sqrt(ubound(:,1))),'color','b');
    hold('on')
    % plot(nList, 10*log10(sqrt(lbound(:,1))),'color','r');
    plot(nList, 10*log10(err(:,1)),'s','color','b');
    set(gca,'xtick',nList)
    % plot(nList, 10*log10(err(:,2)),'d','color','b');
    xticklabels(num2cell(2.^nList));
    legend('ubound', 'BFP-scaling', 'No scaling')
    xlabel('FFT Points')
    ylabel('dB')
    % ylim([-55 -10])
end