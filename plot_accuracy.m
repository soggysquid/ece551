load([dir, '/acc_s12w0b32scale1'])
figure(1);
scatter(Nmat,20*log10(ratio1),'*');
hold('on')
plot(Nmat,20*log10(lubound))
legend('error', 'lower-upper bound', ...
    'Location','northwest');
title('Rms(FFT err)/rms(result) No Scaling');
xlabel('M')
set(gca, 'XTick', Nmat)
ylabel('rms(err)/rms(result) dB')
print([dir, '/fft_noscale_accuracy.png'], '-dpng');

figure(2)
load([dir,'/acc_s12w0b32scale3'])
scatter(Nmat,log10(ratio1),'*');
hold('on')
plot(Nmat,log10(uubound))
legend('error', 'upper-upper bound', ...
    'Location','northwest');
clear title;
title('Rms(FFT err)/rms(result)');
xlabel('M')
set(gca, 'XTick', Nmat)
ylabel('rms(err)/rms(result) dB')
print([dir, '/fft_scale_accuracy.png'], '-dpng');

figure(3)
load([dir, '/acc_s12w0b32scale1'])
scatter(Nmat,10*log10(ratio2),'*');
hold('on')
load([dir,'/acc_s12w0b32scale3'])
scatter(Nmat,10*log10(ratio2),'*');
legend('error w/o scaling', 'error with scaling','Location','northwest')
xlabel('M')
ylabel('(PSD err)/(noise power) dB')
set(gca, 'XTick', Nmat)
title('(PSD err)/(noise power)');
print([dir, '/psd_accuracy.png'], '-dpng');

figure(4)
load([dir, '/acc_s12w0b32scale1'])
scatter(Nmat,lat./(2.^Nmat)');
hold('on')
load([dir, '/acc_s12w0b32scale3'])
scatter(Nmat,lat./(2.^Nmat)');
legend('w/o scaling', 'with scaling')
title('Ratio of latency to FFT points');
xlabel('M');
set(gca, 'XTick', Nmat)
ylabel('cycles/N');
print([dir, '/latency.png'], '-dpng');
