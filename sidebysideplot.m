function sidebysideplot(ax1, ax2, title1, title2, yspan, xspan1, xspan2);
% Copying a figure into another
% get figure by selecting and the command ax1=gca;
% test1.fig and test2.fig are the names of the figure files which you would % like to copy into multiple subplots
h3 = figure; %create new figure
s1 = subplot(1,2,1); %create and get handle to the subplot axes
s2 = subplot(1,2,2);
fig1 = get(ax1,'children'); %get handle to all the children in the figure
fig2 = get(ax2,'children');
copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
copyobj(fig2,s2);
set(s1, 'Ylim', yspan, 'Xlim', xspan1);
ylabel(s1, 'PSD dBm'); xlabel(s1, 'bins'); title(s1, title1);
set(s2, 'Ylim', yspan, 'Xlim', xspan2);
xlabel(s2, 'bins'); title(s2, title2);
legend('Matlab','Simulink')

