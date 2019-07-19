function stem_sections(eptrials)

%maze section plots. Better if posplot or evtposplot is already be in a figure.

hold on

%plot rewards
%plot(mean(eptrials(eptrials(:,5)==1,2)), mean(eptrials(eptrials(:,5)==1,3)), 'k.', 'markersize', 30);
%plot(mean(eptrials(eptrials(:,5)==2,2)), mean(eptrials(eptrials(:,5)==2,3)), 'k.', 'markersize', 30);

%center of maze
comx = mean([mean(eptrials(eptrials(:,5)==1,2)) mean(eptrials(eptrials(:,5)==2,2))]);
comy = mean([mean(eptrials(eptrials(:,5)==1,3)) mean(eptrials(eptrials(:,5)==2,3))])-100;


%stem section boundaries [xlow xhigh ylow yhigh]
stem1 = [comx-30 comx+30 comy-90 comy-52.5];
stem2 = [comx-30 comx+30 comy-52.5 comy-15];
stem3 = [comx-30 comx+30 comy-15 comy+22.5];
stem4 = [comx-30 comx+30 comy+22.5 comy+60];


%rectangle plots determined by above boundaries
recstem1 = rectangle('Position', [stem1(1,1), stem1(1,3), (stem1(1,2) - stem1(1,1)),  (stem1(1,4) - stem1(1,3))], 'linestyle', '-');
recstem2 = rectangle('Position', [stem2(1,1), stem2(1,3), (stem2(1,2) - stem2(1,1)),  (stem2(1,4) - stem2(1,3))], 'linestyle', '-');
recstem3 = rectangle('Position', [stem3(1,1), stem3(1,3), (stem3(1,2) - stem3(1,1)),  (stem3(1,4) - stem3(1,3))], 'linestyle', '-');
recstem4 = rectangle('Position', [stem4(1,1), stem4(1,3), (stem4(1,2) - stem4(1,1)),  (stem4(1,4) - stem4(1,3))], 'linestyle', '-');


%hold off

%axis([150 650 60 455])
%set(gca, 'Xtick',(150:75:600), 'Ytick',(30:75:435), 'fontsize', 10)