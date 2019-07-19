function sections(eptrials, section_boundaries)
%plots section boundaries

hold on

%center of maze
comx = mean([mean(eptrials(eptrials(:,5)==1,2)) mean(eptrials(eptrials(:,5)==2,2))]);
comy = mean([mean(eptrials(eptrials(:,5)==1,3)) mean(eptrials(eptrials(:,5)==2,3))])-100;

%plot rewards
plot(mean(eptrials(eptrials(:,5)==1,2)), mean(eptrials(eptrials(:,5)==1,3)), 'k.', 'markersize', 30);
plot(mean(eptrials(eptrials(:,5)==2,2)), mean(eptrials(eptrials(:,5)==2,3)), 'k.', 'markersize', 30);

    
stem = [comx-30 comx+30 comy-90 comy+60]; %common stem

%rectangle plots determined by above boundaries
recstrt = rectangle('Position', [section_boundaries(1,1), section_boundaries(1,3), (section_boundaries(1,2) - section_boundaries(1,1)),  (section_boundaries(1,4) - section_boundaries(1,3))]);
recstem = rectangle('Position', [stem(1,1), stem(1,3), (stem(1,2) - stem(1,1)),  (stem(1,4) - stem(1,3))]);
recchce = rectangle('Position', [section_boundaries(3,1), section_boundaries(3,3), (section_boundaries(3,2) - section_boundaries(3,1)),  (section_boundaries(3,4) - section_boundaries(3,3))]);
recchmL = rectangle('Position', [section_boundaries(4,1), section_boundaries(4,3), (section_boundaries(4,2) - section_boundaries(4,1)),  (section_boundaries(4,4) - section_boundaries(4,3))]);
recchmR = rectangle('Position', [section_boundaries(5,1), section_boundaries(5,3), (section_boundaries(5,2) - section_boundaries(5,1)),  (section_boundaries(5,4) - section_boundaries(5,3))]);
recrwdL = rectangle('Position', [section_boundaries(6,1), section_boundaries(6,3), (section_boundaries(6,2) - section_boundaries(6,1)),  (section_boundaries(6,4) - section_boundaries(6,3))]);
recrwdR = rectangle('Position', [section_boundaries(7,1), section_boundaries(7,3), (section_boundaries(7,2) - section_boundaries(7,1)),  (section_boundaries(7,4) - section_boundaries(7,3))]);
%reciti = rectangle('Position', [section_boundaries(8,1), section_boundaries(8,3), (section_boundaries(8,2) - section_boundaries(8,1)),  (section_boundaries(8,4) - section_boundaries(8,3))], 'linestyle', '--');


%{
%maze section boundaries [xlow xhigh ylow yhigh]
strt = [comx-50 comx+50  comy-200 comy-90]; %start area
stem = [comx-30 comx+30 comy-90 comy+60]; %common stem
chce = [comx-45 comx+45 comy+60 comy+155]; %choice area
chmL = [comx+45 comx+135 comy+75 comy+155]; %approach arm left
chmR = [comx-135 comx-45 comy+75 comy+155]; %approach arm right
rwdL = [comx+135 comx+200 comy+60 comy+155]; %reward area left
rwdR = [comx-200 comx-135 comy+60 comy+155]; %reward area right
iti = [comx+30 comx+180  comy-170 comy-40]; %intertrial interval platform


%rectangle plots determined by above boundaries
recstrt = rectangle('Position', [strt(1,1), strt(1,3), (strt(1,2) - strt(1,1)),  (strt(1,4) - strt(1,3))]);
recstem = rectangle('Position', [stem(1,1), stem(1,3), (stem(1,2) - stem(1,1)),  (stem(1,4) - stem(1,3))]);
recchce = rectangle('Position', [chce(1,1), chce(1,3), (chce(1,2) - chce(1,1)),  (chce(1,4) - chce(1,3))]);
recchmL = rectangle('Position', [chmL(1,1), chmL(1,3), (chmL(1,2) - chmL(1,1)),  (chmL(1,4) - chmL(1,3))]);
recchmR = rectangle('Position', [chmR(1,1), chmR(1,3), (chmR(1,2) - chmR(1,1)),  (chmR(1,4) - chmR(1,3))]);
recrwdL = rectangle('Position', [rwdL(1,1), rwdL(1,3), (rwdL(1,2) - rwdL(1,1)),  (rwdL(1,4) - rwdL(1,3))]);
recrwdR = rectangle('Position', [rwdR(1,1), rwdR(1,3), (rwdR(1,2) - rwdR(1,1)),  (rwdR(1,4) - rwdR(1,3))]);
reciti = rectangle('Position', [iti(1,1), iti(1,3), (iti(1,2) - iti(1,1)),  (iti(1,4) - iti(1,3))], 'linestyle', '--');
%}

axis([(section_boundaries(7,1) - 20) (section_boundaries(6,2) + 20) (section_boundaries(1,3) - 20) (section_boundaries(5,4) + 20)])
set(gca, 'Xtick',((section_boundaries(7,1) - 20):75:(section_boundaries(6,2) + 20)), 'Ytick',((section_boundaries(1,3) - 20):75:(section_boundaries(5,4) + 20)), 'fontsize', 10)

stem_sections(eptrials)