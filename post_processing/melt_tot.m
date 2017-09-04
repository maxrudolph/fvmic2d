clear; close all; fclose all;
totmelt = [3.9658 3.4982 2.6493 2.1486 1.8202 0.7506];
config = {'Decoupling no root','Reference case','Small root','Thickened root w/ decoupling','Thickened root','Thickened root deep decoupling'};
% config = {'Narrow 15 km root' 'No thickened root' 'Reference case' '5 km root' '10 km root' '15 km root' 'No decoupling 15 km root' 'Deeper decoupling'};
figure;
bar(totmelt);
% bar(1,totmelt(1),'red'); hold on
% bar(2,totmelt(2),'blue'); hold on
% bar(3,totmelt(3),'black'); hold on
% bar(4,totmelt(4),'green'); hold on
% bar(5,totmelt(5),'yellow'); hold on
% bar(6,totmelt(6),'m');hold on
set(gca,'XTickLabel',config,'FontSize', 12, 'XTickLabelRotation',45)
ylabel('Total Melt production (m^3/m/yr)','FontSize',12);
set(gca,'FontSize',12)
set(gcf, 'Color', 'w');
% export_fig test.png -m2