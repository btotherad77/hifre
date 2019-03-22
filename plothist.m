function [ output_args ] = plothist( a, ca )
%PLOTHIST Summary of this function goes here
%   Detailed explanation goes here
map = brewermap(3,'Set1'); 
figure
histf(a,-1.3:.01:1.3,'facecolor',map(1,:),'facealpha',.5,'edgecolor','none')
hold on
histf(ca,-1.3:.01:1.3,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none')

box off
axis tight
legalpha('H1','H2','location','northwest')
legend boxoff

end

