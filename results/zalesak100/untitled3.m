%a=load('A5.dat');
%scatter(a(:,1), a(:,2),11,'filled', 's','k')
%xlim([0,1])
%ylim([0,1])
%axis square
%box on

a=load('X5.dat');
scatter(a(:,1), a(:,2),11,a(:,3),'filled', 's')
xlim([0,1])
ylim([0,1])
box on