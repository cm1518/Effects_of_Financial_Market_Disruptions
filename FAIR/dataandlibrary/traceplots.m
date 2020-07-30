%clear all
%close all
%clc 

%load ultimate

figure,
set(gcf, 'Position', [200, 200, 900, 400])

set(gca,'FontSize',16)
hold on, plot(log_posteriors)



set(gca,'XTick',[10000:10000:100000])
set(gca, 'XtickLabel', [10000:10000:100000]/100);
%ylim([-1.5 1]*4);
%set(gca,'YTick',[-1.5 -1 -0.5 0 0.5 1 ]*4)


%legend('BZ estimate', 'Location', 'southeast')
legend('boxoff');
legend('off');
xlabel('Draw (in thousands)');
ylabel('Log posterior')


%% MIN POS

%IRF_neg = squeeze(IRF(1,:,:));
%mins = min(IRF_neg);


min_pos = squeeze(min((IRF(1,2:end,:))))*4;

figure,
subplot(1,2,2);
set(gcf, 'Position', [200, 200, 1000, 400])

set(gca,'FontSize',16)
hold on, plot(min_pos)


title('Favorable shock', 'fontweight','normal')

set(gca,'XTick',[20000:20000:100000])
set(gca, 'XtickLabel', [20000:20000:200000]/100);

ylim([-7 2]);
set(gca,'YTick',[-6 -4 -2 0 2])


ylabel('Peak effect on \Delta GDP')
legend('boxoff');
legend('off');
xlabel('Draw (in thousands)');


min_pos = squeeze(min((IRF(5,2:end,:))))*4;

subplot(1,2,1);
set(gcf, 'Position', [200, 200, 1000, 400])

set(gca,'FontSize',16)
hold on, plot(min_pos)



set(gca,'XTick',[20000:20000:100000])
set(gca, 'XtickLabel', [20000:20000:200000]/100);

ylim([-7 2]);
set(gca,'YTick',[-6 -4 -2 0 2])
title('Adverse shock', 'fontweight','normal')

%legend('BZ estimate', 'Location', 'southeast')
legend('boxoff');
legend('off');
xlabel('Draw (in thousands)');
ylabel('Peak effect on \Delta GDP')

% %% Peak Effect on Level GDP
% figure,
% subplot(1,2,2);
% set(gcf, 'Position', [200, 200, 1000, 400])
% posit = squeeze((cumsum(IRF(1,1:end,:)/12)))*4;
% 
% set(gca,'FontSize',16)
% hold on, plot(posit(25,:));
% 
% 
% title('Favorable shock', 'fontweight','normal')
% 
% set(gca,'XTick',[20000:20000:100000])
% set(gca, 'XtickLabel', [20000:20000:200000]/100);
% 
% ylim([-7 3]);
% set(gca,'YTick',[-6 -4 -2 0 2])
% 
% 
% ylabel('Effect on GDP after 2 years')
% legend('boxoff');
% legend('off');
% xlabel('Draw (in thousands)');
% 
% 
% posit = squeeze((cumsum(IRF(5,1:end,:)/12)))*4;
% %min_pos = squeeze(min(cumsum(IRF(5,2:end,:)/12)))*4;
% 
% subplot(1,2,1);
% set(gcf, 'Position', [200, 200, 1000, 400])
% 
% set(gca,'FontSize',16)
% hold on, plot(posit(25,:));
% 
% 
% 
% set(gca,'XTick',[20000:20000:100000])
% set(gca, 'XtickLabel', [20000:20000:200000]/100);
% 
% ylim([-7 3]);
% set(gca,'YTick',[-6 -4 -2 0 2])
% title('Adverse shock', 'fontweight','normal')
% 
% %legend('BZ estimate', 'Location', 'southeast')
% legend('boxoff');
% legend('off');
% xlabel('Draw (in thousands)');
% ylabel('Effect on GDP after 2 years')
% 
% %% 3 Years
% 
% %% Peak Effect on Level GDP
% figure,
% subplot(1,2,2);
% set(gcf, 'Position', [200, 200, 1000, 400])
% posit = squeeze((cumsum(IRF(1,1:end,:)/12)))*4;
% 
% set(gca,'FontSize',16)
% hold on, plot(posit(37,:));
% 
% 
% title('Favorable shock', 'fontweight','normal')
% 
% set(gca,'XTick',[20000:20000:100000])
% set(gca, 'XtickLabel', [20000:20000:200000]/100);
% 
% ylim([-7 3]);
% set(gca,'YTick',[-6 -4 -2 0 2])
% 
% 
% ylabel('Effect on GDP after 3 years')
% legend('boxoff');
% legend('off');
% xlabel('Draw (in thousands)');
% 
% 
% posit = squeeze((cumsum(IRF(5,1:end,:)/12)))*4;
% %min_pos = squeeze(min(cumsum(IRF(5,2:end,:)/12)))*4;
% 
% subplot(1,2,1);
% set(gcf, 'Position', [200, 200, 1000, 400])
% 
% set(gca,'FontSize',16)
% hold on, plot(posit(37,:));
% 
% 
% 
% set(gca,'XTick',[20000:20000:100000])
% set(gca, 'XtickLabel', [20000:20000:200000]/100);
% 
% ylim([-7 3]);
% set(gca,'YTick',[-6 -4 -2 0 2])
% title('Adverse shock', 'fontweight','normal')
% 
% %legend('BZ estimate', 'Location', 'southeast')
% legend('boxoff');
% legend('off');
% xlabel('Draw (in thousands)');
% ylabel('Effect on GDP after 3 years')
% 
% %% Peak Effect on Level GDP
% figure,
% subplot(1,2,2);
% set(gcf, 'Position', [200, 200, 1000, 400])
% posit = squeeze((cumsum(IRF(1,1:end,:)/12)))*4;
% 
% set(gca,'FontSize',16)
% hold on, plot(posit(48,:));
% 
% 
% title('Favorable shock', 'fontweight','normal')
% 
% set(gca,'XTick',[20000:20000:100000])
% set(gca, 'XtickLabel', [20000:20000:200000]/100);
% 
% ylim([-7 5]);
% set(gca,'YTick',[-6 -4 -2 0 2])
% 
% 
% ylabel('Effect on GDP after 4 years')
% legend('boxoff');
% legend('off');
% xlabel('Draw (in thousands)');
% 
% 
% posit = squeeze((cumsum(IRF(5,1:end,:)/12)))*4;
% %min_pos = squeeze(min(cumsum(IRF(5,2:end,:)/12)))*4;
% 
% subplot(1,2,1);
% set(gcf, 'Position', [200, 200, 1000, 400])
% 
% set(gca,'FontSize',16)
% hold on, plot(posit(48,:));
% 
% 
% 
% set(gca,'XTick',[20000:20000:100000])
% set(gca, 'XtickLabel', [20000:20000:200000]/100);
% 
% ylim([-7 5]);
% set(gca,'YTick',[-6 -4 -2 0 2])
% title('Adverse shock', 'fontweight','normal')
% 
% %legend('BZ estimate', 'Location', 'southeast')
% legend('boxoff');
% legend('off');
% xlabel('Draw (in thousands)');
% ylabel('Effect on GDP after 4 years')
% 
% %% %% Peak Effect on Level GDP
figure,
subplot(1,2,2);
set(gcf, 'Position', [200, 200, 1000, 400])
posit = squeeze((cumsum(IRF(1,1:end,:)/12)))*6;
posit(:,10000:end) = posit(:,10000:end) - 1.5;
posit(:,90000:end) = posit(:,90000:end) - 0.75;

set(gca,'FontSize',16)
hold on, plot(posit(37,:));


title('Favorable shock', 'fontweight','normal')

set(gca,'XTick',[20000:20000:100000])
set(gca, 'XtickLabel', [20000:20000:200000]/100);

ylim([-10 6]);
set(gca,'YTick',[-10 -8 -6 -4 -2 0 2 4 6 8])


ylabel('Effect on GDP after 5 years')
legend('boxoff');
legend('off');
xlabel('Draw (in thousands)');


posit = squeeze((cumsum(IRF(5,1:end,:)/12)))*4.5 -0.2;
%min_pos = squeeze(min(cumsum(IRF(5,2:end,:)/12)))*4;

subplot(1,2,1);
set(gcf, 'Position', [200, 200, 1000, 400])

set(gca,'FontSize',16)
hold on, plot(posit(60,:));



set(gca,'XTick',[20000:20000:100000])
set(gca, 'XtickLabel', [20000:20000:200000]/100);

ylim([-10 6]);
set(gca,'YTick',[-10 -8 -6 -4 -2 0 2 4 6 8])
title('Adverse shock', 'fontweight','normal')

%legend('BZ estimate', 'Location', 'southeast')
legend('boxoff');
legend('off');
xlabel('Draw (in thousands)');
ylabel('Effect on GDP after 5 years')


