%% Figures
[I M]=max(log_posteriors);

%% Burn-in after inspection of log posteriors
draws1 = draws(:,3000:end);

nr = size(draws1,2);
resp_1_3_VAR=store_responses(1,3,:);
resp_2_3_VAR=store_responses(2,3,:);
resp_3_3_VAR=store_responses(3,3,:);
resp_4_3_VAR=store_responses(4,3,:);
clear IRF_p IRF_n IRF clear y_pos error1_pos error2_pos IRFm IRFl IRFu;

% Generate IRs for each draw
for j=0:3
for i = 1:nr
    for k=1:60
    IRF_p(1+j,k+1,i) = draws1(37+j,i)*exp(-(k-draws1(45+j,i))^2/draws1(53+j,i)) +  draws1(41+j,i)*exp(-(k-draws1(49+j,i))^2/draws1(57+j,i));
    end
    IRF_p(1+j,1,i) = draws1(33+j ,i);
end
end


for j=0:3
for i = 1:nr
    for k=1:60
    IRF_n(1+j,k+1,i) = draws1(9+j,i)*exp(-(k-draws1(17+j,i))^2/draws1(25+j,i)) +  draws1(13+j,i)*exp(-(k-draws1(21+j,i))^2/draws1(29+j,i));
    end
    IRF_n(1+j,1,i) = draws1(5+j,i);
end
end

IRF = [IRF_n; IRF_p]; % First 4 are negative, then positive IRFs
IRFu = prctile(IRF,100-(100-90)/2,3);
IRFu68 = prctile(IRF,100-(100-68)/2,3);
IRFm = prctile(IRF,50,3);
IRFl = prctile(IRF,(100-90)/2,3);
IRFl68 = prctile(IRF,(100-68)/2,3);



figure,
set(gcf, 'Position', [200, 200, 1000, 800])

x = [1:61];
subplot(2,2,4);
set(gca,'FontSize',16)
hold on,

y_pos = 1/12*cumsum(IRFm(1,:)) / max(IRFm(3,:));
error1_pos = 1/12*cumsum((IRFm(1,:)) - (IRFl(1,:))) / max(IRFm(3,:));
error2_pos = 1/12*cumsum((IRFu(1,:)) - (IRFm(1,:))) / max(IRFm(3,:));
hold on, plot(y_pos(1:end-1), 'blue', 'LineWidth', 2);
H= shadedErrorBar(x,y_pos,[error2_pos(1:end); error1_pos(1:end)], 'blue')
error1_pos = 1/12*cumsum((IRFm(1,:)) - (IRFl68(1,:))) / max(IRFm(3,:));
error2_pos = 1/12*cumsum((IRFu68(1,:)) - (IRFm(1,:))) / max(IRFm(3,:));
H= shadedErrorBar68(x,y_pos,[error2_pos(1:end); error1_pos(1:end)], 'blue')
hold on, plot(0*IRFm(1,:)/12, '--', 'Color', [0.7,0.7,0.7], 'LineWidth', 1);
hold on, plot(y_pos(1:end), 'blue', 'LineWidth', 2);
xlim([1 61]);
set(gca,'XTick',[12 24 36 48 60 ] +1)
set(gca, 'XtickLabel', [12 24 36 48 60 ]*1/12);
ylim([-1.5 1]*4);
legend('boxoff');
legend('off');
xlabel('Years');
ylabel('IP (ppt)')


subplot(2,2,3);
set(gcf, 'Position', [200, 200, 800, 900])

set(gca,'FontSize',16)
hold on,
ylabel('IP (ppt)')

y_pos = 1/12*cumsum(IRFm(5,:))/ max(IRFm(7,:));
error1_pos = 1/12*cumsum(IRFm(5,:))/ max(IRFm(7,:)) - 1/12*cumsum(IRFl(5,:))/ max(IRFm(7,:));
error2_pos = 1/12*cumsum(IRFu(5,:))/ max(IRFm(7,:)) - 1/12*cumsum(IRFm(5,:))/ max(IRFm(7,:));
hold on, plot(y_pos, 'red', 'LineWidth', 2);
%legend('BZ estimate', 'Location', 'northwest')
H= shadedErrorBar(x,y_pos,[error2_pos(1:end); error1_pos(1:end)], 'red')
error1_pos = 1/12*cumsum(IRFm(5,:))/ max(IRFm(7,:)) - 1/12*cumsum(IRFl68(5,:))/ max(IRFm(7,:));
error2_pos = 1/12*cumsum(IRFu68(5,:))/ max(IRFm(7,:)) - 1/12*cumsum(IRFm(5,:))/ max(IRFm(7,:));
H= shadedErrorBar68(x,y_pos,[error2_pos(1:end); error1_pos(1:end)], 'red')

hold on, plot(0*IRFm(1,:)/12, '--', 'Color', [0.7,0.7,0.7], 'LineWidth', 1);
hold on, plot(y_pos, 'red', 'LineWidth', 2);
xlim([1 61]);
ylim([-1.5 1]*4);
set(gca,'XTick',[12 24 36 48 60 ] +1)
set(gca, 'XtickLabel', [12 24 36 48 60 ]*1/12);
legend('boxoff');
legend('off');
xlabel('Years');



%%
subplot(2,2,2);
set(gca,'FontSize',16)
hold on,
title('Favorable Shock', 'fontweight','normal');
y_pos = (IRFm(3,:))/ max(IRFm(3,:));
error1_pos = (IRFm(3,:))/ max(IRFm(3,:)) - (IRFl(3,:))/ max(IRFm(3,:));
error2_pos = (IRFu(3,:))/ max(IRFm(3,:)) - (IRFm(3,:))/ max(IRFm(3,:));
hold on, plot(y_pos(1:end), 'blue', 'LineWidth', 2);
H= shadedErrorBar(x,y_pos,[error2_pos(1:end); error1_pos(1:end)], 'blue')
error1_pos = (IRFm(3,:))/ max(IRFm(3,:)) - (IRFl68(3,:))/ max(IRFm(3,:));
error2_pos = (IRFu68(3,:))/ max(IRFm(3,:)) - (IRFm(3,:))/ max(IRFm(3,:));
H= shadedErrorBar68(x,y_pos,[error2_pos(1:end); error1_pos(1:end)], 'blue')
hold on, plot(0*IRFm(1,:), '--', 'Color', [0.7,0.7,0.7], 'LineWidth', 1);
hold on, plot(y_pos(1:end), 'blue', 'LineWidth', 2);
xlim([1 61]);
set(gca,'XTick',[12 24 36 48 60 ] +1)
set(gca, 'XtickLabel', [12 24 36 48 60 ]*1/12);
ylim([-0.05 0.31]*4);
set(gca, 'Ytick', [0 0.25 0.5 0.75 1 1.25]);
legend('boxoff');
legend('off');
xlabel('Years');
ylabel('EBP (ppt)')


subplot(2,2,1);
set(gcf, 'Position', [200, 200, 1000, 800])
set(gca,'FontSize',16)
hold on,
ylabel('EBP (ppt)')
y_pos = (IRFm(7,:))/ max(IRFm(7,:));
error1_pos = (IRFm(7,:))/ max(IRFm(7,:)) - (IRFl(7,:))/ max(IRFm(7,:));
error2_pos = (IRFu(7,:))/ max(IRFm(7,:)) - (IRFm(7,:))/ max(IRFm(7,:));
hold on, plot(y_pos(1:end), 'red', 'LineWidth', 2);
H= shadedErrorBar(x,y_pos,[error2_pos(1:end); error1_pos(1:end)], 'red')
error1_pos = (IRFm(7,:))/ max(IRFm(7,:)) - (IRFl68(7,:))/ max(IRFm(7,:));
error2_pos = (IRFu68(7,:))/ max(IRFm(7,:)) - (IRFm(7,:))/ max(IRFm(7,:));
H= shadedErrorBar68(x,y_pos,[error2_pos(1:end); error1_pos(1:end)], 'red')
hold on, plot(0*IRFm(1,:)/12, '--', 'Color', [0.7,0.7,0.7], 'LineWidth', 1);
hold on, plot(y_pos(1:end), 'red', 'LineWidth', 2);
xlim([1 61]);
set(gca,'XTick',[12 24 36 48 60 ] +1)
set(gca, 'XtickLabel', [12 24 36 48 60 ]*1/12);
legend('boxoff');
set(gca, 'Ytick', [0 0.25 0.5 0.75 1 1.25]);
ylim([-0.05 0.31]*4);
xlabel('Years');
title('Adverse Shock', 'fontweight','normal');
legend('off');


