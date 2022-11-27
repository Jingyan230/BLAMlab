%% probability density of reach direction error from one participant

subj = 3; % subject to be analyzed

figure(1); clf; hold on
title(subj_name{subj})
pts = -180:.1:180;
for i = 1:2
    for j = 1:2
        for k = 1:2
            s_error = data{subj}.(hands{i}).(mental{j}).(corsi{k}).initDir_error * 180 / pi;
            fks = ksdensity(s_error, pts,'Bandwidth',10);
            if i == 1
                plot(pts, fks, 'LineWidth', 2, 'color', col(2*(j-1)+k,:));
            else
                plot(pts, fks, '--', 'LineWidth', 2, 'color', col(2*(j-1)+k,:));
            end
        end
    end
end
xlabel('Reach direction error')
ylabel('Probability density')
xticks(-180:90:180)
axis([-180 180 0 .05])
legend(t_name)
set(gca,'TickDir','out')
hold off

%% reach direction vs target direction from one participant

idx = 1;
figure(2); clf; hold on
for i = 1:2
    for j = 1:2
        for k = 1:2
            target_angle = data{subj}.(hands{i}).(mental{j}).(corsi{k}).targetAngle * 180/pi;
            reach_angle = data{subj}.(hands{i}).(mental{j}).(corsi{k}).initDir * 180/pi;
            subplot(2,4,idx); hold on
            plot([-185 185], [-185 185], 'k')
            plot(target_angle,reach_angle,'.r','MarkerSize',10)
            if i == 2
                xlabel('Target direction')
            end
            
            if idx == 1
                title('No concurrent task')
                ylabel('Initial cursor direction')
            elseif idx == 2
                title('Corsi')
            elseif idx == 3
                title('MR')
            elseif idx == 4
                title('Corsi + MR')
            elseif idx == 5
                ylabel('Initial cursor direction')
            end
            
            daspect([1 1 1])
            axis([-185 185 -185 185])
            xticks(-180:90:180)
            yticks(-180:90:180)
            set(gca,'TickDir','out')
            box off
            
            idx = idx + 1;
        end
    end
end
hold off

%% reach direction vs target direction pooled across all participants

bins = -180:30:180;
Nbins = length(bins) - 1;
counts = NaN(Nbins);
idx = 1;
clims = [0 1];

col1 = [1 1 1];
col2 = [1 0 0];
Nstep = 100;
map = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2)...
    ,Nstep)', linspace(col1(3),col2(3),Nstep)'];

figure(10); clf
for i = 1:2
    for j = 1:2
        for k = 1:2
            Ntrials = length(data{1}.(hands{i}).(mental{j}).(corsi{k}).targetAngle);
            targetAngle = NaN(Ntrials, Nsubj);
            initDir = NaN(Ntrials, Nsubj);
            for m = 1:Nsubj
                dat = data{m}.(hands{i}).(mental{j}).(corsi{k});
                
                targetAngle(:,m) = dat.targetAngle * 180/pi;
                initDir(:,m) = dat.initDir * 180/pi;
                
            end
            for m = 1:Nbins
                idx2 = (targetAngle >= bins(m)) == (targetAngle < bins(m+1));
                total = sum(idx2,'all'); % sum the number of trials in bin
                
                counts(:,m) = histcounts(initDir(idx2), bins) ./ total;
            end
            
            subplot(2,4,idx); hold on
            imagesc(counts, clims)
            plot([0 13],[0 13],'k')
            
            colormap(map)
            
            if idx == 1
                title('No concurrent task')
                ylabel('Initial cursor direction')
            elseif idx == 2
                title('Corsi')
            elseif idx == 3
                title('MR')
            elseif idx == 4
                title('Corsi + MR')
            elseif idx == 5
                ylabel('Initial cursor direction')
            end
            
            if i == 2
                xlabel('Target direction')
            end
            
            xticks(0.5:6:12.5)
            yticks(0.5:6:12.5)
            xticklabels(-180:180:180)
            yticklabels(-180:180:180)
            axis([0.5 Nbins+0.5 0.5 Nbins+0.5])
            axis square
            set(gca,'TickDir','out')
            box off
            idx = idx + 1;
        end
    end
end

print('C:\Users\Chris\Dropbox\Conferences\SFN 2022\heatmap','-dpdf','-painters')

%%

idx = 1;
figure(20); clf
for i = 1:2
    for j = 1:2
        for k = 1:2
            Ntrials = length(data{1}.(hands{i}).(mental{j}).(corsi{k}).targetAngle);
            targetAngle = NaN(Ntrials, Nsubj);
            initDir = NaN(Ntrials, Nsubj);
            for m = 1:Nsubj
                dat = data{m}.(hands{i}).(mental{j}).(corsi{k});
                
                targetAngle(:,m) = dat.targetAngle * 180/pi;
                initDir(:,m) = dat.initDir * 180/pi;
                
            end
%             for m = 1:Nbins
%                 idx2 = (targetAngle >= bins(m)) == (targetAngle < bins(m+1));
%                 total = sum(idx2,'all'); % sum the number of trials in bin
%                 
%                 counts(:,m) = histcounts(initDir(idx2), bins) ./ total;
%             end
            
            subplot(2,4,idx); hold on
            plot([-180 180],[-180 180],'k')
            scatter(targetAngle(:), initDir(:), 5, 'r', 'filled', 'MarkerFaceAlpha', 0.4)
                        
            if idx == 1
                title('No concurrent task')
                ylabel('Initial cursor direction')
            elseif idx == 2
                title('Corsi')
            elseif idx == 3
                title('MR')
            elseif idx == 4
                title('Corsi + MR')
            elseif idx == 5
                ylabel('Initial cursor direction')
            end
            
            if i == 2
                xlabel('Target direction')
            end
            
            xticks(-180:180:180)
            yticks(-180:180:180)
            axis([-180 180 -180 180])
            axis square
            set(gca,'TickDir','out')
            box off
            idx = idx + 1;
        end
    end
end

print('C:\Users\Chris\Dropbox\Conferences\SFN 2022\raster','-dpdf','-painters')

%%
idx = 1;
figure(11); clf
for i = 1:2
    for j = 1:2
        for k = 1:2
            Ntrials = length(data{1}.(hands{i}).(mental{j}).(corsi{k}).targetAngle);
            initDir_error = NaN(Ntrials, Nsubj);
            for m = 1:Nsubj
                dat = data{m}.(hands{i}).(mental{j}).(corsi{k});
                
                initDir_error(:,m) = dat.initDir_error * 180/pi;
            end
            
            subplot(2,4,idx); hold on
            histogram(initDir_error, -180:10:180)
            
            if idx == 1
                title('No concurrent task')
                ylabel('Counts')
            elseif idx == 2
                title('Sequence')
            elseif idx == 3
                title('MR')
            elseif idx == 4
                title('Sequence + MR')
            elseif idx == 5
                ylabel('Counts')
            end
            
            if i == 2
                xlabel('Reach direction error')
            end
            
            xticks(-180:90:180)
            yticks(0:100:300)
            axis([-180 180 0 350])
            set(gca,'TickDir','out')
            box off
            idx = idx + 1;
        end
    end
end

print('C:\Users\Chris\Dropbox\Conferences\SFN 2022\histogram','-dpdf','-painters')

%% collect kinematic data into separate variables

clear endError_all initError_all pathLength_all reachTime_all
for i = 1:2
    for j = 1:2
        for k = 1:2
            % collecting the data from all subj into dataset
            for m = 1:Nsubj
                dat = data{m}.(hands{i}).(mental{j}).(corsi{k});
                
                endError_all.(hands{i}).(mental{j}).(corsi{k})(m) = mean(abs(dat.endpoint_error)) * 180/pi;
                initError_all.(hands{i}).(mental{j}).(corsi{k})(m) = mean(abs(dat.initDir_error),'omitnan') * 180/pi;
                
                pathLength_all.(hands{i}).(mental{j}).(corsi{k})(m) = mean(dat.pathLength);
                reachTime_all.(hands{i}).(mental{j}).(corsi{k})(m) = mean(dat.timeReach) / 1000;
            end
        end
    end
end


map(1:64, 1) = "average";
map(65:128, 1) = "bimanual";
rot([1:32 65:96],1) = "no_rot"; 
rot([33:64 97:128],1) = "rot";
cor([1:16 33:48 65:80 97:112], 1) = "no_corsi";
cor([17:32 49:64 81:96 113:128], 1) = "corsi";
subject = repmat((1:16)', [8 1]);

y = [endError_all.average.no_rotation.no_corsi'; endError_all.average.no_rotation.corsi';
    endError_all.average.rotation.no_corsi'; endError_all.average.rotation.corsi';
    endError_all.bimanual.no_rotation.no_corsi'; endError_all.bimanual.no_rotation.corsi';
    endError_all.bimanual.rotation.no_corsi'; endError_all.bimanual.rotation.corsi'];
T = table(map, rot, cor, subject, y, 'VariableNames', {'mapping','rotation','corsi','subject','sd'});
writetable(T,'C:/Users/Chris/Documents/R/corsirotation/data/endpoint.csv')
% 
% 
% y = [pathLength_all.average.no_rotation.no_corsi'; pathLength_all.average.no_rotation.corsi';
%     pathLength_all.average.rotation.no_corsi'; pathLength_all.average.rotation.corsi';
%     pathLength_all.bimanual.no_rotation.no_corsi'; pathLength_all.bimanual.no_rotation.corsi';
%     pathLength_all.bimanual.rotation.no_corsi'; pathLength_all.bimanual.rotation.corsi'];
% T = table(map, rot, cor, subject, y, 'VariableNames', {'mapping','rotation','corsi','subject','path_length'});
% writetable(T,'C:/Users/Chris/Documents/R/corsirotation/data/path_length.csv')

y = [initError_all.average.no_rotation.no_corsi'; initError_all.average.no_rotation.corsi';
    initError_all.average.rotation.no_corsi'; initError_all.average.rotation.corsi';
    initError_all.bimanual.no_rotation.no_corsi'; initError_all.bimanual.no_rotation.corsi';
    initError_all.bimanual.rotation.no_corsi'; initError_all.bimanual.rotation.corsi'];
T = table(map, rot, cor, subject, y, 'VariableNames', {'mapping','rotation','corsi','subject','sd'});
writetable(T,'C:/Users/Chris/Documents/R/corsirotation/data/initError.csv')
%%

% plot out
figure(3); clf
subplot(1,2,2); hold on
idx = 1;
for j = 1:2
    for k = 1:2
        avg = initError_all.average.(mental{j}).(corsi{k});
        bim = initError_all.bimanual.(mental{j}).(corsi{k});
        error = [avg; bim];
        
        color = col((j-1)*2+k,:);
        plot(idx:idx+1, error, 'Color', [color 0.5], 'HandleVisibility', 'off')
        plot(idx:idx+1, mean(error,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
        
        idx = idx + 2;
    end
end
xticks(1:8)
xticklabels({'Average','Bimanual','Average','Bimanual','Average','Bimanual','Average','Bimanual'})
xtickangle(45)
axis([0.5 8.5 0 120])
ylabel(['Initial cursor direction error (' char(0176) ')'])
legend({'No concurrent task','Seq','MR','Seq+MR'}, 'Location', 'northwest')

% subplot(1,2,2); hold on
% idx = 1;
% for j = 1:2
%     for k = 1:2
%         avg = initError_all.average.(mental{j}).(corsi{k});
%         bim = initError_all.bimanual.(mental{j}).(corsi{k});
%         error = [avg; bim];
%         
%         color = col((j-1)*2+k,:);
%         plot(idx:idx+1, error, 'Color', [color 0.5])
%         plot(idx:idx+1, mean(error,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
%         
%         idx = idx + 2;
%     end
% end
% xticks(1:8)
% xticklabels({'Average','Bimanual','Average','Bimanual','Average','Bimanual','Average','Bimanual'})
% xtickangle(45)
% ylabel(['Initial direction error (' char(0176) ')'])

subplot(1,2,1); hold on
idx = 1;
for j = 1:2
    for k = 1:2
        no_cor = initError_all.(hands{j}).(mental{k}).no_corsi;
        cor = initError_all.(hands{j}).(mental{k}).corsi;
        error = [no_cor; cor];
        
        color = col((j-1)*2+k,:);
        plot(idx:idx+1, error, 'Color', [color 0.5], 'HandleVisibility', 'off')
        plot(idx:idx+1, mean(error,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
        
        idx = idx + 2;
    end
end
xticks(1:8)
xticklabels({'No seq', 'Seq', 'No seq', 'Seq', 'No seq', 'Seq', 'No seq', 'Seq'})
xtickangle(45)
axis([0.5 8.5 0 120])
ylabel(['Initial cursor direction error (' char(0176) ')'])
legend({'Avg','Avg + MR','Bim','Bim + MR'}, 'Location', 'northwest')

print('C:\Users\Chris\Dropbox\Conferences\SFN 2022\endpoint','-dpdf','-painters')

%% 

figure(8); clf


subplot(1,2,2); hold on
idx = 1;
for j = 1:2
    for k = 1:2
        avg = initError_all.average.(mental{j}).(corsi{k});
        bim = initError_all.bimanual.(mental{j}).(corsi{k});
        error = [avg; bim];
        
        color = col((j-1)*2+k,:);
        plot(idx:idx+1, error, 'Color', [color 0.5])
        plot(idx:idx+1, mean(error,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
        
        idx = idx + 2;
    end
end
xticks(1:8)
xticklabels({'Average','Bimanual','Average','Bimanual','Average','Bimanual','Average','Bimanual'})
xtickangle(45)
ylabel(['Initial direction error (' char(0176) ')'])

% print('C:\Users\Chris\Dropbox\Conferences\SFN 2022\endpoint','-dpdf','-painters')

%%

figure(4); clf;
subplot(1,2,1); hold on
idx = 1;
for j = 1:2
    for k = 1:2
        avg = pathLength_all.average.(mental{j}).(corsi{k});
        bim = pathLength_all.bimanual.(mental{j}).(corsi{k});
        pathLength = [avg; bim];
        
        color = col((j-1)*2+k,:);
        plot(idx:idx+1, pathLength, 'Color', [color 0.5])
        plot(idx:idx+1, mean(pathLength,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
        
        idx = idx + 2;
    end
end
xticks(1:8)
xticklabels({'Average','Bimanual','Average','Bimanual','Average','Bimanual','Average','Bimanual'})
xtickangle(45)
ylabel('Path length (m)')

subplot(1,2,2); hold on
idx = 1;
for j = 1:2
    for k = 1:2
        avg = reachTime_all.average.(mental{j}).(corsi{k});
        bim = reachTime_all.bimanual.(mental{j}).(corsi{k});
        reachTime = [avg; bim];
        
        color = col((j-1)*2+k,:);
        plot(idx:idx+1, reachTime, 'Color', [color 0.5])
        plot(idx:idx+1, mean(reachTime,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
        
        idx = idx + 2;
    end
end
xticks(1:8)
xticklabels({'Average','Bimanual','Average','Bimanual','Average','Bimanual','Average','Bimanual'})
xtickangle(45)
ylabel('Movement duration (s)')