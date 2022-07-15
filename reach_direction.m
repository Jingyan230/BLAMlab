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

% plot out
figure(3); clf
subplot(1,2,1); hold on
idx = 1;
for j = 1:2
    for k = 1:2
        avg = endError_all.average.(mental{j}).(corsi{k});
        bim = endError_all.bimanual.(mental{j}).(corsi{k});
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
ylabel(['Endpoint error (' char(0176) ')'])
legend({'No dual task','Corsi','MR','Corsi+MR'}, 'Location', 'northwest')

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