%% minimum angle between endpoint and any target in corsi task

for i = 1:2
    for j = 1:2
        % collecting the data from all subj into dataset
        for m = 1:Nsubj
            dat = data{m}.(hands{i}).(mental{j}).corsi;
            minError_all.(hands{i}).(mental{j})(m) = mean(abs(dat.min_error)) * 180/pi;
        end
    end
end

figure(1); clf; hold on
idx = 1;
for j = 1:2
    avg = minError_all.average.(mental{j});
    bim = minError_all.bimanual.(mental{j});
    error = [avg; bim];
    
    color = col((j-1)*2+2,:);
    plot(idx:idx+1, error, 'Color', [color 0.5], 'HandleVisibility', 'off')
    plot(idx:idx+1, mean(error,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
    
    idx = idx + 2;
end
xticks(1:4)
xticklabels({'Average','Bimanual','Average','Bimanual'})
xtickangle(45)
ylabel(['Minimum direction error (' char(0176) ')'])
legend({'Corsi','Corsi+MR'})

%% error between end angle and target in opposite direction

for i = 1:2
    for j = 1:2
        % collecting the data from all subj into dataset
        for m = 1:Nsubj
            dat = data{m}.(hands{i}).rotation.(corsi{j});
            oppError_all.(hands{i}).(corsi{j})(m) = mean(abs(dat.opposite_error)) * 180/pi;
        end
    end
end

figure(2); clf; hold on
idx = 1;
for j = 1:2
    avg = oppError_all.average.(corsi{j});
    bim = oppError_all.bimanual.(corsi{j});
    error = [avg; bim];
    
    color = col(j+2,:);
    plot(idx:idx+1, error, 'Color', [color 0.5], 'HandleVisibility', 'off')
    plot(idx:idx+1, mean(error,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
    
    idx = idx + 2;
end
xticks(1:4)
xticklabels({'Average','Bimanual','Average','Bimanual'})
xtickangle(45)
ylabel(['Error from opposite target (' char(0176) ')'])
legend({'MR','Corsi+MR'})