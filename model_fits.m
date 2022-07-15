%% model fits plotted by handedness condition

clear mu_all sd_all unif_wt_all initError_all
for i = 1:2
    for j = 1:2
        for k = 1:2
            % collecting the data from all subj into dataset
            for m = 1:Nsubj
                dat = data{m}.(hands{i}).(mental{j}).(corsi{k});
                
                mu_all.(hands{i}).(mental{j}).(corsi{k})(m) = dat.mu;
                sd_all.(hands{i}).(mental{j}).(corsi{k})(m) = dat.sd;
                unif_wt_all.(hands{i}).(mental{j}).(corsi{k})(m) = dat.unif_weight;
            end
        end
    end
end

figure(1); clf
subplot(1,2,1); hold on
idx = 1;
for j = 1:2
    for k = 1:2
        avg = sd_all.average.(mental{j}).(corsi{k});
        bim = sd_all.bimanual.(mental{j}).(corsi{k});
        sd = [avg; bim];
        
        color = col((j-1)*2+k,:);
        plot(idx:idx+1, sd, 'Color', [color 0.5], 'HandleVisibility', 'off')
        plot(idx:idx+1, mean(sd,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
        
        idx = idx + 2;
    end
end
xticks(1:8)
xticklabels({'Average','Bimanual','Average','Bimanual','Average','Bimanual','Average','Bimanual'})
xtickangle(45)
ylabel(['Circular standard deviation (' char(0176) ')'])
legend({'No dual task','Corsi','MR','Corsi+MR'},'Location','northwest')

subplot(1,2,2); hold on
idx = 1;
for j = 1:2
    for k = 1:2
        avg = unif_wt_all.average.(mental{j}).(corsi{k});
        bim = unif_wt_all.bimanual.(mental{j}).(corsi{k});
        unif_wt = [avg; bim];
        
        color = col((j-1)*2+k,:);
        plot(idx:idx+1, unif_wt, 'Color', [color 0.5])
        plot(idx:idx+1, mean(unif_wt,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
        
        idx = idx + 2;
    end
end
xticks(1:8)
xticklabels({'Average','Bimanual','Average','Bimanual','Average','Bimanual','Average','Bimanual'})
xtickangle(45)
ylabel('Proportion of random reaches')

%% model fits plotted by mental rotation condition

figure(2); clf
subplot(1,2,1); hold on
idx = 1;
for i = 1:2
    for k = 1:2
        no_mr = sd_all.(hands{i}).no_rotation.(corsi{k});
        mr = sd_all.(hands{i}).rotation.(corsi{k});
        sd = [no_mr; mr];
        
        color = col((i-1)*2+k,:);
        plot(idx:idx+1, sd, 'Color', [color 0.5], 'HandleVisibility', 'off')
        plot(idx:idx+1, mean(sd,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
        
        idx = idx + 2;
    end
end
xticks(1:8)
xticklabels({'No MR','MR','No MR','MR','No MR','MR','No MR','MR'})
xtickangle(45)
ylabel(['Circular standard deviation (' char(0176) ')'])
legend({'avg','avg Corsi','bim','bim Corsi'})

subplot(1,2,2); hold on
idx = 1;
for i = 1:2
    for k = 1:2
        no_mr = unif_wt_all.(hands{i}).no_rotation.(corsi{k});
        mr = unif_wt_all.(hands{i}).rotation.(corsi{k});
        unif_wt = [no_mr; mr];
        
        color = col((i-1)*2+k,:);
        plot(idx:idx+1, unif_wt, 'Color', [color 0.5])
        plot(idx:idx+1, mean(unif_wt,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
        
        idx = idx + 2;
    end
end
xticks(1:8)
xticklabels({'No MR','MR','No MR','MR','No MR','MR','No MR','MR'})
xtickangle(45)
ylabel('Proportion of random reaches')

%% model fits plotted by corsi condition

figure(3); clf
subplot(1,2,1); hold on
idx = 1;
for i = 1:2
    for j = 1:2
        no_cor = sd_all.(hands{i}).(mental{j}).no_corsi;
        cor = sd_all.(hands{i}).(mental{j}).corsi;
        sd = [no_cor; cor];
        
        color = col((i-1)*2+j,:);
        plot(idx:idx+1, sd, 'Color', [color 0.5], 'HandleVisibility', 'off')
        plot(idx:idx+1, mean(sd,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
        
        idx = idx + 2;
    end
end
xticks(1:8)
xticklabels({'No Corsi','Corsi','No Corsi','Corsi','No Corsi','Corsi','No Corsi','Corsi'})
xtickangle(45)
ylabel(['Circular standard deviation (' char(0176) ')'])
legend({'avg','avg MR','bim','bim MR'})

subplot(1,2,2); hold on
idx = 1;
for i = 1:2
    for j = 1:2
        no_cor = unif_wt_all.(hands{i}).(mental{j}).no_corsi;
        cor = unif_wt_all.(hands{i}).(mental{j}).corsi;
        unif_wt = [no_cor; cor];
        
        color = col((i-1)*2+j,:);
        plot(idx:idx+1, unif_wt, 'Color', [color 0.5])
        plot(idx:idx+1, mean(unif_wt,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
        
        idx = idx + 2;
    end
end
xticks(1:8)
xticklabels({'No Corsi','Corsi','No Corsi','Corsi','No Corsi','Corsi','No Corsi','Corsi'})
xtickangle(45)
ylabel('Proportion of random reaches')

%% model fits after excluding non-recoverable kappa and weights
threshold = [NaN NaN NaN NaN 6 4 3 2 2 2 1];

idx = 1;
for i = 1:2
    for j = 1:2
        for k = 1:2
            % collecting the data from all subj into dataset
            for m = 1:Nsubj
                dat = data{m}.(hands{i}).(mental{j}).(corsi{k});
                
                vm_weight = 1 - dat.unif_weight;
                
                good_kappa = double(dat.kappa > threshold);
                good_weight = floor(vm_weight*10) + 1;
                
                good_kappa(good_weight) = good_kappa(good_weight) + 1;
                
                % data{10}.average.rotation.corsi.sd has bad fit but this
                % procedure already excludes this data
                
                if sum(good_kappa == 2)
                    sd_test(m,idx) = dat.sd;
                else
                    sd_test(m,idx) = NaN;
                end
                
                if vm_weight < 0.4
                    vm_wt_test(m,idx) = NaN;
                else
                    vm_wt_test(m,idx) = vm_weight;
                end
            end
            idx = idx + 1;
            
        end
    end
end

figure(4); clf
subplot(1,2,1); hold on
for i = 1:4
    plot(2*(i-1)+1:2*(i-1)+2, sd_test(:,[i i+4])', '.', 'MarkerSize', 15, 'Color', col(i,:), 'HandleVisibility', 'off')
    plot(2*(i-1)+1:2*(i-1)+2, sd_test(:,[i i+4])', 'Color', col(i,:), 'HandleVisibility', 'off')
    plot(2*(i-1)+1:2*(i-1)+2, mean(sd_test(:,[i i+4]), 1, 'omitnan'), '-o', 'Color', col(i,:), 'MarkerFaceColor', col(i,:), 'MarkerSize', 8, 'LineWidth', 3)
end
xlim([0.5 8.5])
xticks(1:8)
xticklabels({'Average','Bimanual','Average','Bimanual','Average','Bimanual','Average','Bimanual'})
ylabel('Circular standard deviation')
legend({'No dual task','Corsi','MR','Corsi+MR'},'Location','northwest')

subplot(1,2,2); hold on
for i = 1:4
    plot(2*(i-1)+1:2*(i-1)+2, 1 - vm_wt_test(:,[i i+4])', '.', 'MarkerSize', 15, 'Color', col(i,:), 'HandleVisibility', 'off')
    plot(2*(i-1)+1:2*(i-1)+2, 1 - vm_wt_test(:,[i i+4])', 'Color', col(i,:), 'HandleVisibility', 'off')
    plot(2*(i-1)+1:2*(i-1)+2, 1 - mean(vm_wt_test(:,[i i+4]), 1, 'omitnan'), '-o', 'Color', col(i,:), 'MarkerFaceColor', col(i,:), 'MarkerSize', 8, 'LineWidth', 3)
end
xlim([0.5 8.5])
xticks(1:8)
xticklabels({'Average','Bimanual','Average','Bimanual','Average','Bimanual','Average','Bimanual'})
ylabel('Proportion of random reaches')

%% check model fits for each subject

subj = 10;

delt = pi/32; % interval for plotting PDF
points = -pi:delt:pi-delt; % vector of points to plot PDF
vmPDF = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa))); % PDF of von Mises distribution

figure(5); clf
idx = 1;
for i = 1:2
    for j = 1:2
        for k = 1:2
            subplot(2,4,idx); hold on
            dat = data{subj}.(hands{i}).(mental{j}).(corsi{k});
            mu = dat.mu;
            kappa = dat.kappa;
            weight = 1 - dat.unif_weight;

            initError = dat.initDir_error;
            
            pdf = vmPDF(points, mu, kappa); % PDF of von Mises distribution
            mixPDF = weight * pdf + (1-weight) * (1 / (2*pi)); % weight von Mises with uniform distribution

            histogram(initError, points,'Normalization','pdf');
            plot(points, mixPDF, 'LineWidth', 2)
            ylim([0 5])
            title(t_name{idx})
            idx = idx + 1;
        end
    end
end