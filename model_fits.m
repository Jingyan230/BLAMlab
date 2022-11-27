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

%%
rot(1:32,1) = "no_rot";
rot(33:64,1) = "rot";
cor([1:16 33:48], 1) = "no_corsi";
cor([17:32 49:64], 1) = "corsi";
subject = repmat((1:16)', [4 1]);

y = [sd_all.average.no_rotation.no_corsi'; sd_all.average.no_rotation.corsi';
    sd_all.average.rotation.no_corsi'; sd_all.average.rotation.corsi'];
T = table(rot, cor, subject, y, 'VariableNames', {'rotation','corsi','subject','sd'});
writetable(T,'C:/Users/Chris/Documents/R/corsirotation/data/sd1.csv')

%%

map(1:32,1) = "average";
map(33:64,1) = "bimanual";
cor([1:16 33:48], 1) = "no_corsi";
cor([17:32 49:64], 1) = "corsi";
subject = repmat((1:16)', [4 1]);

y = [sd_all.average.no_rotation.no_corsi'; sd_all.average.no_rotation.corsi';
    sd_all.bimanual.no_rotation.no_corsi'; sd_all.bimanual.no_rotation.corsi'];
T = table(map, cor, subject, y, 'VariableNames', {'mapping','corsi','subject','sd'});
writetable(T,'C:/Users/Chris/Documents/R/corsirotation/data/sd2.csv')

%%

map(1:32,1) = "average";
map(33:64,1) = "bimanual";
rot([1:16 33:48], 1) = "no_rotation";
rot([17:32 49:64], 1) = "rotation";
subject = repmat((1:16)', [4 1]);

y = [sd_all.average.no_rotation.no_corsi'; sd_all.average.rotation.no_corsi';
    sd_all.bimanual.no_rotation.no_corsi'; sd_all.bimanual.rotation.no_corsi'];
T = table(map, rot, subject, y, 'VariableNames', {'mapping','rotation','subject','sd'});
writetable(T,'C:/Users/Chris/Documents/R/corsirotation/data/sd3.csv')
%% 
figure(10); clf; hold on
plot(1, sd_all.bimanual.no_rotation.no_corsi, '.k')
plot(2, sd_all.average.rotation.no_corsi, '.k')
plot(3, sd_all.bimanual.rotation.no_corsi, '.k')

plot(5, sd_all.bimanual.no_rotation.no_corsi, '.k')
plot(6, sd_all.average.no_rotation.corsi, '.k')
plot(7, sd_all.bimanual.no_rotation.corsi, '.k')

plot(9, sd_all.average.rotation.no_corsi, '.k')
plot(10, sd_all.average.no_rotation.corsi, '.k')
plot(11, sd_all.average.rotation.corsi, '.k')

% plot(1, sd_all.average.rotation.no_corsi, '.k')
% plot(2, sd_all.average.no_rotation.corsi, '.k')
% plot(3, sd_all.bimanual.no_rotation.no_corsi, '.k')
% 
% plot(5, sd_all.bimanual.rotation.no_corsi, '.k')
% plot(6, sd_all.bimanual.no_rotation.corsi, '.k')
% plot(7, sd_all.average.rotation.corsi, '.k')

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
subplot(2,1,1); hold on
for i = 1:4
    plot(2*(i-1)+1:2*(i-1)+2, sd_test(:,[i i+4])', '.', 'MarkerSize', 15, 'Color', col(i,:), 'HandleVisibility', 'off')
    plot(2*(i-1)+1:2*(i-1)+2, sd_test(:,[i i+4])', 'Color', col(i,:), 'HandleVisibility', 'off')
    plot(2*(i-1)+1:2*(i-1)+2, mean(sd_test(:,[i i+4]), 1, 'omitnan'), '-ko', 'MarkerFaceColor', col(i,:), 'MarkerSize', 8, 'LineWidth', 1)
end
xlim([0.5 8.5])
xticks(1:8)
xticklabels({'Average','Bimanual','Average','Bimanual','Average','Bimanual','Average','Bimanual'})
ylabel('Circular standard deviation')
legend({'No dual task','Corsi','MR','Corsi+MR'},'Location','northwest')

subplot(2,1,2); hold on
for i = 1:4
    plot(2*(i-1)+1:2*(i-1)+2, 1 - vm_wt_test(:,[i i+4])', '.', 'MarkerSize', 15, 'Color', col(i,:), 'HandleVisibility', 'off')
    plot(2*(i-1)+1:2*(i-1)+2, 1 - vm_wt_test(:,[i i+4])', 'Color', col(i,:), 'HandleVisibility', 'off')
    plot(2*(i-1)+1:2*(i-1)+2, 1 - mean(vm_wt_test(:,[i i+4]), 1, 'omitnan'), '-ko', 'MarkerFaceColor', col(i,:), 'MarkerSize', 8, 'LineWidth', 1)
end
xlim([0.5 8.5])
xticks(1:8)
xticklabels({'Average','Bimanual','Average','Bimanual','Average','Bimanual','Average','Bimanual'})
ylabel('Proportion of random reaches')

print('C:\Users\Chris\Dropbox\Conferences\SFN 2022\model_initDir','-dpdf','-painters')

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

%% model comparison by BIC

corsi_bic = data{1}.average.no_rotation.corsi.unif_bic;
no_corsi_bic = data{1}.average.no_rotation.no_corsi.unif_bic;

figure(6); clf; hold on
plot([0.5 4.5], [no_corsi_bic no_corsi_bic], 'k')
plot([4.5 8.5], [corsi_bic corsi_bic], 'k')
idx = 1;
for k = 1:2
    for i = 1:2
        for j = 1:2
            % collecting the data from all subj into dataset
            for m = 1:Nsubj
                dat = data{m}.(hands{i}).(mental{j}).(corsi{k});
                vm_bic(m) = dat.vm_bic;
            end
            
            unif_bic = dat.unif_bic;
            plot(idx, vm_bic, '.r', 'MarkerSize', 10)
            idx = idx + 1;
        end
    end
end
t_name2 = {'avg no con task','avg MR','bi MR','bi no con task','avg corsi','avg corsi+MR','bi corsi','bi corsi+MR'};
xticks(1:8)
xticklabels(t_name2)
xtickangle(45)
ylabel('BIC')
