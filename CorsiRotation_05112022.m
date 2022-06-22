clear all
close all
rng(2);
tic
delayStep = 20;
space = 6;

% loading data

% The data files were stored under the folder 'Data/subject/subject/'.
% The data of each subject were stored under different folder, the name of
% the folder is listed as 'subj_name'. Under the folder of each subject,
% the data from each task were stored in different folder, named as
% 'avg-nC-nR', 'avg-nC-R', 'avg-C-nR', 'avg-C-R', 'bi-nC-nR', 'bi-nC-R',
% 'bi-C-nR', 'bi-C-R' for different tasks. 'avg' as the average mapping,
% 'bi' as the bimanual mapping, 'C''nC' refers to Corsi or no-Corsi,
% 'R''nR' refers to Rotation or no-Rotation.

folder = 'Data/';
subj_name = {'subj03','subj04','subj05','subj06','subj07','subj08','subj09','subj10','subj11','subj12','subj13','subj14','subj15','subj16','subj17','subj18'};
data = loadDataR(folder,subj_name,delayStep);
toc
% save data data
%%
% load data
colors = [0 0 0
          128 128 128
          255 69 0
          255 165 0]./255;
Nsubj = length(subj_name);

hands = {'average','bimanual'};
mental = {'no_rotation','rotation'};
corsi = {'no_corsi','corsi'};

%
% t = {'avg_nC_nR','avg_nC_R','avg_C_nR','avg_C_R','bi_nC_nR','bi_nC_R','bi_C_nR','bi_C_R'};
% Mental Rotation (MR), Spatial Working Memory (SWM)
t_name = {'avg Baseline','avg corsi','avg MR','avg corsi+MR','bi Baseline','bi corsi','bi MR','bi corsi+MR'};

%%
subj = 1; % subject to be analyzed

% probability density of reach direction error from one participant
figure(1); clf; hold on
title(subj_name{subj})
col = lines;
pts = -180:.1:180;
for i = 1:2
    for j = 1:2
        for k = 1:2
            s_error = data{subj}.(hands{i}).(mental{j}).(corsi{k}).error * 180 / pi;
            fks = ksdensity(s_error, pts,'Bandwidth',10);
            if i == 1
                plot(pts, fks, 'LineWidth', 2, 'color', col(2*(j-1)+k,:));
            else
                plot(pts, fks, '--', 'LineWidth', 2, 'color', col(2*(j-1)+k,:));
            end
        end
    end
end
plot(-90*[1,1],[0,.04]);
plot(90*[1,1],[0,.04]);
plot(0*[1,1],[0,.04]);
xlabel('Reach direction error')
xticks(-180:90:180)
axis([-180 180 0 .05])
legend(t_name)
set(gca,'TickDir','out')
hold off

idx = 1;
% reach direction vs target direction from one participant
figure(2); clf; hold on
for i = 1:2
    for j = 1:2
        for k = 1:2
            target_angle = data{subj}.(hands{i}).(mental{j}).(corsi{k}).targetAngle * 180/pi;
            reach_angle = data{subj}.(hands{i}).(mental{j}).(corsi{k}).reachAngle * 180/pi;
            subplot(2,4,idx); hold on
            plot([-185 185], [-185 185], 'k')
            plot(target_angle,reach_angle,'.r','MarkerSize',10)
            xlabel('Target direction')
            ylabel('Initial cursor direction')
            daspect([1 1 1])
            axis([-185 185 -185 185])
            xticks(-180:90:180)
            yticks(-180:90:180)
            title(t_name{idx})
            set(gca,'TickDir','out')
            box off
            
            idx = idx + 1;
        end
    end
end
hold off
%% New figure 6. Reaching direction, path length, movement duration

clear error_all pathLength_all reachTime_all
for i = 1:2
    for j = 1:2
        for k = 1:2
            % collecting the data from all subj into dataset
            for m = 1:Nsubj
                error_all.(hands{i}).(mental{j}).(corsi{k})(m) = mean(data{m}.(hands{i}).(mental{j}).(corsi{k}).error) * 180/pi;
                pathLength_all.(hands{i}).(mental{j}).(corsi{k})(m) = mean(data{m}.(hands{i}).(mental{j}).(corsi{k}).pathLength);
                reachTime_all.(hands{i}).(mental{j}).(corsi{k})(m) = mean(data{m}.(hands{i}).(mental{j}).(corsi{k}).timeReach) / 1000;
                sd_all.(hands{i}).(mental{j}).(corsi{k})(m) = data{m}.(hands{i}).(mental{j}).(corsi{k}).sd;
                unif_wt_all.(hands{i}).(mental{j}).(corsi{k})(m) = data{m}.(hands{i}).(mental{j}).(corsi{k}).unif_weight;
            end
        end
    end
end

% plot out
figure(3); clf
subplot(1,3,1); hold on
idx = 1;
for j = 1:2
    for k = 1:2
        avg = error_all.average.(mental{j}).(corsi{k});
        bim = error_all.bimanual.(mental{j}).(corsi{k});
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
ylabel(['Reach direction error (' char(0176) ')'])

subplot(1,3,2); hold on
idx = 1;
for j = 1:2
    for k = 1:2
        avg = pathLength_all.average.(mental{j}).(corsi{k});
        bim = pathLength_all.bimanual.(mental{j}).(corsi{k});
        pathLength = [avg; bim];
        
        color = col((j-1)*2+k,:);
        plot(idx:idx+1, pathLength, 'Color', [color 0.5], 'HandleVisibility', 'off')
        plot(idx:idx+1, mean(pathLength,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
        
        idx = idx + 2;
    end
end
xticks(1:8)
xticklabels({'Average','Bimanual','Average','Bimanual','Average','Bimanual','Average','Bimanual'})
xtickangle(45)
ylabel('Path length (m)')
legend({'No dual task','Corsi','MR','Corsi+MR'})

subplot(1,3,3); hold on
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

figure(4); clf
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

figure(5); clf
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


%%

figure(6); clf
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
legend({'avg','avg Corsi','bim','bim Corsi'})

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

%%

figure(7); clf
subplot(1,2,1); hold on
idx = 1;
for i = 1:2
    for j = 1:2
        for k = 1:2
            sd = sd_all.(hands{i}).(mental{j}).(corsi{k});
            
            color = col((j-1)*2+k,:);
            plot(idx + (rand(Nsubj,1)-0.5) * 0.5, sd, '.', 'Color', [color 0.5], 'MarkerSize', 20, 'HandleVisibility', 'off')
            plot(idx, mean(sd,2),'-o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
            
            idx = idx + 1;
        end
    end
end
xticks(1:8)
xticklabels(t_name)
xtickangle(45)
ylabel(['Circular standard deviation (' char(0176) ')'])

subplot(1,2,2); hold on
idx = 1;
for i = 1:2
    for j = 1:2
        for k = 1:2
            unif_wt = unif_wt_all.(hands{i}).(mental{j}).(corsi{k});
            
            color = col((j-1)*2+k,:);
            plot(idx + (rand(Nsubj,1)-0.5) * 0.5, unif_wt, '.', 'Color', [color 0.5], 'MarkerSize', 20, 'HandleVisibility', 'off')
            plot(idx, mean(unif_wt,2),'o', 'Color', color, 'MarkerFaceColor', color,'MarkerSize',8,'LineWidth',3)
            
            idx = idx + 1;
        end
    end
end
xticks(1:8)
xticklabels(t_name)
xtickangle(45)
ylabel('Proportion of random reaches')

%%
function data = loadDataR(folder, subj_name, ~)

disp('Analyzing...')
task = {'avg-nC-nR','avg-nC-R','avg-C-nR','avg-C-R','bi-nC-nR','bi-nC-R','bi-C-nR','bi-C-R'};
center = [.6 .25];

for l = 1:length(subj_name) % loop over subjects
    disp(['    ',subj_name{l}])
    for k = 1:8 % loop over tasks
        
        if k <= 4
            hands = 'average';
        else
            hands = 'bimanual';
        end
        
        if mod(k,2) == 0
            mental = 'rotation';
        else
            mental = 'no_rotation';
        end
        
        if k == 1 || k == 2 || k == 5 || k == 6
            corsi = 'no_corsi';
        else
            corsi = 'corsi';
        end
        
        % initialized the direction errors
        error = [];
        r_angle = [];
        t_angle = [];
        pathLength = [];
        timeReach = [];
        path = [folder subj_name{l} '/' task{k} '/']; % set data path
        fnames = dir(path); % get filenames in path
        fnames = fnames(not([fnames.isdir])); % get rid of directories in fnames
        
        % import target order file
        opts = detectImportOptions([path 'tOrder.txt'], 'FileType', 'text');
        opts.Delimiter = {' '};
        opts.VariableNames = {'trial','target1','target2','target3','target4','target5'};
        opts.DataLines = [1 Inf];
        opts = setvartype(opts,{'trial','target1','target2','target3','target4','target5'},'double');
        T = readtable([path 'tOrder.txt'],opts);
     
        % get rid of trial number information
        tOrder = T{:,2:6};
       
        Ntrials = size(fnames,1)-1; % number of trials
        for j = 1:Ntrials % loop over trials
            d = dlmread([path fnames(j).name],' ',11,0); % read data into d
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % identify correct answers
            answer = tOrder(j,:); % extract the correct answer for trial j
            answer = answer(~isnan(answer)); % get rid of 0's when sequence length<5
            Ntargets = length(answer); % number of targets in a trial
            
            time = d(:,8)-d(1,8); % time vector
            trajRaw = d(:,5:6); % raw trajectory
            trajFilt = sgolayfilt(trajRaw,3,7); % savitzky-golay filtered trajectories
%             vel = diff(trajFilt)./diff(time/1000); % compute velocity of movements
            state = d(:,7);
            start = [];
            stop = [];
            for i = 1:length(state)
                if (state(i) == 4) && (state(i-1) ~= 4) % find times of reach start
                    start = [start i];
                elseif (state(i) == 4) && (state(i+1) ~= 4) % find times of reach end
                    stop = [stop i];
                end
            end
            
            traj = {};
            
            for i = 1:Ntargets
                
                % pull out trial information between reach start and end
                trajRaw_window = trajRaw(start(i):stop(i),:);
                time_window = time(start(i):stop(i),:);
                
                r = 0.12; % distance to targets
                min_diff = r; % smallest distance between cursor and outer circle
                reach = trajRaw_window(1,:); % current cursor position
                reach_index = 1;
                for ti = 1:length(trajRaw_window)
                    cdistance = sqrt(sum((trajRaw_window(ti,:)-center).^2)); % distance between cursor and center
                    ti_diff = abs(cdistance - r); % distance to outer circle
                    
                    % if distance to outer circle is smaller than smallest
                    % distance...
                    if ti_diff < min_diff
                        reach = trajRaw_window(ti,:);
                        min_diff = ti_diff; % record new smallest distance
                        reach_index = ti;
                    end
                end
                
%                 if abs(reach_index - (stop(i) - start(i))) > 1
%                     1;
%                 end
                trajFilt_p = sgolayfilt(trajRaw_window(1:reach_index,:),3,7); % savitzky-golay filtered trajectories
                dpath = diff(trajFilt_p,1); 
                dL = sqrt(sum(dpath.^2,2)); % difference in cursor position at each time
                pathLength_p = sum(dL); % total path length 
                reach_time = time_window(reach_index)-time_window(1);
                target_angle = answer(i); % angle of target
                reach_vector = reach - center;
                reach_angle = atan2(reach_vector(2),reach_vector(1)); % angle of reach
                
                % adjust target angle for trials with mental rotation
                if k == 2 || k ==4 || k==6 || k==8
                    target_angle = target_angle - pi/2;
                end
                
                % unwrap target and reach angles
                while target_angle > pi
                    target_angle = target_angle - 2*pi;
                end
                while target_angle <= -pi
                    target_angle = target_angle + 2*pi;
                end
                while reach_angle > pi
                    reach_angle = reach_angle - 2*pi;
                end
                while reach_angle <= -pi
                    reach_angle = reach_angle + 2*pi;
                end
                
                % compute reach direction error
                error_pi = target_angle - reach_angle;
                while error_pi > pi
                    error_pi = error_pi - 2*pi;
                end
                while error_pi <= -pi
                    error_pi = error_pi + 2*pi;
                end
                
                error = [error error_pi];
                r_angle = [r_angle reach_angle];
                t_angle = [t_angle target_angle];
                pathLength = [pathLength pathLength_p];
                timeReach = [timeReach reach_time];
                traj{i} = trajFilt_p;
                
            end
            % store variables from each trial in data
            data{l}.(hands).(mental).(corsi).trajFilt{j} = traj;
            data{l}.(hands).(mental).(corsi).answer{j} = answer;
        end
        
        samples = error;
        if isrow(samples)
            samples = samples';
        end
        
        % set values for initial parameters for optimization
        muInit = 0; % mean of von Mises distribution
        kappaInit = 1; % concentration parameter of von Mises
        weightInit = 0.9; % relative weight of von Mises and uniform distributions

        % fit model
        log_likelihood = @(params) calc_likelihood(params, samples);
        paramsInit = [muInit kappaInit weightInit]; % set parameters to current values of mu and kappa
        params_opt = fmincon(log_likelihood, paramsInit, [], [], [], [], [-pi 0 0], [pi 500 1]);
        
        % store fitted parameter values
%         mu_opt = params_opt(1);
        kappa_opt = params_opt(2);
        weight_opt = params_opt(3);

        % compute circular standard deviation
        R = besseli(1,kappa_opt)/besseli(0,kappa_opt);
        sd = sqrt(-2 * log(R)) * 180 / pi; % circular standard deviation
        unif_weight = 1-weight_opt;
       
        % save variables from across trials in data
        data{l}.(hands).(mental).(corsi).error = error; 
        data{l}.(hands).(mental).(corsi).sd = sd;
        data{l}.(hands).(mental).(corsi).unif_weight = unif_weight;
        data{l}.(hands).(mental).(corsi).reachAngle = r_angle;
        data{l}.(hands).(mental).(corsi).targetAngle = t_angle;
        data{l}.(hands).(mental).(corsi).pathLength = pathLength;
        data{l}.(hands).(mental).(corsi).timeReach = timeReach;
    end
end
end

% function for computing log-likelihood
function neg_log_likelihood = calc_likelihood(params,samples)
    mu = params(1);
    kappa = params(2);
    weight = params(3);
    
    likelihood_unif = ones(size(samples)) .* ((1 - weight) / (2*pi));
    likelihood_vm = weight * exp(kappa * cos(samples-mu)) / (2 * pi * besseli(0,kappa));
    
    likelihood_all = sum([likelihood_unif likelihood_vm],2);
    neg_log_likelihood = -sum(log(likelihood_all));
end
