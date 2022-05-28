clear all
close all
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

% folder = 'Data/subject/subject/';
% subj_name = {'subj03','subj04','subj05','subj06','subj07','subj08','subj09','subj10','subj11','subj12','subj13','subj14','subj15','subj16','subj17','subj18'};
% data = loadDataR(folder,subj_name,delayStep);
% toc
% save ThesisData
%%
load ThesisData
colors = [0 0 0
          128 128 128
          255 69 0
          255 165 0]./255;
Nsubj = length(subj_name);
disp('Done')

%
t = {'avg_nC_nR','avg_nC_R','avg_C_nR','avg_C_R','bi_nC_nR','bi_nC_R','bi_C_nR','bi_C_R'};
% Mental Rotation (MR), Spatial Working Memory (SWM)
t_name = {'avg Baseline','avg MR','avg SWM','avg SWM+MR','bi Baselie','bi MR','bi SWM','bi SWM+MR'};

% Figure 4. the histogram of reaching direction from one participant
for i = 1
    figure, hold on
    title(subj_name{i})
    xlim([-180 180]), ylim([0 .05]); 
    for j = 1:4
        ck = linspace(.5,1,4);
        s_error = data{1, 1}.(t{j}).error * 180 / pi;
        pts = [-180:.1:180];
        fks = ksdensity(s_error, pts,'Bandwidth',10);
        plot(pts, fks, 'LineWidth', 2, 'color',[0 ck(j) 0]);
    end
    for j = 5:8
        ck = linspace(.5,1,4);
        s_error = data{1, 1}.(t{j}).error * 180 / pi;
        pts = [-180:.1:180];
        fks = ksdensity(s_error, pts,'Bandwidth',10);
        plot(pts, fks, 'LineWidth', 2, 'color',[1 ck(j-4) 0]);
    end
    % ylabel('Error of initial direction')
    plot(-90*[1,1],[0,.04]);
    plot(90*[1,1],[0,.04]);
    plot(0*[1,1],[0,.04]);
    xlabel('Reaching direction error')
    xticks([-180, -90, 0, 90, 180])
    legend(t_name);
    hold off
end

% Figure 3. the scater plot of reaching angle from one participant
for i = 1
    figure, hold on
    for j = 1:8
        ideal_angle = data{1, i}.(t{j}).idealAngle;
        reach_angle = data{1, i}.(t{j}).reachAngle;
        subplot(2,4,j), plot(ideal_angle,reach_angle,'.','MarkerSize',5)
        xlabel('Ideal reaching angle')
        ylabel('Actual reaching angle')
        daspect([1 1 1])
        xlim([-4 4]),ylim([-4 4]); 
        legend(subj_name{i});
        title(t_name{j})
    end
    hold off
end
%% New figure 6. Reaching direction, path length, movement duration

ept_sd = {}; 
ept_wt = {}; 
ept_se = {};
ept_mn = {};
for j = 1:8
    error_sum = [];
    pathLength_sum = [];
    reachTime_sum = [];
    ept_task_sd = zeros(1,3);
    ept_task_wt = zeros(1,3);
    % collecting the data from all subj into dataset
    for i = 1:16
        error_subj = mean(data{1, i}.(t{j}).error);
        error_sum = [error_sum error_subj]; % collecting the mean reaching direction error from all subj into dataset
        pathLength_subj = mean(data{1, i}.(t{j}).pathLength);
        pathLength_sum = [pathLength_sum pathLength_subj]; % collecting the mean path length from all subj into dataset
        reachTime_subj = mean(data{1, i}.(t{j}).timeReach);
        reachTime_sum = [reachTime_sum reachTime_subj]; % collecting the mean movement duration from all subj into dataset
    end
    ept = {error_sum, pathLength_sum, reachTime_sum}; % collecting the three dataset into the cell 'ept'
    % calculating the sigma of reaching direction error
    for k = 1:3
        samples = ept{k};
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
        mu_opt = params_opt(1);
        kappa_opt = params_opt(2);
        weight_opt = params_opt(3);

        % compute circular standard deviation
        R = besseli(1,params_opt(2))/besseli(0,params_opt(2));
        sd = sqrt(-2 * log(R)); % circular standard deviation
        
        % save the data of std of the error, pathLength, time into ept_task_sd
        if k ==1
            ept_task_sd(k) = sd * 180 / pi;
            ept_task_wt(k) = 1-weight_opt;
        else
            ept_task_sd(k) = std(samples);
        end
        
    end
    % calculate the mean among all subject for the three parameters
    ept_task_mn = [mean(error_sum)* 180 / pi, mean(pathLength_sum), mean(reachTime_sum)]; 
    % save the std of different tasks
    ept_sd{end+1} = ept_task_sd;
    % calculate and save the se of different tasks
    ept_se{end+1} = ept_task_sd/(sqrt(16)); 
    % save the mean of different tasks
    ept_mn{end+1} = ept_task_mn;
    ept_wt{end+1} = ept_task_wt;
end

% plot out
title_ept = {'error of reaching direction','path length of each reach','movement duration of each reach'};
ylable_ept = {'error (deg)', 'path length (m)', 'movement duration (ms)'};
figure
for ept_i = 1:3
    subplot(3,1,ept_i)
    hold on
    error_n_y = [ept_mn{1}(ept_i) ept_mn{5}(ept_i)];
    error_r_y = [ept_mn{2}(ept_i) ept_mn{6}(ept_i)];
    error_c_y = [ept_mn{3}(ept_i) ept_mn{7}(ept_i)];
    error_cr_y = [ept_mn{4}(ept_i) ept_mn{8}(ept_i)];
    error_n_err = [ept_se{1}(ept_i) ept_se{5}(ept_i)];
    error_r_err = [ept_se{2}(ept_i) ept_se{6}(ept_i)];
    error_c_err = [ept_se{3}(ept_i) ept_se{7}(ept_i)];
    error_cr_err = [ept_se{4}(ept_i) ept_se{8}(ept_i)];
    x1 = [0.97 1.97];x2 = [0.99 1.99];x3 = [1.01 2.01];x4 = [1.03 2.03];
    errorbar(x1,error_n_y,error_n_err, 'LineWidth', 2)
    errorbar(x2,error_r_y,error_r_err, 'LineWidth', 2)
    errorbar(x3,error_c_y,error_c_err, 'LineWidth', 2)
    errorbar(x4,error_cr_y,error_cr_err, 'LineWidth', 2)
    if ept_i == 1
        ylim([-8 8])
    elseif ept_i == 2
        ylim([0 0.2])
    elseif ept_i == 3
        ylim([0 4000])
    end
    xlim([0.8 2.5])
    ylabel(ylable_ept{ept_i})
    title(title_ept{ept_i})
    legend('Baseline','Mental Rotation','Spatial Working Memory','MR+SWM')
    xticks([1 2])
    xticklabels({'average','bimanual'})
    hold off
end

%% 
% figure 5. Alpha and sigma under Average mapping
e_sigma = [];
e_alpha = [];
gR = {};
gC = {};
gB = {};
e_sigma_mn = [];
e_alpha_mn = [];
e_sigma_se = [];
e_alpha_se = [];
for j = 1:8
    error_sigma_sum = [];
    error_alpha_sum = [];
    if j == 1
        gR_sub = 'n';
        gC_sub = 'n';
        gB_sub = 'Average';
    elseif j == 2
        gR_sub = 'MR';
        gC_sub = 'n';
        gB_sub = 'Average';
    elseif j == 3
        gR_sub = 'n';
        gC_sub = 'SWM';
        gB_sub = 'Average';
    elseif j == 4
        gR_sub = 'MR';
        gC_sub = 'SWM';
        gB_sub = 'Average';
    elseif j == 5
        gR_sub = 'n';
        gC_sub = 'n';
        gB_sub = 'Bimanual';
    elseif j == 6
        gR_sub = 'MR';
        gC_sub = 'n';
        gB_sub = 'Bimanual';
    elseif j == 7
        gR_sub = 'n';
        gC_sub = 'SWM';
        gB_sub = 'Bimanual';
    elseif j == 8
        gR_sub = 'MR';
        gC_sub = 'SWM';
        gB_sub = 'Bimanual';
    end
    for i = 1:16
        error_sigma_subj = data{1, i}.(t{j}).sd;
        error_sigma_sum = [error_sigma_sum error_sigma_subj];
        e_sigma = [e_sigma error_sigma_subj];
        error_alpha_subj = data{1, i}.(t{j}).sdwt;
        error_alpha_sum = [error_alpha_sum error_alpha_subj];
        e_alpha = [e_alpha error_alpha_subj];
        gR{end+1} = gR_sub;
        gC{end+1} = gC_sub;
        gB{end+1} = gB_sub;
    end
    e_sigma_mn = [e_sigma_mn mean(error_sigma_sum)];
    e_alpha_mn = [e_alpha_mn mean(error_alpha_sum)];
    e_sigma_se = [e_sigma_se std(error_sigma_sum)/sqrt(15)];
    e_alpha_se = [e_alpha_se std(error_alpha_sum)/sqrt(15)];
end

figure
subplot(1,2,1)
hold on
sigma_n_mn_a = [e_sigma_mn(1) e_sigma_mn(2)];
sigma_c_mn_a = [e_sigma_mn(3) e_sigma_mn(4)];
sigma_n_se_a = [e_sigma_se(1) e_sigma_se(2)];
sigma_c_se_a = [e_sigma_se(3) e_sigma_se(4)];
x = 1:2;
errorbar(x,sigma_n_mn_a,sigma_n_se_a, 'LineWidth', 2)
errorbar(x,sigma_c_mn_a,sigma_c_se_a, 'LineWidth', 2)
xlim([0.5 2.5]),ylim([0 30])
ylabel('sigma')
title('sigma of error on reaching direction (degrees) - Average')
legend('Baseline','Spatial Working Memory')
xticks([1 2])
xticklabels({'Baseline','Mental Rotation'})
hold off
subplot(1,2,2)
hold on
sigma_n_mn_b = [e_sigma_mn(5) e_sigma_mn(6)];
sigma_c_mn_b = [e_sigma_mn(7) e_sigma_mn(8)];
sigma_n_se_b = [e_sigma_se(5) e_sigma_se(6)];
sigma_c_se_b = [e_sigma_se(7) e_sigma_se(8)];
x = 1:2;
errorbar(x,sigma_n_mn_b,sigma_n_se_b, 'LineWidth', 2)
errorbar(x,sigma_c_mn_b,sigma_c_se_b, 'LineWidth', 2)
xlim([0.5 2.5]), ylim([0 30])
ylabel('sigma')
title('sigma of error on reaching direction (degrees) - Bimanual')
legend('Baseline','Spatial Working Memory')
xticks([1 2])
xticklabels({'Baseline','Mental Rotation'})
hold off

figure
subplot(1,2,1)
hold on
alpha_n_mn_a = [e_alpha_mn(1) e_alpha_mn(2)];
alpha_c_mn_a = [e_alpha_mn(3) e_alpha_mn(4)];
alpha_n_se_a = [e_alpha_se(1) e_alpha_se(2)];
alpha_c_se_a = [e_alpha_se(3) e_alpha_se(4)];
x = 1:2;
errorbar(x,alpha_n_mn_a,alpha_n_se_a, 'LineWidth', 2)
errorbar(x,alpha_c_mn_a,alpha_c_se_a, 'LineWidth', 2)
xlim([0.5 2.5]),ylim([0 0.45])
ylabel('alpha')
title('alpha of error on reaching direction - Average')
legend('Baseline','Spatial Working Memory')
xticks([1 2])
xticklabels({'Baseline','Mental Rotation'})
hold off
subplot(1,2,2)
hold on
alpha_n_mn_b = [e_alpha_mn(5) e_alpha_mn(6)];
alpha_c_mn_b = [e_alpha_mn(7) e_alpha_mn(8)];
alpha_n_se_b = [e_alpha_se(5) e_alpha_se(6)];
alpha_c_se_b = [e_alpha_se(7) e_alpha_se(8)];
x = 1:2;
errorbar(x,alpha_n_mn_b,alpha_n_se_b, 'LineWidth', 2)
errorbar(x,alpha_c_mn_b,alpha_c_se_b, 'LineWidth', 2)
xlim([0.5 2.5]), ylim([0 0.45])
ylabel('alpha')
title('alpha of error on reaching direction - Bimanual')
legend('Baseline','Spatial Working Memory')
xticks([1 2])
xticklabels({'Baseline','Mental Rotation'})
hold off
% figure 7. alpha and sigma of all tasks
figure
hold on
sigma_n_mn = [e_sigma_mn(1) e_sigma_mn(5)];
sigma_r_mn = [e_sigma_mn(2) e_sigma_mn(6)];
sigma_c_mn = [e_sigma_mn(3) e_sigma_mn(7)];
sigma_cr_mn = [e_sigma_mn(4) e_sigma_mn(8)];
sigma_n_se = [e_sigma_se(1) e_sigma_se(5)];
sigma_r_se = [e_sigma_se(2) e_sigma_se(6)];
sigma_c_se = [e_sigma_se(3) e_sigma_se(7)];
sigma_cr_se = [e_sigma_se(4) e_sigma_se(8)];
x1 = [0.97 1.97];x2 = [0.99 1.99];x3 = [1.01 2.01];x4 = [1.03 2.03];
errorbar(x1,sigma_n_mn,sigma_n_se, 'LineWidth', 2)
errorbar(x2,sigma_r_mn,sigma_r_se, 'LineWidth', 2)
errorbar(x3,sigma_c_mn,sigma_c_se, 'LineWidth', 2)
errorbar(x4,sigma_cr_mn,sigma_cr_se, 'LineWidth', 2)
xlim([0.5 2.5])
ylabel('sigma')
title('sigma of error on reaching direction (degrees)')
legend('Baseline','MR','SWM','SWM+MR')
xticks([1 2])
xticklabels({'average','bimanual'})
hold off

figure
hold on
alpha_n_mn = [e_alpha_mn(1) e_alpha_mn(5)];
alpha_r_mn = [e_alpha_mn(2) e_alpha_mn(6)];
alpha_c_mn = [e_alpha_mn(3) e_alpha_mn(7)];
alpha_cr_mn = [e_alpha_mn(4) e_alpha_mn(8)];
alpha_n_se = [e_alpha_se(1) e_alpha_se(5)];
alpha_r_se = [e_alpha_se(2) e_alpha_se(6)];
alpha_c_se = [e_alpha_se(3) e_alpha_se(7)];
alpha_cr_se = [e_alpha_se(4) e_alpha_se(8)];
x1 = [0.97 1.97];x2 = [0.99 1.99];x3 = [1.01 2.01];x4 = [1.03 2.03];
errorbar(x1,alpha_n_mn,alpha_n_se, 'LineWidth', 2)
errorbar(x2,alpha_r_mn,alpha_r_se, 'LineWidth', 2)
errorbar(x3,alpha_c_mn,alpha_c_se, 'LineWidth', 2)
errorbar(x4,alpha_cr_mn,alpha_cr_se, 'LineWidth', 2)
xlim([0.5 2.5])
ylabel('alpha')
title('alpha of error on reaching direction')
legend('Baseline','MR','SWM','SWM+MR')
xticks([1 2])
xticklabels({'average','bimanual'})
hold off


%%
function data = loadDataR(folder, subj_name, ~)

disp('Analyzing...')
task = {'avg-nC-nR','avg-nC-R','avg-C-nR','avg-C-R','bi-nC-nR','bi-nC-R','bi-C-nR','bi-C-R'};
task_data = {'avg_nC_nR','avg_nC_R','avg_C_nR','avg_C_R','bi_nC_nR','bi_nC_R','bi_C_nR','bi_C_R'};
center = [.6 .25];

for l = 1:length(subj_name) % loop over subjects
    disp(['    ',subj_name{l}])
    for k = 1:8 % loop over tasks
        % initialized the direction errors
        error = [];
        r_angle = [];
        i_angle = [];
        pathLength = [];
        timeReach = [];
        path = [folder subj_name{l} '/' task{k} '/']; % set data path
        fnames = dir(path); % get filenames in path
        fnames = fnames(not([fnames.isdir])); % get rid of directories in fnames
        
        % import target order file
        opts = detectImportOptions([path 'tOrder.txt'], 'FileType', 'text');
        opts.Delimiter = {' '};
        opts.VariableNames = {'trial','a1','a2','a3','a4','a5'};
        opts.DataLines = [1 Inf];
        opts = setvartype(opts,{'trial','a1','a2','a3','a4','a5'},'double');
        T = readtable([path 'tOrder.txt'],opts);
     
        % get rid of trial number information
        tOrder = T{:,2:end};
       
        Ntrials = size(fnames,1)-1; % number of trials
        for j = 1:Ntrials % loop over trials
            d = dlmread([path fnames(j).name],' ',11,0); % read data into d
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % identify correct answers
            answer = tOrder(j,:); % extract the correct answer for trial j
            answer = answer(~isnan(answer)); % get rid of 0's when sequence length<5
            time = d(:,8)-d(1,8); % time vector
            trajRaw = d(:,5:6); % raw trajectory
            trajFilt = sgolayfilt(trajRaw,3,7); % savitzky-golay filtered trajectories
            vel = diff(trajFilt)./diff(time/1000); % compute velocity of movements
            state = d(:,7);
            reachPoint = {};
            idealReach = [];
            traj_s = [];
            traj_e = [];
            for si = 1:length(state)
                if (state(si) == 4) && (state(si-1) ~= 4)
                    traj_s = [traj_s si];
                elseif (state(si) == 4) && (state(si+1) ~= 4)
                    traj_e = [traj_e si];
                end
            end
            for tri = 1:length(traj_s)
                trajRaw_p = trajRaw(traj_s(tri):traj_e(tri),:);
                time_p = time(traj_s(tri):traj_e(tri),:);
                r = 0.12;
                min_diff = r;
                reach = trajRaw_p(1,:);
                cdmax = 0;
                reach_index = 1;
                for ti = 1:length(trajRaw_p)
                    cdistance = sqrt(sum((trajRaw_p(ti,:)-center).^2));
                    ti_diff = abs(cdistance - r);
                    if ti_diff < min_diff
                        reach = trajRaw_p(ti,:);
                        min_diff = ti_diff;
                        reach_index = ti;
                    end
                end
                trajFilt_p = sgolayfilt(trajRaw_p(1:reach_index,:),3,7); % savitzky-golay filtered trajectories
                dpath = diff(trajFilt_p,1);
                dL = sqrt(sum(dpath.^2,2));
                pathLength_p = sum(dL);
                reach_time = time_p(reach_index)-time_p(1);
                iangle = answer(tri);
                reach_vector = reach - center;
                rangle = atan2(reach_vector(2),reach_vector(1));
                if k == 2 || k ==4 || k==6 || k==8
                    iangle = iangle - pi/2;
                end
                
                while iangle > pi
                    iangle = iangle - 2*pi;
                end
                while iangle <= -pi
                    iangle = iangle + 2*pi;
                end
                while rangle > pi
                    rangle = rangle - 2*pi;
                end
                while rangle <= -pi
                    rangle = rangle + 2*pi;
                end
                
                error_pi = iangle - rangle;
                while error_pi > pi
                    error_pi = error_pi - 2*pi;
                end
                while error_pi <= -pi
                    error_pi = error_pi + 2*pi;
                end
                % error_pi = error_pi * 180 / pi;
                error = [error error_pi];
                r_angle = [r_angle rangle];
                i_angle = [i_angle iangle];
                reachPoint = {reachPoint trajRaw_p(reach_index,:)};
                idealReach = [idealReach i_angle];
                pathLength = [pathLength pathLength_p];
                timeReach = [timeReach reach_time];
            end
            % store variables from each trial in data
            data{l}.(task_data{k}).trajRaw{j} = trajRaw;
            data{l}.(task_data{k}).trajFilt{j} = trajFilt;
            data{l}.(task_data{k}).vel{j} = vel;
            data{l}.(task_data{k}).answer{j} = answer;
            data{l}.(task_data{k}).reachPoint{j} = reachPoint;
            data{l}.(task_data{k}).idealReach{j} = idealReach;
        end
        
        rng(2);
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
        mu_opt = params_opt(1);
        kappa_opt = params_opt(2);
        weight_opt = params_opt(3);

            % compute circular standard deviation
        R = besseli(1,params_opt(2))/besseli(0,params_opt(2));
        sd = sqrt(-2 * log(R)); % circular standard deviation
        sd_e = sd * 180 / pi;
        wt_sde = 1-weight_opt;
            
            % delt = pi/64; % interval for plotting PDF
            % x = -pi:delt:pi-delt; % vector of points to plot PDF
            % vmPDF = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa))); % PDF of von Mises distribution
            % pdf = vmPDF(x, mu_opt, kappa_opt); % PDF of von Mises distribution
            % mixPDF = weight_opt * pdf + (1-weight_opt) * (1 / (2*pi)); % weight von Mises with uniform distribution

            % figure; clf; hold on
            % histogram(samples*180/pi,x*180/pi,'Normalization','probability');
            % plot(x*180/pi, mixPDF./sum(mixPDF), 'LineWidth', 2)
            % xlim([-180 180])
            % title(task{k},subj_name{l})
       
        % save variables from across trials in data
        data{l}.(task_data{k}).error = error; 
        data{l}.(task_data{k}).variance = var(error * 180 / pi); 
        data{l}.(task_data{k}).sd = sd_e;
        data{l}.(task_data{k}).sdwt = wt_sde;
        data{l}.(task_data{k}).std = std(error * 180 / pi); % standard deviation of error
        data{l}.(task_data{k}).reachAngle = r_angle;
        data{l}.(task_data{k}).idealAngle = i_angle;
        data{l}.(task_data{k}).pathLength = pathLength;
        data{l}.(task_data{k}).timeReach = timeReach;
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
