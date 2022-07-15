clear all
close all
rng(2);
tic

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
data = loadDataR(folder,subj_name);
toc

% save data data
% load data

col = lines;
Nsubj = length(subj_name);

hands = {'average','bimanual'};
mental = {'no_rotation','rotation'};
corsi = {'no_corsi','corsi'};

t_name = {'avg no con task','avg corsi','avg MR','avg corsi+MR','bi no con task','bi corsi','bi MR','bi corsi+MR'};

function data = loadDataR(folder, subj_name)

disp('Analyzing...')
task = {'avg-nC-nR','avg-nC-R','avg-C-nR','avg-C-R','bi-nC-nR','bi-nC-R','bi-C-nR','bi-C-R'};
center = [.6 .25];
r = 0.12; % distance to targets

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
            Nreaches = 40;
        else
            corsi = 'corsi';
            Nreaches = 50;
        end
        
        % initialized the direction errors
        endpoint_error = NaN(Nreaches,1);
        initDir_error = NaN(Nreaches,1);
        opposite_error = NaN(Nreaches,1);
        min_error = NaN(Nreaches,1);
        end_angle = NaN(Nreaches,1);
        target_angle = NaN(Nreaches,1);
        pathLength = NaN(Nreaches,1);
        timeReach = NaN(Nreaches,1);
        initDir = NaN(Nreaches,1);
        initTime = NaN(Nreaches,1);
        idx = 1;

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
            trajR = {};
            first_trial = idx;
            
            for i = 1:Ntargets
                
                % pull out trial information between reach start and end
                trajRaw_window = trajRaw(start(i):stop(i),:);
                time_window = time(start(i):stop(i),:);
                
                cdistance = sqrt(sum((trajRaw_window-center).^2,2));
                ti_diff = abs(cdistance - r);
                [~, reach_index] = min(ti_diff, [], 1);
                reach = trajRaw_window(reach_index,:);
                
                leave_start = find(cdistance > 0.01,1);
                init_index = leave_start + 20;
                
                trajFilt_p = sgolayfilt(trajRaw_window(1:reach_index,:),3,7); % savitzky-golay filtered trajectories
                dpath = diff(trajFilt_p,1); 
                dL = sqrt(sum(dpath.^2,2)); % difference in cursor position at each time
                pathLength(idx) = sum(dL); % total path length 
                timeReach(idx) = time_window(reach_index)-time_window(1);
                target_angle(idx) = answer(i); % angle of target
                reach_vector = reach - center;
                end_angle(idx) = atan2(reach_vector(2),reach_vector(1)); % angle of reach
                
                % trying 150 ms after movement initiation
                vel = diff(trajFilt_p,1)./(diff(time_window(1:reach_index))/1000);
                if init_index > length(vel)
                    initDir(idx) = NaN;
                    initTime(idx) = NaN;
                else
                    initDir(idx) = atan2(vel(init_index,2), vel(init_index,1));
                    initTime(idx) = init_index;
                end
                
                % rotate cursor trajectory into common reference frame
                if mod(k,2) == 0 % add 90 degrees if there was a mental rotation
                    rotAngle = target_angle(idx) - pi;
                else
                    rotAngle = target_angle(idx) - pi / 2;
                end
                
                rotMat = [cos(rotAngle) sin(rotAngle); -sin(rotAngle) cos(rotAngle)];
                trajRot = (rotMat*(trajFilt_p' - repmat(center', [1 length(trajFilt_p)])))';
                
                % adjust target angle for trials with mental rotation
                if strcmp(mental,'rotation')
                    target_angle(idx) = target_angle(idx) - pi/2;
                    opposite_angle = target_angle(idx) + pi;
                    opposite_error(idx) = angdiff(opposite_angle, end_angle(idx));
                end
                
                % unwrap target and reach angles
                while target_angle(idx) <= -pi
                    target_angle(idx) = target_angle(idx) + 2*pi;
                end
                
                endpoint_error(idx) = angdiff(target_angle(idx), end_angle(idx));
                initDir_error(idx) = angdiff(target_angle(idx), initDir(idx));
                traj{i} = trajFilt_p;
                trajR{i} = trajRot;
                
                idx = idx + 1;
            end
            
            if strcmp(corsi, 'corsi')
                
                [A, B] = meshgrid(target_angle(first_trial:idx-1), end_angle(first_trial:idx-1));
                c=cat(2,A',B');
                pairs = reshape(c,[],2);
                error_pairs = angdiff(pairs(:,1), pairs(:,2));
                
                for i = 1:Ntargets
                    [~, min_idx] = min(abs(error_pairs(Ntargets*(i-1)+1:Ntargets*(i-1)+5)));
                    min_error(first_trial+i-1) = error_pairs(Ntargets*(i-1) + min_idx);
                end
            end
            
            % store variables from each trial in data
            data{l}.(hands).(mental).(corsi).trajFilt{j} = traj;
            data{l}.(hands).(mental).(corsi).trajRot{j} = trajR;
            data{l}.(hands).(mental).(corsi).answer{j} = answer;
        end
        
        samples = initDir_error;
        
        % set values for initial parameters for optimization
        muInit = 0; % mean of von Mises distribution
        kappaInit = 1; % concentration parameter of von Mises
        weightInit = 0.9; % relative weight of von Mises and uniform distributions

        % fit model
        log_likelihood = @(params) calc_likelihood(params, samples);
        paramsInit = [muInit kappaInit weightInit]; % set parameters to current values of mu and kappa
        [params_opt,~,~,~,~,~,hessian] = fmincon(log_likelihood, paramsInit, [], [], [], [], [-pi 0 0], [pi 500 1]);
        
        % store fitted parameter values
        mu_opt = params_opt(1);
        kappa_opt = params_opt(2);
        weight_opt = params_opt(3);

        % compute circular standard deviation
        R = besseli(1,kappa_opt)/besseli(0,kappa_opt);
        sd = sqrt(-2 * log(R)) * 180 / pi; % circular standard deviation
        unif_weight = 1-weight_opt;
        
        % save variables from across trials in data
        data{l}.(hands).(mental).(corsi).endpoint_error = endpoint_error;
        data{l}.(hands).(mental).(corsi).initDir_error = initDir_error;
        data{l}.(hands).(mental).(corsi).opposite_error = opposite_error;
        data{l}.(hands).(mental).(corsi).min_error = min_error;
        data{l}.(hands).(mental).(corsi).mu = mu_opt;
        data{l}.(hands).(mental).(corsi).kappa = kappa_opt;
        data{l}.(hands).(mental).(corsi).sd = sd;
        data{l}.(hands).(mental).(corsi).unif_weight = unif_weight;
        data{l}.(hands).(mental).(corsi).hessian = hessian;
        data{l}.(hands).(mental).(corsi).endAngle = end_angle;
        data{l}.(hands).(mental).(corsi).targetAngle = target_angle;
        data{l}.(hands).(mental).(corsi).pathLength = pathLength;
        data{l}.(hands).(mental).(corsi).timeReach = timeReach;
        data{l}.(hands).(mental).(corsi).initDir = initDir;
        data{l}.(hands).(mental).(corsi).initTime = initTime;
    end
end
end

% function for computing log-likelihood
function neg_log_likelihood = calc_likelihood(params,samples)

    dat = samples(~isnan(samples));
    if isrow(dat)
        dat = dat';
    end

    mu = params(1);
    kappa = params(2);
    weight = params(3);
    
    likelihood_unif = ones(size(dat)) .* ((1 - weight) / (2*pi));
    likelihood_vm = weight * exp(kappa * cos(dat-mu)) / (2 * pi * besseli(0,kappa));
    
    likelihood_all = sum([likelihood_unif likelihood_vm],2);
    neg_log_likelihood = -sum(log(likelihood_all));
end
