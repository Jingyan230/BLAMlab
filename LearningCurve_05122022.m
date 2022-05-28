clear all
close all
tic
delayStep = 20;
space = 6;

% loading data

% The data files were stored under the folder 'Data/subject/subject/'.
% The data of each subject were stored under different folder, the name of
% the folder is listed as 'subj_name'. Under the folder of each subject,
% the data from each task were stored in different folder.

% folder = 'Data/subject/subject/';
% subj_name = {'subj03'};
% data = loadDataB(folder,subj_name,delayStep);
% toc
% save LearningCurve
%%
load LearningCurve
colors = [0 0 0
          128 128 128
          255 69 0
          255 165 0]./255;
Nsubj = length(subj_name);
disp('Done')

bi_task = {'d1_prac_avg1','d1_prac_bi1','d1_prac_bi2','d1_prac_bi3','d1_prac_bi4','d1_prac_bi5','d1_prac_bi6','d2_prac_bi7'};

% Figure 3. the scater plot of reaching angle from one participant
for i = 1
    figure, hold on
    y = [];
    err = [];
    for j = 1:8
        y_nPL = mean(data{1, i}.(bi_task{j}).normalPL);
        err_nPL = std(data{1, i}.(bi_task{j}).normalPL);
        y = [y y_nPL];
        err = [err err_nPL];
    end
    x = 1:8;
    errorbar(x,y,err, 'LineWidth', 2)
    xlabel('Practice phase')
    ylabel('normalized path length')
    xlim([0 9])
    xticks(1:8)
    xticklabels({'d1-avg','d1-bi1','d1-bi2','d1-bi3','d1-bi4','d1-bi5','d1-bi6','d2-bi'})
    legend(subj_name{i});
    title('Learning Curve')
    hold off
end

%%
function data = loadDataB(folder, subj_name, ~)

disp('Analyzing...')
task = {'d1_prac_avg1','d1_prac_bi1','d1_prac_bi2','d1_prac_bi3','d1_prac_bi4','d1_prac_bi5','d1_prac_bi6','d2_prac_bi7'};
bi_data = {'d1_prac_avg1','d1_prac_bi1','d1_prac_bi2','d1_prac_bi3','d1_prac_bi4','d1_prac_bi5','d1_prac_bi6','d2_prac_bi7'};
center = [.6 .25];

for l = 1:length(subj_name) % loop over subjects
    disp(['    ',subj_name{l}])
    for k = 1:8 % loop over tasks
        % initialized the direction errors
        pathLength = [];
        timeReach = [];
        normalPL = [];
        path = [folder subj_name{l} '/' task{k} '/']; % set data path
        fnames = dir(path); % get filenames in path
        fnames = fnames(not([fnames.isdir])); % get rid of directories in fnames
        
        % import target order file
        opts = detectImportOptions([path 'tFile.tgt'], 'FileType', 'text');
        opts.Delimiter = {' '};
        opts.VariableNames = {'trial','tx','ty','time'};
        opts.DataLines = [1 Inf];
        opts = setvartype(opts,{'trial','tx','ty','time'},'double');
        T = readtable([path 'tFile.tgt'],opts);
     
        % get rid of trial number information
        % tOrder = T{:,2:end};
       
        Ntrials = size(fnames,1)-1; % number of trials
        for j = 1:Ntrials % loop over trials
            d = dlmread([path fnames(j).name],' ',8,0); % read data into d
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % identify correct answers
            time = d(:,9)-d(1,9); % time vector
            trajRaw = d(:,5:6); % raw trajectory
            trajFilt = sgolayfilt(trajRaw,3,7); % savitzky-golay filtered trajectories
            vel = diff(trajFilt)./diff(time/1000); % compute velocity of movements
            state = d(:,7);
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
                trajFilt_p = sgolayfilt(trajRaw_p,3,7); % savitzky-golay filtered trajectories
                dpath = diff(trajFilt_p,1);
                dL = sqrt(sum(dpath.^2,2));
                pathLength_p = sum(dL);
                normalizePL = pathLength_p/0.12;
                normalPL = [normalPL normalizePL];
                pathLength = [pathLength pathLength_p];
            end
            % store variables from each trial in data
            data{l}.(bi_data{k}).trajRaw{j} = trajRaw;
            data{l}.(bi_data{k}).trajFilt{j} = trajFilt;
            data{l}.(bi_data{k}).vel{j} = vel;
        end
       
        % save variables from across trials in data
        data{l}.(bi_data{k}).pathLength = pathLength;
        data{l}.(bi_data{k}).normalPL = normalPL;
    end
end
end
