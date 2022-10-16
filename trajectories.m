%% plot trajectories rotated into common reference frame

subj = 2;

hands = {'average','bimanual'};
mental = {'no_rotation','rotation'};
corsi = {'no_corsi','corsi'};
idx = 1;

angles = -pi:pi/32:pi;

figure(1); clf
for i = 1:length(hands)
    for j = 1:length(mental)
        for k = 1:length(corsi)
            subplot(2,4,idx); hold on
            
            block = data{subj}.(hands{i}).(mental{j}).(corsi{k});
            traj = block.trajRot;
            Ntrials = length(traj);
            idx2 = 1;
            
            plot(0.12*cos(angles), 0.12*sin(angles), 'k')
            for m = 1:Ntrials
                Ntargets = length(traj{m});
                
                for n = 1:Ntargets
                    time = block.initTime(idx2);
                    
                    plot(traj{m}{n}(:,1), traj{m}{n}(:,2), 'Color', [0 0 0 0.2])
%                     plot(traj{m}{n}(time,1), traj{m}{n}(time,2), '.g', 'MarkerSize', 10)
                    
                    initDir = block.initDir(idx2) - block.targetAngle(idx2) + pi/2;
                    
%                     x = 0.01*cos(initDir);
%                     y = 0.01*sin(initDir);
%                     
%                     x2 = [traj{m}{n}(time,1) traj{m}{n}(time,1)+x];
%                     y2 = [traj{m}{n}(time,2) traj{m}{n}(time,2)+y];
%                     plot(x2,y2,'g')
                    
                    idx2 = idx2 + 1;
                end
            end

            plot(0, 0.12, 'ok', 'MarkerFaceColor', 'y', 'MarkerSize', 5)

            if idx == 1
                title('No concurrent task')
                ylabel('Average')
            elseif idx == 2
                title('Corsi')
            elseif idx == 3
                title('MR')
            elseif idx == 4
                title('Corsi + MR')
            elseif idx == 5
                ylabel('Bimanual')
            end
            set(gca, 'Tickdir', 'out')
            axis([-0.15 0.15 -0.15 0.15])
            axis square
            idx = idx + 1;
        end
    end
end

print('C:\Users\Chris\Dropbox\Conferences\SFN 2022\traj2','-dpdf','-painters')
%% plot unrotated trajectories

idx = 1;
figure(2); clf
for i = 1:length(hands)
    for j = 1:length(mental)
        for k = 1:length(corsi)
            subplot(2,4,idx); hold on
            
            block = data{subj}.(hands{i}).(mental{j}).(corsi{k});
            traj = block.trajFilt;
            targ = 0.12 * [cos(block.targetAngle) sin(block.targetAngle)];
            Ntrials = length(traj);
            idx2 = 1;
            
            plot(0.12*cos(angles), 0.12*sin(angles), 'k')
            for m = 1:Ntrials
                Ntargets = length(traj{m});
                
                for n = 1:Ntargets
                    time = block.initTime(idx2);
                    trial = [traj{m}{n}(:,1) - 0.6 traj{m}{n}(:,2) - 0.25];
                    
                    plot(trial(:,1), trial(:,2), 'Color', [0 0 0 0.3])
                    plot(trial(time,1), trial(time,2), '.g', 'MarkerSize', 10)
                    
                    initDir = block.initDir(idx2);
                    
                    x = 0.01*cos(initDir);
                    y = 0.01*sin(initDir);
                    
                    x2 = [trial(time,1) trial(time,1)+x];
                    y2 = [trial(time,2) trial(time,2)+y];
                    plot(x2,y2,'g')
                    
                    idx2 = idx2 + 1;
                end
            end
            
            if idx == 1
                title('No concurrent task')
                ylabel('Average')
            elseif idx == 2
                title('Corsi')
            elseif idx == 3
                title('MR')
            elseif idx == 4
                title('Corsi + MR')
            elseif idx == 5
                ylabel('Bimanual')
            end
            axis equal
            idx = idx + 1;
        end
    end
end

%%
col = lines;
idx = 1;

subj = 3;
trialNum = 2;

figure(3); clf
for i = 1:length(hands)
    for j = 1:length(mental)
        for k = 1:length(corsi)
            subplot(2,4,idx); hold on
            
            block = data{subj}.(hands{i}).(mental{j}).(corsi{k});
            traj = block.trajFilt{trialNum};
            targ = 0.12 * [cos(block.targetAngle) sin(block.targetAngle)];
            idx2 = (trialNum - 1) * 5;
            
            plot(0.12*cos(angles), 0.12*sin(angles), 'k')
            Ntargets = length(traj);
            
            for n = 1:Ntargets
                
                trial = [traj{n}(:,1) - 0.6 traj{n}(:,2) - 0.25];
                
                plot(trial(:,1), trial(:,2), 'Color', col(n, :), 'LineWidth', 2)
                if k == 1
                    plot(targ(trialNum, 1), targ(trialNum, 2), 'ok', 'MarkerFaceColor', col(n, :), 'MarkerSize', 5)
                else
                    plot(targ(idx2 + n, 1), targ(idx2 + n, 2), 'ok', 'MarkerFaceColor', col(n, :), 'MarkerSize', 5)
                end
            end
            
            if idx == 1
                title('No concurrent task')
                ylabel('Average')
            elseif idx == 2
                title('Corsi')
            elseif idx == 3
                title('MR')
            elseif idx == 4
                title('Corsi + MR')
            elseif idx == 5
                ylabel('Bimanual')
            end
            axis([-0.15 0.15 -0.15 0.15])
            axis square
            idx = idx + 1;
        end
    end
end

print('C:\Users\Chris\Dropbox\Conferences\SFN 2022\traj1','-dpdf','-painters')