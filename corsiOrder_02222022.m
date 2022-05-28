clear all
rng(8); % random seed
practice = false; % flags whether tFiles are for practice blocks
                 % true: practice tFiles
                 % false: experiment tFiles

% text for file names
if practice
    practiceText = '-practice';
else
    practiceText = '';
end

%% generate tFiles for Corsi
for m = 1:2 % m = 1: no rotation; m = 2: rotation
    
    % vector of sequence lengths for each trial
    sLength = [5*ones(20,1)];
    
    if practice
        Ntrials = 10; % number of trials in block
        Nblocks = 1; % number of blocks to divide trials into
    else
        Ntrials = length(sLength);
        Nblocks = 1;
    end
    
    sLength = datasample(sLength,Ntrials,'Replace',false); % shuffle sequence lengths

    threshold = pi/6; % threshold to limit how close targets can be
    orderAll = cell(1,46);
    
    % loop to generate target locations
    for i = 1:Ntrials

        proceed = 1; % flag to run while loop
        small = 0; % flag for bad randomization
        
        while proceed
            order = 2*pi*(rand(1,sLength(i))-0.5); % randomly generate target directions

            combinations = nchoosek(order,2); % generate pairwise combinations of targets

            % check to see if any consecutive targets have angle difference of
            % less than threshold
            for j = 1:size(combinations,1)
                ang = combinations(j,1) - combinations(j,2); % angle difference
                
                % unwrap angles into [-pi, pi)
                while ang > pi
                    ang = ang - 2*pi;
                end
                while ang <= -pi
                    ang = ang + 2*pi;
                end
                
                if abs(ang) < threshold
                    small = 1; % indicate bad randomization
                end
            end

            % if no bad trials found, then end while loop
            if small == 0
                proceed = 0;

            % else, generate another set of targets
            else
                small = 0;
            end
        end

        orderAll{i} = order; % store all trials
    end
    
    % split trials based on idx
    if practice
        idx{1} = 1:Ntrials;
    else
        idx{1} = 1:Ntrials/2;
        idx{2} = Ntrials/2+1:Ntrials;
    end
    
    for k = 1:Nblocks
        
        if m == 1 
            tFile = fopen(['C-nR' num2str(k) practiceText '.txt'],'w');
        else
            tFile = fopen(['C-R' num2str(k) practiceText '.txt'],'w');
        end
        
        count = 1; % trial counter
        for i = idx{k}
            
            % store target angles as floats
            if i == length(idx{k}) % print for last trial in file
                fprintf(tFile,repmat('%1.4f ', [1 length(orderAll{i})+1]), [count orderAll{i}]);
            else % print for every other trial
                fprintf(tFile,[repmat('%1.4f ', [1 length(orderAll{i})+1]) '\r\n'], [count orderAll{i}]);
            end
            count = count + 1;
        end
        fclose(tFile);
    end
end

%% generate tFiles for no Corsi

for m = 1:2 % m = 1: no rotation; m = 2: rotation
    if practice
        Ntrials = 20;
    else
        Ntrials = 40;
    end
    
    % randomize target locations
    order = 2*pi*(rand(1,Ntrials)-0.5);
    
    if m == 1
        tFile = fopen(['nC-nR' practiceText '.txt'],'w');
    else
        tFile = fopen(['nC-R' practiceText '.txt'],'w');
    end
        
    count = 1; % trial counter
    for i = 1:Ntrials
        % store target angles as floats
        if i == Ntrials
            fprintf(tFile,repmat('%1.4f ', [1 length(order(i))+1]), [i order(i)]);
        else
            fprintf(tFile,[repmat('%1.4f ', [1 length(order(i))+1]) '\r\n'], [i order(i)]);
        end
        count = count + 1;
    end
    fclose(tFile);
end