%% fit models for parameter recovery analysis

rng(7);

N_vm = 50;

kappa = 0:10;
weight = 0:.1:1;

R = besseli(1,kappa)./besseli(0,kappa);
sd = sqrt(-2 * log(R)) * 180 / pi;

% set values for initial parameters for optimization
muInit = 0; % mean of von Mises distribution
kappaInit = 1; % concentration parameter of von Mises
weightInit = 0.9; % relative weight of von Mises and uniform distributions
Ntrials = 50;

tic
mu_opt = NaN(Ntrials, length(kappa));
kappa_opt = NaN(Ntrials, length(kappa));
weight_opt = NaN(Ntrials, length(kappa));
sd_opt = NaN(Ntrials, length(kappa));

for j = 1:length(weight)
    for i = 1:length(kappa)
        for m = 1:Ntrials
            N_vm = round(50 * weight(j));
            
            if N_vm == 0
                samples_vm = [];
            else
                samples_vm = vmrand(0,kappa(i),[N_vm 1]);
            end
            
            if N_vm == 1
                samples_unif = [];
            else
                samples_unif = (rand([50-N_vm 1]) - 0.5) * 2 * pi;
            end
            
            samples = [samples_vm; samples_unif];
            
            % fit model
            log_likelihood = @(params) calc_likelihood(params, samples);
            paramsInit = [muInit kappaInit weightInit]; % set parameters to current values of mu and kappa
            params_opt = fmincon(log_likelihood, paramsInit, [], [], [], [], [-pi 0 0], [pi 500 1]);
            
            % store fitted parameter values
            mu_opt(m,i,j) = params_opt(1) * 180/pi;
            kappa_opt(m,i,j) = params_opt(2);
            weight_opt(m,i,j) = params_opt(3);
            
            % compute circular standard deviation
            R = besseli(1,params_opt(2))/besseli(0,params_opt(2));
            sd_opt(m,i,j) = sqrt(-2 * log(R)) * 180 / pi; % circular standard deviation
        end
    end
end
toc

%% recovery of kappa for fixed value of weight w

w = 1; % index of weight to be used

figure(1); clf; 
subplot(1,3,1); hold on
plot([0 kappa(end)], [0 kappa(end)], 'k');
for i = 1:length(kappa)
    plot(kappa(i) + 0.4*(rand([Ntrials 1]) - 0.5), kappa_opt(:,i,w), '.r', 'MarkerSize', 10)
    plot(kappa(i), mean(kappa_opt(:,i,w)), 'ok', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', 12)
end
ylim([0 25])
xlabel('Actual \kappa')
ylabel('Fitted \kappa')

sd_sd = squeeze(std(sd_opt, [], 1));
subplot(1,3,2); hold on
plot([0 sd(2)], [0 sd(2)], 'k');
for i = 1:length(sd)
    plot(sd(i) + 0.4*(rand([Ntrials 1]) - 0.5), sd_opt(:,i,w), '.r', 'MarkerSize', 10)
    errorbar(sd(i), mean(sd_opt(:,i,w)), sd_sd(i,w), 'ok', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', 12)
end
title(['Weight of von Mises = ' num2str(weight(w))])
xlabel('Actual circ stdev')
ylabel('Fitted circ stdev')

subplot(1,3,3); hold on
plot(kappa,sd_sd(:,w), '-or', 'MarkerFaceColor', 'r')
xlabel('Actual \kappa')
ylabel('Stdev of circ stdev')

%% recovery of weight for fixed value of kappa

k = 1; % index of kappa

figure(2); clf; hold on
plot([0 1], [0 1], 'k')
for i = 1:length(weight)
    plot(weight(i), weight_opt(:,k,i), '.r', 'MarkerSize', 10)
    plot(weight(i), mean(weight_opt(:,k,i)), 'ok', 'MarkerFaceColor', 'r', 'MarkerSize', 12)
end
title(['\kappa = ' num2str(kappa(k)) ])
xlabel('Actual weight')
ylabel('Fitted weight')

%% heatmap of parameter values where kappa (circ st dev) is recoverable

threshold = 15; % threshold below which circ st dev is recoverable

clims = [0 1];
col1 = [1 1 1];
col2 = [1 0 0];
Nstep = 100;
map = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2)...
    ,Nstep)', linspace(col1(3),col2(3),Nstep)'];

sd_sd = squeeze(std(sd_opt, [], 1));
test = sd_sd < threshold;
figure(3); clf; hold on
imagesc(test, clims)
colormap(map)
title(['Std < ' num2str(threshold)])
xlabel('Weight (vm)')
ylabel('\kappa')
xticks(1:length(weight))
yticks(1:length(kappa))
xticklabels(weight)
yticklabels(kappa)
axis([0.5 length(weight)+0.5 0.5 length(kappa)+0.5])
axis square
set(gca,'TickDir','out')
box off

%% kappa and weight of participant data

idx = 1;
figure(4); clf
subplot(1,2,1); hold on
for i = 1:length(hands)
    for j = 1:length(mental)
        for k = 1:length(corsi)
            for m = 1:16
                kappa(m) = data{m}.(hands{i}).(mental{j}).(corsi{k}).kappa;
            end
            plot(idx, kappa, '.k', 'MarkerSize', 10)
            plot(idx, mean(kappa), 'ok', 'MarkerFaceColor', 'r', 'MarkerSize', 8)
            idx = idx+1;
        end
    end
end
xticks(1:8)
xticklabels(t_name)
xtickangle(45)
ylabel('\kappa')
ylim([0 120])

idx = 1;
subplot(1,2,2); hold on
for i = 1:length(hands)
    for j = 1:length(mental)
        for k = 1:length(corsi)
            for m = 1:16
                vm_weight(m) = 1 - data{m}.(hands{i}).(mental{j}).(corsi{k}).unif_weight;
            end
            plot(idx, vm_weight, '.k', 'MarkerSize', 10)
            plot(idx, mean(vm_weight), 'ok', 'MarkerFaceColor', 'r', 'MarkerSize', 8)
            idx = idx+1;
        end
    end
end
xticks(1:8)
xticklabels(t_name)
xtickangle(45)
ylabel('von Mises weight')
ylim([0 1])

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