
subj = 1;

kappa = 0:5:100;
weight = 0:0.01:0.99;

figure(1); clf
idx = 1;

for a = 1:2
    for b = 1:2
        for c = 1:2
            dat = data{subj}.(hands{a}).(mental{b}).(corsi{c});
            samples = dat.initDir_error;
            
            mu_opt = dat.mu;
            kappa_opt = dat.kappa;
            weight_opt = 1 - dat.unif_weight;
            
            log_likelihood = NaN(length(weight), length(kappa));
            
            for j = 1:length(kappa)
                for i = 1:length(weight)
                    log_likelihood(i,j) = calc_likelihood(samples, mu_opt, kappa(j), weight(i));
                end
            end
            
            [X, Y] = meshgrid(weight, kappa);
            subplot(2,4,idx)
            surf(X', Y', log_likelihood); hold on
            plot3(weight_opt, kappa_opt, calc_likelihood(samples, mu_opt, kappa_opt, weight_opt), '.r', 'MarkerSize', 20);
            title(t_name{idx})
            xlabel('weight')
            ylabel('kappa')
            zlabel('neg log likelihood')
            shading interp
            axis square
            
            idx = idx + 1;
        end
    end
end

function neg_log_likelihood = calc_likelihood(samples, mu, kappa, weight)

    dat = samples(~isnan(samples));
    if isrow(dat)
        dat = dat';
    end
    
    likelihood_unif = ones(size(dat)) .* ((1 - weight) / (2*pi));
    likelihood_vm = weight * exp(kappa * cos(dat-mu)) / (2 * pi * besseli(0,kappa));
    
    likelihood_all = sum([likelihood_unif likelihood_vm],2);
    neg_log_likelihood = -sum(log(likelihood_all));
end