clear; close all;

% this is where results are stored
load sims_for_power.mat

for nSubs = 1:length(no_subjs)
    
    FAR_pN(nSubs) = sum(allhs_logregress_fisher_fp(:,nSubs))/no_tests;
    HR_pN(nSubs) = sum(allhs_logregress_fisher(:,:,nSubs))/no_tests;
    
end

%% Plot FPR
figure, hold on;
   
plot(no_subjs,FAR_pN);
ylim([0 .1])
xlim([no_subjs(1)-10 no_subjs(end)+10])
xlabel('Sample Size')
ylabel('False Positive Rate')
yline(.05)
yline(.9)
yline(.8)
legend({'Fishers Method'})

%% Plot Power
figure, hold on;
    
plot(no_subjs,HR_pN);
ylim([0 1])
xlim([no_subjs(1)-10 no_subjs(end)+10])
xlabel('Sample Size')
ylabel('Power')
yticks([0 .5 .8 .9 1])
yline(.05)
yline(.9)
yline(.8)
legend({'Fishers Method'})
