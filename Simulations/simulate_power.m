%% Adapted from Zoefel et al. (2019)

% This script simulates N experiments (no_tests, defined below) in which a
% phase effect is present and N experiments in which a phase effect is
% absent.

%% this script requires the circular statistics toolbox
clear;

% if the toolbox is is in the current directory
addpath([cd '/circstat-matlab-master'])

rng(1996)

%% general parameters
% number of experiments to run per sample size
no_tests = 10000;

% number of phase bins used to chatacterise the "true" effect
no_phasebins_sampl = 192;

%% experimental parameters
% 'baseline' in Zoefel et al.
offset = 0.5;

% number of phase bins
no_phasebins = 16;

% number of trials
no_trials_bin = 192;

% phase indices
allphases = 1:no_phasebins_sampl;

% 1 = randomly sampled, 2 = imposed externally
mode = 2;

%% neural parameters
% 'total width' of the deflections
w_total = 0.5;

% widths of positive and negative deflections
% 'asymmetry' in Zoefel et al.
w_pos = 0.5;
w_neg = 0;

% Area parameter or effect size
AP = 5.5 ;
AP_norm = AP/100

% number of sample sizes to be tested
% In 16 increments for full counterbalancing
no_subjs = 160:16:320;

% set custom cluster pool
P=cbupool(11); 
P.SubmitArguments = '--ntasks=11 --mem-per-cpu=12G --time=192:00:00';
parpool(P,11)

% preallocate variables
allhs_logregress_fisher_fp = nan(no_tests,length(no_subjs));

%% test false positives (no effect present)    
% for each experiment, these vectors contain 1 if null hypothesis is
% rejected and NaN if otherwise.

parfor nSubs = 1:length(no_subjs)
    no_subjs(nSubs)    

    for t = 1:no_tests
        % matrices to save data needed for the different tests
        % these are combined later for group-level statistics
        logregress_ps_fp = nan(no_subjs(nSubs),1);
        
        for subj = 1:no_subjs(nSubs)
               
            %% simulate the data
            %% input variables: effect present or not (1/0), *number* of phase bin with peak, *number* of phase bin with trough,
            %% width of positive peak, width of negative peak, total width, effect size, baseline, number of trials (/192), number of phase bins,
            %% number of phase bins to characterize the true effect, study design
            [phasetrials, ~]=generate_phase_effect(0,0,0,0,0,0,0,offset,no_trials_bin,no_phasebins,no_phasebins_sampl,mode);

            %% "phasetrials" is a matrix with one line for each trial and the following columns:
            %% 1: single-trial response
            %% 2: phase that corresponds to the phase bin from column 4
            %% 3: number of the phase bin that phase from column 4 falls into
            %% 4: single-trial phase (actual phase out of "no_phasebins_sampl" possible phases)
            
            %% "totalperf" is response averaged across trials and divided into phase bins
            
            %% METHOD: Logistic regression (sine + cosine)
            %% determine beta values for circular predictors
            pred = nan(size(phasetrials,1),2);
            pred(:,1) = sin(phasetrials(:,4));
            pred(:,2) = cos(phasetrials(:,4));
            
            %% get p-value for individual participants (F-test)
            g1 = fitglm(pred,phasetrials(:,1),'distribution','binomial','link','logit');
            logregress_ps_fp(subj) = coefTest(g1);
        end
        
        %% Fisher's method to combine p-values from logistic regression
        comb_p = pval_Fisher(logregress_ps_fp',1);
        allhs_logregress_fisher_fp(t,nSubs) = comb_p <= 0.05;
    end
end


%% Test for true positive - effect present
% preallocate variables
allhs_logregress_fisher = nan(no_tests,1,length(no_subjs));

parfor nSubs = 1:length(no_subjs)
    no_subjs(nSubs)
    
    for t = 1:no_tests      
        %% matrices to save data needed for the different tests
        logregress_ps = nan(no_subjs(nSubs),1);
        
        for subj = 1:no_subjs(nSubs)
            %% determine "best" phase for each participant
            randvec = allphases(randperm(no_phasebins_sampl));
            whichphase_pos = randvec(1);
            
            %% phase with trough is always 180 degrees away from best phase (phase with peak)
            % although negative deflection does not exist for our purposes
            whichphase_neg = whichphase_pos+(no_phasebins_sampl/2);
            if whichphase_neg > no_phasebins_sampl
                whichphase_neg = whichphase_neg-no_phasebins_sampl;
            end
            
            %% now with effect present
            [phasetrials,~]=generate_phase_effect(1,whichphase_pos,whichphase_neg,w_pos,w_neg,w_total,AP_norm,offset,no_trials_bin,no_phasebins,no_phasebins_sampl,mode);
            
            %% METHOD: Logistic regression (sine + cosine)
            %% determine beta values for circular predictors
            pred = nan(size(phasetrials,1),2);
            pred(:,1) = sin(phasetrials(:,4));
            pred(:,2) = cos(phasetrials(:,4));

            %% get p-value for individual participants (F-test)
            g1 = fitglm(pred,phasetrials(:,1),'distribution','binomial','link','logit');
            logregress_ps(subj) = coefTest(g1);            
            
        end
        
        %% Fisher's method to combine p-values from logistic regression
        comb_p = pval_Fisher(logregress_ps',1);
        allhs_logregress_fisher(t,1,nSubs) = comb_p <= 0.05;
    end
end

save all_simulations_binary.mat