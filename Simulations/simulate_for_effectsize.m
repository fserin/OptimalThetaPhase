%% Adapted from Zoefel et al. (2019)
% Simulations for calculating an area effect across phase lags that
% estimates the pairwise synchronous vs asynchronous theta contrast with
% the effect size of d = 0.4

% The simulation runs 10000 subjects to make the estimation.
% Area parameters are tested manually in 0.5 increments to see which
% parameter gives d = 0.4

% We chose to use the 0- vs. 180-degree contrast as it parallels
% Clouter et al. the best as well as having a better suited distribution

% the simulations require circstats toolbox
addpath circstat-matlab-master

% to suppress fitglm warning
% should be run after an instance of warning for it to work
% w = warning('query','last')
% id = w.identifier
% warning('off',id)

%% this script requires the circular statistics toolbox
clear;

rng(1996)

%% general parameters
% number of experiments to run
no_tests = 1;

% number of phase bins used to chatacterize the "true" effect
no_phasebins_sampl = 360;
allphases = 1:no_phasebins_sampl;

%% experimental parameters
% 'baseline' in paper
offset = 0.5;

% number of phase bins
no_phasebins = 16

% bin used for alignment
centerbin = no_phasebins/2+1;

% number of trials
no_trials_bin = 192;

% 1 = randomly sampled, 2 = imposed externally
mode = 2;

% number of subjects per experiment
no_subjs = 10000

% 'total width' in the paper. tested values are 1, 0.75, 0.5
w_total = 0.5;

% 'asymmetry' in paper. w_pos + w_neg = 1; example shown corresponds to symmetric effect.
w_pos = 0.5;
w_neg = 0;

% Area Parameter
AP = 5.5
AP_norm = AP/100;

% Determine how to select optimal phase for each subject

% Uncomment this for uniform randomness within -180 to 180
% PhaseSelection = 'UniformRandom'
% Uncomment this for uniform randomness within -45 to 45
% PhaseSelection = 'RestrictedUniform'
% Uncomment this for triangular randomness within -45 to 45
PhaseSelection = 'Triangular'

% Triangular distribution to restrict optimal phase
TriangleDist = makedist('Triangular',-45,0,45);
% histogram(random(TriangleDist,10000000,1))

%% NO EFFECT
nullR_sub = nan(no_trials_bin,no_subjs);
nullC_sub = nan(no_trials_bin,no_subjs);

for nSubs = 1:length(no_subjs)
    for t = 1:no_tests
        for subj = 1:no_subjs(nSubs)
            subj
            % simulate the data
            % input variables: effect present or not (1/0), *number* of phase bin with peak, *number* of phase bin with trough,
            % width of positive peak, width of negative peak, total width, effect size, baseline, number of trials (/192), number of phase bins,
            % number of phase bins to characterize the true effect, study design
            [phasetrials,~]=generate_phase_effect2(0,0,0,0,0,0,0,offset,no_trials_bin,no_phasebins,no_phasebins_sampl,mode);
            
            pred = nan(size(phasetrials,1),2);
            pred(:,1) = sin(phasetrials(:,4));
            pred(:,2) = cos(phasetrials(:,4));
            
            % logistic regression
            g1 = fitglm(pred,phasetrials(:,1),'distribution','binomial','link','logit');
            
            
            nullR_sub(:,subj) = phasetrials(:,1); % responses per sub
            nullC_sub(:,subj) = phasetrials(:,4); % phase lag of each response per sub
        end
    end
end

%% Calculate and plot cohens d for 0 vs 180 phase bins
NullDiff = nan(1,no_subjs);
for iSub = 1:no_subjs
    
    NullDiff(iSub) = mean(nullR_sub(nullC_sub(:,iSub) == 0,iSub))- mean(nullR_sub(nullC_sub(:,iSub) == -pi,iSub));
    
end

NullDiff_d = mean(NullDiff)/std(NullDiff)
figure
histogram(NullDiff)
xline(mean(NullDiff),'--r',['d = ' num2str(NullDiff_d)],'LabelOrientation','horizontal')
xlabel('0  vs. -180 phase lag difference scores (Effect Absent)')


%% EFFECT PRESENT
R_sub = nan(no_trials_bin,no_subjs);
C_sub = nan(no_trials_bin,no_subjs);

% for permutation
no_it = 1000

for nSubs = 1:length(no_subjs)
    for t = 1:no_tests
        for subj = 1:no_subjs(nSubs)
            subj
            % determine "best" phase for each participant
            
            if strcmp('UniformRandom',PhaseSelection)
                % Uncomment these for uniform randomness within -180 to 180
                randvec = allphases(randperm(no_phasebins_sampl));
                whichphase_pos = randvec(1);
            elseif strcmp('RestrictedUniform',PhaseSelection)
                % Uncomment these for uniform randomness within -45 to 45
                randvec = allphases(randperm(no_phasebins_sampl/4)+no_phasebins_sampl*3/8);
                whichphase_pos = randvec(1);
            elseif strcmp('Triangular',PhaseSelection)
                % Uncomment these for uniform randomness within -45 to 45
                whichphase_pos = round(random(TriangleDist,1,1)) + no_phasebins_sampl/2;
            else
                error('Check phase selection method')
            end
            
            % phase with trough is always 180 degrees away from best phase (phase with peak)
            whichphase_neg = whichphase_pos+(no_phasebins_sampl/2);
            if whichphase_neg > no_phasebins_sampl
                whichphase_neg = whichphase_neg-no_phasebins_sampl;
            end
            
            % now with effect present
            [phasetrials,~]=generate_phase_effect2(1,whichphase_pos,whichphase_neg,w_pos,w_neg,w_total,AP_norm,offset,no_trials_bin,no_phasebins,no_phasebins_sampl,mode);
            
            
            R_sub(:,subj) = phasetrials(:,1); % responses per sub
            C_sub(:,subj) = phasetrials(:,4); % phase lag of each response per sub
            
            %% determine beta values for circular predictors and calculate RMS
            pred = nan(size(phasetrials,1),2);
            pred(:,1) = sin(phasetrials(:,4));
            pred(:,2) = cos(phasetrials(:,4));
            
            %% get p-value for individual participants (F-test)
            g1 = fitglm(pred,phasetrials(:,1),'distribution','binomial','link','logit');
            
            % permutate to get
            Trials = phasetrials(:,1);
            permDiffs = nan(1,no_it);
            
            for it = 1:no_it
                inds = randperm(no_trials_bin);
                permTrials = Trials(inds);
                
                permDiffs(it) = mean(permTrials(C_sub(:,subj) == 0))- mean(permTrials(C_sub(:,subj) == -pi));
            end
            permDiff(subj) = mean(permDiffs);
        end
    end
end

%% Calculate and plot cohens d for 0 vs 180 phase bins
Diff = nan(1,no_subjs);
for iSub = 1:no_subjs
    
    Diff(iSub) = mean(R_sub(C_sub(:,iSub) == 0,iSub))- mean(R_sub(C_sub(:,iSub) == -pi,iSub));
    
end

% get cohen's d
Diff_d = mean(Diff-permDiff)/std(Diff-permDiff)

% plot
figure
histogram(Diff-permDiff)
xline(mean(Diff-permDiff),'--r',['d = ' num2str(Diff_d)],'LabelOrientation','horizontal')
xlabel('0  vs. -180 phase lag difference scores (Effect Present)')

save sims_for_effectsize.mat