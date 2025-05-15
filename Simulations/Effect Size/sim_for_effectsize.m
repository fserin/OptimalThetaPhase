%% methods tested in Zoefel, Davis, Valente, Riecke, 2019

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
no_phasebins = 2;

% bin used for alignment
centerbin = no_phasebins/2+1;

% number of trials
no_trials_bin = 192;

% 1 = randomly sampled, 2 = imposed externally
mode = 2;

%% neural parameters

% 'total width' in the paper. tested values are 1, 0.75, 0.5
w_total = 0.5;

% 'asymmetry' in paper. w_pos + w_neg = 1; example shown corresponds to symmetric effect.
w_pos = 0.5;
w_neg = 0;

% Effect size
ESs = 5;
a_norm = ESs/100;

%% other variables needed below

% number of subjects per experiment
no_subjs = 10000;

%% NO EFFECT
NullAUC = nan(1,no_subjs);
for nSubs = 1:length(no_subjs)
    for t = 1:no_tests
        for subj = 1:no_subjs(nSubs)
            
            % simulate the data
            % input variables: effect present or not (1/0), *number* of phase bin with peak, *number* of phase bin with trough,
            % width of positive peak, width of negative peak, total width, effect size, baseline, number of trials (/192), number of phase bins,
            % number of phase bins to characterize the true effect, study design
            [phasetrials,totalperf, NullAUC(subj)]=generate_phase_effect(0,0,0,0,0,0,0,offset,no_trials_bin,no_phasebins,no_phasebins_sampl,mode);
            
        end
    end
end

NoiseAUC = mean(NullAUC);

%% EFFECT PRESENT
AUC = nan(1,no_subjs);
R_sub = nan(no_trials_bin,no_subjs);
C_sub = nan(no_trials_bin,no_subjs);
for nSubs = 1:length(no_subjs)
    for t = 1:no_tests
        for subj = 1:no_subjs(nSubs)
            
            % determine "best" phase for each participant
%             randvec = allphases(randperm(no_phasebins_sampl));
            randvec = allphases(randperm(no_phasebins_sampl/4)+no_phasebins_sampl*3/4);
            whichphase_pos = randvec(1);

            % phase with trough is always 180 degrees away from best phase (phase with peak)
            whichphase_neg = whichphase_pos+(no_phasebins_sampl/2);
            if whichphase_neg > no_phasebins_sampl
                whichphase_neg = whichphase_neg-no_phasebins_sampl;
            end
            
            % now with effect present
            [phasetrials,totalperf, AUC(subj)]=generate_phase_effect(1,whichphase_pos,whichphase_neg,w_pos,w_neg,w_total,a_norm,offset,no_trials_bin,no_phasebins,no_phasebins_sampl,mode);
            
            R_sub(:,subj) = phasetrials(:,1);
            C_sub(:,subj) = phasetrials(:,4);
        end
    end
end

CohensD = mean(AUC-NoiseAUC)/std(AUC)

figure
histogram((AUC-NoiseAUC)/std(AUC))
xline(CohensD,'--r',['d = ' num2str(CohensD)],'LabelOrientation','horizontal')


for iSub = 1:size(R_sub,2)
    
    Diff(iSub) = mean(R_sub(C_sub(:,iSub) == 0,iSub))- mean(R_sub(C_sub(:,iSub) == -pi,iSub));
    
end

Diff_D = mean(Diff)/std(Diff)
figure
histogram(Diff/std(Diff))
xline(Diff_D,'--r',['d = ' num2str(Diff_D)],'LabelOrientation','horizontal')

% save all_simulations_binary.mat