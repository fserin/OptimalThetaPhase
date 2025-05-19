
%% function to generate a phasic modulation of dichotomous ('hit' vs 'miss') trial outcomes as described in Zoefel, Davis, Valente, Riecke, 2019.

%% required input: effect present or not (1/0), number of phase bin with peak, number of phase bin with trough,
%% width of positive peak, width of negative peak, total width, effect size, baseline, number of trials (/192), number of phase bins,
%% number of phase bins to characterize the true effect, study design

%% output 1: matrix with one line for each trial and the following columns:
%% 1: single-trial response
%% 2: phase that corresponds to the phase bin from column 4
%% 3: number of the phase bin that phase from column 4 falls into
%% 4: single-trial phase (actual phase out of "no_phasebins_sampl" possible phases)

%% output 2: response averaged across trials and divided into phase bins

function [phasetrials,total_effect_analysed, AUC] = generate_phase_effect(effectpresent,peak_pos,peak_neg,w_pos,w_neg,w_total,a_norm,offset,no_trials_bin,nbins_analysis,nbins_total,mode)

total_effect = zeros(nbins_total,1);

if effectpresent == 1
    
    %% in the following, "pos" and "neg" refer to the positive and negative deflection, respectively
    t = 0:1/nbins_total:1-1/nbins_total;
    
    % number of bins covered by the deflections
    nbins_pos = round(w_pos * nbins_total * w_total);
    nbins_neg = round(w_neg * nbins_total * w_total);
    
    % frequency that corresponds to these numbers of bins
    freq_pos = (nbins_total/2)/nbins_pos;
    freq_neg = (nbins_total/2)/nbins_neg;
    
    % amplitude necessary to achieve given area under curve
    a_pos = freq_pos * a_norm;
    a_neg = freq_neg * a_norm;
    
    % make the deflections, based on two sine waves
    sine_pos = a_pos*sin(2*pi*t*freq_pos);
    sine_neg = a_neg*sin(2*pi*t*freq_neg+pi);
    cycle_pos = sine_pos(1:nbins_pos);
    area_pos = sum(cycle_pos);
    cycle_neg = sine_neg(1:nbins_neg);
    area_neg = sum(cycle_pos);
    
    % now combine the two waves to make the phase effect
    if mod(nbins_pos,2) == 0
        pos_start = peak_pos-(nbins_pos/2);
        pos_end = peak_pos+(nbins_pos/2)-1;
    else
        pos_start = peak_pos-(nbins_pos-1)/2;
        pos_end = peak_pos+(nbins_pos-1)/2;
    end
    
    %% the following complicated stuff is to make sure everything is phase-wrapped
    if pos_start <= 0 || pos_end > nbins_total
        if pos_start <= 0
            pos_start = pos_start+nbins_total;
        end
        if pos_end > nbins_total
            pos_end = pos_end-nbins_total;
        end
        total_effect(pos_start:end) = cycle_pos(1:length(pos_start:nbins_total));
        total_effect(1:pos_end) = cycle_pos(length(pos_start:nbins_total)+1:end);
    else
        total_effect(pos_start:pos_end) = cycle_pos;
    end
    
    if mod(nbins_neg,2) == 0
        neg_start = peak_neg-(nbins_neg/2);
        neg_end = peak_neg+(nbins_neg/2)-1;
    else
        neg_start = peak_neg-(nbins_neg-1)/2;
        neg_end = peak_neg+(nbins_neg-1)/2;
    end
    
    if neg_start <= 0 || neg_end > nbins_total
        if neg_start <= 0
            neg_start = neg_start+nbins_total;
        end
        if neg_end > nbins_total
            neg_end = neg_end-nbins_total;
        end
        total_effect(neg_start:end) = cycle_neg(1:length(neg_start:nbins_total));
        total_effect(1:neg_end) = cycle_neg(length(neg_start:nbins_total)+1:end);
    else
        total_effect(neg_start:neg_end) = cycle_neg;
    end
end



%% total_effect is the hit probability in each bin
total_effect = total_effect+offset;

% figure
% plot(total_effect)
% ylim(0:1)
%% total_effect is the "true" effect. in the following, we sample this effect using the experimental parameters provided

allphases_total = -pi:2*pi/nbins_total:pi-2*pi/nbins_total; % phase "bins" used to make true effect
allphases_analysis = -pi:2*pi/nbins_analysis:pi-2*pi/nbins_analysis; % phase bins used for analysis
phasecounter = zeros(nbins_analysis,1);

phasetrials = zeros(no_trials_bin,4);
total_effect_analysed = zeros(nbins_analysis,1);

if mode == 1 %% randomly selected phases
    
    for trial = 1:no_trials_bin
        
        % randomly select phase for trial
        currp = randi(nbins_total);
        
        whichdiff = 1000;
        % select bin that this phase falls into
        for ph = 1:nbins_analysis
            currr = abs(circ_dist(allphases_total(currp),allphases_analysis(ph)));
            if currr < whichdiff
                whichdiff = currr;
                whichphase = ph;
            end
        end
        
        % for each trial, make decision (correct/incorrect) with the
        % corresponding probability
        
        phasecounter(whichphase) = phasecounter(whichphase)+1;
        currd = rand(1);
        if currd <= total_effect(currp)
            phasetrials(trial,1) = 1; % hit
            total_effect_analysed(whichphase,1) = total_effect_analysed(whichphase,1)+1;
        end
        
        phasetrials(trial,2) = allphases_analysis(whichphase);
        phasetrials(trial,3) = whichphase;
        phasetrials(trial,4) = allphases_total(currp);
    end
    
elseif mode == 2 %% externally imposed phase
    
    whichtrial = 0;
    
    for trial = 1:no_trials_bin/nbins_analysis
        
        %% only test phases corresponding to the bins used for analysis
        for whichphase = 1:nbins_analysis
            
            whichtrial = whichtrial+1;
            currp = allphases_analysis(whichphase);
            
            % select the phase (out of all possible 192 phases) that corresponds to the imposed phase
            whichdiff = 1000;
            for ph = 1:nbins_total
                currr = abs(circ_dist(allphases_total(ph),currp));
                if currr < whichdiff
                    whichdiff = currr;
                    whichphase2 = ph;
                end
            end
            
            %% for each trial, make decision (correct/incorrect) with the
            %% corresponding probability
            
            phasecounter(whichphase) = phasecounter(whichphase)+1;
            currd = rand(1);
            if currd <= total_effect(whichphase2)
                phasetrials(whichtrial,1) = 1;
                total_effect_analysed(whichphase,1) = total_effect_analysed(whichphase,1)+1;
            end
            
            phasetrials(whichtrial,2) = currp;
            phasetrials(whichtrial,3) = whichphase;
            phasetrials(whichtrial,4) = currp;
        end
    end
end

%% divide by number of phases tested in each bin to get detection
%% probability for each bin

for p = 1:nbins_analysis
    total_effect_analysed(p) = total_effect_analysed(p)./phasecounter(p);
end

%% Calculate area under curve (AUC)
% set scores below 0.5 to 0.5 so that areas below the offset are not
% included in AUC calculation.
% Use trapz function to calculate AUC
A = total_effect_analysed;
A(A < .5) = 0.5;
AUC = trapz(A);