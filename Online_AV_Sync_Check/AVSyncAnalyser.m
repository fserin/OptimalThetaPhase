% Audio-Visual Synchronization Analyser - Cross-Correlation Method
% Designed for: 192 trials, 3000ms each, alternating white/black every 2 frames @ 60Hz
% Uses cross-correlation for maximum precision

clear; clc; close all;

% Experimental parameters
num_trials = 192;
trial_duration_ms = 3000;
screen_refresh_rate = 60; % Hz
frames_per_flash = 2;
flash_duration_ms = (frames_per_flash / screen_refresh_rate) * 1000; % ~33.3 ms

% Load the audio file
[file, path] = uigetfile('*.wav', 'Select a WAV file');
if isequal(file, 0)
    disp('User canceled file selection');
    return;
end

fullPath = fullfile(path, file);
[audio, fs] = audioread(fullPath);

% Check if stereo
[numSamples, numChannels] = size(audio);
if numChannels < 2
    error('Audio file must be stereo (2 channels)');
end

% Define channels (adjust if needed)
photodiode = audio(:, 1);  % Visual stimulus (photodiode signal)
audioSignal = audio(:, 2);  % Audio stimulus (1000 Hz beep)

% Time vector
t = (0:numSamples-1)' / fs;

%% Detect Trials using AUDIO signal (photodiode has more noise and tend to differ in duration)
% Bandpass filter audio around 1000 Hz for FULL RECORDING
[b, a] = butter(4, [900 1100]/(fs/2), 'bandpass');
audioSignal_filtered = filtfilt(b, a, audioSignal);

% Calculate audio envelope for FULL RECORDING
audioSignal_envelope = abs(hilbert(audioSignal_filtered));

% Use adaptive threshold based on envelope distribution
sorted_env = sort(audioSignal_envelope);
noise_floor = sorted_env(round(0.25 * length(sorted_env))); % 25th percentile
signal_level = sorted_env(round(0.75 * length(sorted_env))); % 75th percentile
trial_threshold = noise_floor + 0.001 * (signal_level - noise_floor); % Adjust multiplier if needed

fprintf('Audio Noise floor: %.6f, Signal level: %.6f, Threshold: %.6f\n', ...
    noise_floor, signal_level, trial_threshold);

% Detect periods above threshold
in_trial = audioSignal_envelope > trial_threshold;

% Find trial start and end points
trial_changes = diff([0; in_trial; 0]);
trial_starts_raw = find(trial_changes == 1);
trial_ends_raw = find(trial_changes == -1) - 1;

% Merge trials separated by less than inter-trial interval (1000ms)
min_gap = round(0.5 * 1000 / 1000 * fs); % 500ms
trial_starts = trial_starts_raw(1);
trial_ends = [];

% NEW APPROACH: Use only trial END points, then go back 3100ms and forward 50ms
fprintf('\nDetecting trial END points...\n');

for i = 1:length(trial_starts_raw)
    if i < length(trial_starts_raw)
        gap = trial_starts_raw(i+1) - trial_ends_raw(i);
        if gap > min_gap
            trial_ends = [trial_ends; trial_ends_raw(i)];
            trial_starts = [trial_starts; trial_starts_raw(i+1)];
        end
    else
        trial_ends = [trial_ends; trial_ends_raw(i)];
    end
end

fprintf('Found %d trial end points\n', length(trial_ends));

% Epoch trials: go back 3100ms from end, forward 50ms from end
epoch_back_ms = 3100;
epoch_forward_ms = 100;
epoch_back_samples = round(epoch_back_ms / 1000 * fs);
epoch_forward_samples = round(epoch_forward_ms / 1000 * fs);

trial_starts = zeros(length(trial_ends), 1);
trial_ends_filtered = zeros(length(trial_ends), 1);

for i = 1:length(trial_ends_filtered)
    end_point = trial_ends(i);
    
    % Calculate start and end of epoch
    trial_starts(i) = max(1, end_point - epoch_back_samples);
    trial_ends_filtered(i) = min(numSamples, end_point + epoch_forward_samples);
end

num_detected_trials = length(trial_starts);

% Display trial durations
fprintf('\nDetected trial durations (ms): ');
for i = 1:num_detected_trials
    fprintf('%.0f ', (trial_ends_filtered(i) - trial_starts(i)) / fs * 1000);
end
fprintf('\n');

fprintf('\n===== Audio-Visual Synchronization Analysis =====\n');
fprintf('File: %s\n', file);
fprintf('Sample Rate: %d Hz\n', fs);
fprintf('Duration: %.2f seconds\n', numSamples/fs);
fprintf('\n--- Trial Detection ---\n');
fprintf('Expected trials: %d\n', num_trials);
fprintf('Detected trials: %d\n', num_detected_trials);

%% Analyze each trial using cross-correlation
trial_delays = zeros(num_detected_trials, 1);
trial_correlations = zeros(num_detected_trials, 1);
trial_stats = struct();

% Maximum lag to search ( should be more than enough)
max_lag_samples = round(0.2 * fs); % 

fprintf('\n--- Cross-Correlation Analysis per Trial ---\n');

for trial = 1:num_detected_trials
    % Extract trial segment
    trial_start_idx = trial_starts(trial);
    trial_end_idx = trial_ends_filtered(trial);
    
    pd_trial = photodiode(trial_start_idx:trial_end_idx);
    audio_trial = audioSignal(trial_start_idx:trial_end_idx);
    
    % Preprocess audio for THIS TRIAL: bandpass filter around 1000 Hz
    if fs/2 > 1100
        [b, a] = butter(4, [900 1100]/(fs/2), 'bandpass');
        audio_trial_filtered = filtfilt(b, a, audio_trial);
    else
        [b, a] = butter(4, 900/(fs/2), 'high');
        audio_trial_filtered = filtfilt(b, a, audio_trial);
    end
    
    % Get envelope of audio for the current
    audio_trial_envelope = abs(hilbert(audio_trial_filtered));
    
    % Normalize signals for cross-correlation
    pd_norm = (pd_trial - mean(pd_trial)) / std(pd_trial);
    audio_norm = (audio_trial_envelope - mean(audio_trial_envelope)) / std(audio_trial_envelope);
    
    % Perform cross-correlation
    [corr_values, lags] = xcorr(pd_norm, audio_norm, max_lag_samples, 'coeff');
    
    % Find peak correlation
    [max_corr, max_idx] = max(corr_values);
    optimal_lag = lags(max_idx);
    
    % Convert lag to milliseconds (positive = audio lags behind visual)
    delay_ms = optimal_lag / fs * 1000;
    
    % Store results
    trial_delays(trial) = delay_ms;
    trial_correlations(trial) = max_corr;
    
    trial_stats(trial).delay_ms = delay_ms;
    trial_stats(trial).correlation = max_corr;
    trial_stats(trial).lag_samples = optimal_lag;
    
    fprintf('Trial %d: Start=%.2fs, End=%.2fs, Delay=%+.3fms, Correlation=%.4f\n', ...
        trial, trial_start_idx / fs, trial_end_idx / fs, delay_ms, max_corr);
end

%% Overall Statistics
mean_delay = mean(trial_delays);
std_delay = std(trial_delays);
variance_delay = var(trial_delays);
median_delay = median(trial_delays);
min_delay = min(trial_delays);
max_delay = max(trial_delays);
range_delay = max_delay - min_delay;

mean_correlation = mean(trial_correlations);
min_correlation = min(trial_correlations);

fprintf('\n--- Overall Timing Analysis (Across Trials) ---\n');
fprintf('Mean delay (Audio - Visual): %.3f ms\n', mean_delay);
fprintf('Median delay: %.3f ms\n', median_delay);
fprintf('Standard Deviation (SD): %.3f ms\n', std_delay);
fprintf('Variance: %.3f ms²\n', variance_delay);
fprintf('Min delay: %.3f ms\n', min_delay);
fprintf('Max delay: %.3f ms\n', max_delay);
fprintf('Range: %.3f ms\n', range_delay);
fprintf('\n--- Correlation Quality ---\n');
fprintf('Mean correlation coefficient: %.4f\n', mean_correlation);
fprintf('Min correlation coefficient: %.4f\n', min_correlation);

if mean_delay > 0
    fprintf('\n→ Audio lags behind visual by %.3f ms on average\n', mean_delay);
elseif mean_delay < 0
    fprintf('\n→ Visual lags behind audio by %.3f ms on average\n', abs(mean_delay));
else
    fprintf('\n→ Perfect synchronization detected\n');
end

%% Visualizations
% Main figure with 6 plots
figure('Position', [50, 50, 1400, 900]);

% % Plot 1: Full recording with trial boundaries
% subplot(2, 3, 1);
% yyaxis left
% plot(t, audioSignal, 'r', 'LineWidth', 0.5);
% ylabel('Audio Signal');
% yyaxis right
% plot(t, audioSignal_envelope, 'k', 'LineWidth', 0.5);
% yline(trial_threshold, 'r--', 'Threshold', 'LineWidth', 1.5);
% ylabel('Audio Envelope');
% hold on;
% for trial = 1:min(num_detected_trials, 20) % Show first 20 trial markers to avoid clutter
%     xline(t(trial_starts(trial)), 'g-', 'LineWidth', 1);
% end
% xlabel('Time (s)');
% title('Full Recording with Trial Detection (Audio-Based)');
% grid on;
% 
% % Plot 2: Photodiode signal (for reference, may be noisy)
% subplot(2, 3, 2);
% plot(t, photodiode, 'b', 'LineWidth', 0.5);
% hold on;
% for trial = 1:min(num_detected_trials, 20)
%     xline(t(trial_starts(trial)), 'g--', 'LineWidth', 1);
% end
% xlabel('Time (s)');
% ylabel('Photodiode Signal');
% title('Photodiode Recording (may contain noise)');
% grid on;

% Plot 3: Delay per trial
subplot(2, 1, 1);
plot(1:num_detected_trials, trial_delays, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
yline(mean_delay, 'r--', sprintf('Mean = %.3f ms', mean_delay), 'LineWidth', 2);
yline(mean_delay + std_delay, 'b--', '+1 SD');
yline(mean_delay - std_delay, 'b--', '-1 SD');
xlabel('Trial Number');
ylabel('Delay (ms)');
title('Cross-Correlation Delay per Trial');
grid on;

% Plot 4: Delay histogram
subplot(2, 3, 4);
histogram(trial_delays, 30, 'FaceColor', 'b', 'EdgeColor', 'k');
hold on;
xline(mean_delay, 'r--', 'Mean', 'LineWidth', 2);
xline(median_delay, 'g--', 'Median', 'LineWidth', 2);
xlabel('Delay (ms)');
ylabel('Frequency');
title('Delay Distribution Across Trials');
grid on;

% Plot 5: Correlation coefficient per trial
subplot(2, 3, 5);
plot(1:num_detected_trials, trial_correlations, 'k-', 'LineWidth', 1);
hold on;
yline(mean_correlation, 'r--', sprintf('Mean = %.3f', mean_correlation), 'LineWidth', 2);
xlabel('Trial Number');
ylabel('Correlation Coefficient');
title('Cross-Correlation Quality per Trial');
ylim([0, 1]);
grid on;

% Plot 6: Delay vs Correlation
subplot(2, 3, 6);
scatter(trial_correlations, trial_delays, 100, 'b', 'filled', 'MarkerEdgeColor', 'k');
xlabel('Correlation Coefficient');
ylabel('Delay (ms)');
title('Delay vs Correlation Quality');
grid on;

sgtitle(sprintf('Audio-Visual Sync Analysis (Cross-Correlation): %s | SD = %.3f ms, Variance = %.3f ms²', ...
    file, std_delay, variance_delay), 'Interpreter', 'none');

%% Individual trial figures with combined top subplot

plot_trials = find(((mean_delay-trial_delays)/std_delay > 1.5) | ((mean_delay-trial_delays)/std_delay < -1.5)); %1:num_detected_trials

for i = 1:length(plot_trials)
    trial = plot_trials(i);
    figure;
    
    trial_start_idx = trial_starts(trial);
    trial_end_idx = trial_ends_filtered(trial);
    trial_idx = trial_start_idx:trial_end_idx;
    t_trial = (0:length(trial_idx)-1)' / fs;
    
    pd_trial = photodiode(trial_idx);
    audio_trial = audioSignal(trial_idx);
    
    % Preprocess audio
    if fs/2 > 1100
        [b, a] = butter(4, [900 1100]/(fs/2), 'bandpass');
        audio_trial_filtered = filtfilt(b, a, audio_trial);
    else
        [b, a] = butter(4, 900/(fs/2), 'high');
        audio_trial_filtered = filtfilt(b, a, audio_trial);
    end
    audio_trial_envelope = abs(hilbert(audio_trial_filtered));
    
    % Normalize
    pd_norm = (pd_trial - mean(pd_trial)) / std(pd_trial);
    audio_norm = (audio_trial_envelope - mean(audio_trial_envelope)) / std(audio_trial_envelope);
    
    % Cross-correlation
    [corr_values, lags] = xcorr(pd_norm, audio_norm, max_lag_samples, 'coeff');
    lags_ms = lags / fs * 1000;
        
    % Top subplot: combined photodiode + audio envelope
    subplot(2, 1, 1);
    yyaxis left
    plot(t_trial, pd_trial, 'Color', [0 0 1], 'LineWidth', 1.4);
    ylabel('Photodiode');
    
    yyaxis right
    plot(t_trial, audio_trial_envelope, 'Color', [1 0 0], 'LineWidth', 1.6);
    ylabel('Audio Envelope');
    
    xlabel('Time (s)');
    title(sprintf('Trial %d - Photodiode + Audio Envelope', trial));
    grid on;
    
    % Bottom-left: Cross-correlation function
    subplot(2, 2, 3);
    plot(lags_ms, corr_values, 'Color', [0 0 0], 'LineWidth', 1.8);
    hold on;
    plot(trial_stats(trial).delay_ms, trial_stats(trial).correlation, ...
        'o', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [1 0 0]);
    xlabel('Lag (ms)');
    ylabel('Correlation');
    title('Cross-Correlation');
    grid on;
    xlim([-100 100]);
    
    % Bottom-right: Aligned normalized signals
    subplot(2, 2, 4);
    plot(t_trial, pd_norm, 'Color', [0 0 1], 'LineWidth', 1.6);
    hold on;
    
    lag_samples = trial_stats(trial).lag_samples;
    if lag_samples > 0
        audio_shifted = [zeros(lag_samples, 1); audio_norm(1:end-lag_samples)];
    elseif lag_samples < 0
        audio_shifted = [audio_norm(-lag_samples+1:end); zeros(-lag_samples, 1)];
    else
        audio_shifted = audio_norm;
    end
    
    plot(t_trial, audio_shifted, 'Color', [1 0 0], 'LineWidth', 1.6);
    xlabel('Time (s)');
    ylabel('Normalized');
    title('Aligned Signals (Optimal shift)');
    legend('Photodiode', 'Audio (shifted)');
    grid on;
    
    sgtitle(sprintf('Trial %d: Delay = %+.3f ms, Corr = %.4f', ...
        trial, trial_stats(trial).delay_ms, trial_stats(trial).correlation), ...
        'Interpreter', 'none');
end


fprintf('\nAnalysis complete!\n');