% Magnitude-squared coherence analysis of vessel ROIs
%
% This analysis quantifies frequency-domain synchrony between multiple vessel
% diameter regions of interest (ROIs) within a recording. Each signal
% is converted from percent of baseline to deviation from baseline, optionally
% detrended, bandpass filtered in the vasomotor frequency range, and compared
% pairwise using magnitude-squared coherence.
%
% Coherence estimates how consistently two signals co-vary at each frequency.
% Values range from 0 to 1, where higher values indicate stronger frequency-
% specific coupling between ROI pairs.
%
% Input files (.csv) should be in standardized spreadsheet format:
%   column E = time
%   column F onward = vessel ROI diameter signals
%   row 12 onward = numeric data

inputFolder = 'input\directory\here'; % Replace with path to user input data
outputFolder = fullfile(inputFolder, 'vessel_coherence_output');

% Standardized input file format
dataStartRow = 12;
numHeaderLines = dataStartRow - 1;

timeCol = 5;
roiStartCol = 6;

Fs = 1.951; % Sampling frequency: 1.951 for 2Hz, 1.727 for 1.7Hz, 7.440476 for 7Hz
lowF = 0.025; % Frequency band (lower and upper limits) used to isolate vasomotor-range oscillations
highF = 0.25;

filterOrder = 2;
detrendData = true; % Optional detrending depending on how signals were preprocessed
trimEdgeSec = 5;

% Welch coherence settings
windowSec = 60;
overlapFraction = 0.5;

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Butterworth bandpass filter with filtfilt below for zero-phase filtering
[b, a] = butter(filterOrder, [lowF highF] / (Fs/2), 'bandpass');

pairwiseSummary = table();
fileSummary = table();

% Recursively locate all input .csv files in folders and subfolders
files = dir(fullfile(inputFolder, '**', '*.csv')); 

if isempty(files)
    error('No CSV files found in: %s', inputFolder);
end

for k = 1:numel(files)

    fname = files(k).name;
    fpath = fullfile(files(k).folder, fname);

    M = readmatrix(fpath, 'NumHeaderLines', numHeaderLines);

    if size(M,2) < roiStartCol + 1
        fprintf('Skipping %s: not enough ROI columns.\n', fname);
        continue;
    end

    t = M(:, timeCol);
    roiData = M(:, roiStartCol:end);

    % Keep rows where time and all ROI values are finite so all ROIs remain aligned
    validRows = isfinite(t) & all(isfinite(roiData), 2);
    t = t(validRows);
    roiData = roiData(validRows, :);

    if numel(t) < 20
        fprintf('Skipping %s: too few valid samples.\n', fname);
        continue;
    end

    durationSec = t(end) - t(1);
    nROIs = size(roiData, 2);

    % Assign ROI labels based on column order
    roiNames = strings(1, nROIs);
    for r = 1:nROIs
        roiNames(r) = "ROI_" + r;
    end

    trimMask = t >= (t(1) + trimEdgeSec) & t <= (t(end) - trimEdgeSec);

    roiProcessed = nan(size(roiData));

    % Preprocess each vessel ROI before coherence analysis
    for r = 1:nROIs
        x = roiData(:, r);
        % Input diameter data is percent of baseline, where 100 represents baseline
        % Subtracting 100 then converts it to deviation from baseline
        x_dev = x - 100;

        if detrendData
            x_dev = detrend(x_dev);
        end

        % Zero-phase bandpass filtering avoids introducing artificial phase shifts
        x_f = filtfilt(b, a, x_dev);
        roiProcessed(:, r) = x_f;
    end

    % Convert coherence window length from seconds to samples
    windowSamples = round(windowSec * Fs);

    if windowSamples > length(t)
        windowSamples = floor(length(t) / 2);
    end

    if windowSamples < 8
        fprintf('Skipping %s: recording too short for coherence window.\n', fname);
        continue;
    end

    % Hamming window and 50% overlap for Welch-style coherence estimation
    window = hamming(windowSamples);
    noverlap = round(windowSamples * overlapFraction);
    nfft = max(256, 2^nextpow2(windowSamples));

    [~, nameNoExt, ~] = fileparts(fname);

    allCoherenceCurves = [];
    allMeanCoherence = [];
    allPeakCoherence = [];
    allPeakFreq = [];

    pairCounter = 0;

    % Calculate magnitude-squared coherence for each ROI pair
    for i = 1:nROIs
        for j = i+1:nROIs

            xi = roiProcessed(trimMask, i);
            xj = roiProcessed(trimMask, j);

            [Cxy, F] = mscohere(xi, xj, window, noverlap, nfft, Fs);

            % Summarize coherence within the selected vasomotor band
            bandMask = F >= lowF & F <= highF;

            if ~any(bandMask)
                meanCoherenceInBand = NaN;
                peakCoherenceInBand = NaN;
                freqOfPeakCoherenceHz = NaN;
            else
                Cband = Cxy(bandMask);
                Fband = F(bandMask);

                meanCoherenceInBand = mean(Cband, 'omitnan');
                [peakCoherenceInBand, peakIdx] = max(Cband);
                freqOfPeakCoherenceHz = Fband(peakIdx);
            end

            % Store one summary row per ROI pair
            pairRow = table( ...
                string(fname), durationSec, ...
                string(roiNames(i)), string(roiNames(j)), ...
                lowF, highF, ...
                windowSec, overlapFraction, ...
                meanCoherenceInBand, peakCoherenceInBand, freqOfPeakCoherenceHz, ...
                'VariableNames', { ...
                'File','DurationSec', ...
                'ROI_1','ROI_2', ...
                'LowFreqHz','HighFreqHz', ...
                'WindowSec','OverlapFraction', ...
                'MeanCoherenceInBand','PeakCoherenceInBand','FreqOfPeakCoherenceHz'});

            pairwiseSummary = [pairwiseSummary; pairRow];

            pairCounter = pairCounter + 1;
            allCoherenceCurves(:, pairCounter) = Cxy;
            allMeanCoherence = [allMeanCoherence; meanCoherenceInBand];
            allPeakCoherence = [allPeakCoherence; peakCoherenceInBand];
            allPeakFreq = [allPeakFreq; freqOfPeakCoherenceHz];

            coherenceTable = table(F, Cxy, ...
                'VariableNames', {'FrequencyHz','Coherence'});

            pairName = nameNoExt + "_" + roiNames(i) + "_vs_" + roiNames(j);
            writetable(coherenceTable, fullfile(outputFolder, pairName + "_coherence_curve.xlsx"));

            % Plot pairwise coherence as a function of frequency
            figPair = figure('Visible','off');
            plot(F, Cxy, 'LineWidth', 1.2); hold on;
            xline(lowF, '--');
            xline(highF, '--');
            xlabel('Frequency (Hz)');
            ylabel('Magnitude-squared coherence');
            title(['Coherence: ' char(pairName)], 'Interpreter','none');
            ylim([0 1]);
            xlim([0 max(0.5, highF + 0.05)]);

            saveas(figPair, fullfile(outputFolder, pairName + "_coherence.png"));
            print(figPair, fullfile(outputFolder, pairName + "_coherence.eps"), '-depsc');
            close(figPair);

        end
    end

    % Average coherence curves across all ROI pairs within each recording
    if pairCounter > 0
        meanCoherenceCurve = mean(allCoherenceCurves, 2, 'omitnan');

        meanCurveTable = table(F, meanCoherenceCurve, ...
            'VariableNames', {'FrequencyHz','MeanCoherenceAcrossPairs'});

        writetable(meanCurveTable, fullfile(outputFolder, nameNoExt + "_mean_coherence_curve.xlsx"));

        figMean = figure('Visible','off');
        plot(F, meanCoherenceCurve, 'LineWidth', 1.5); hold on;
        xline(lowF, '--');
        xline(highF, '--');
        xlabel('Frequency (Hz)');
        ylabel('Mean magnitude-squared coherence');
        title(['Mean Coherence Across ROI Pairs: ' nameNoExt], 'Interpreter','none');
        ylim([0 1]);
        xlim([0 max(0.5, highF + 0.05)]);

        saveas(figMean, fullfile(outputFolder, nameNoExt + "_mean_coherence.png"));
        print(figMean, fullfile(outputFolder, nameNoExt + "_mean_coherence.eps"), '-depsc');
        close(figMean);
    end

    fileRow = table( ...
        string(fname), durationSec, nROIs, ...
        lowF, highF, ...
        windowSec, overlapFraction, ...
        mean(allMeanCoherence, 'omitnan'), median(allMeanCoherence, 'omitnan'), ...
        mean(allPeakCoherence, 'omitnan'), median(allPeakCoherence, 'omitnan'), ...
        mean(allPeakFreq, 'omitnan'), median(allPeakFreq, 'omitnan'), ...
        'VariableNames', { ...
        'File','DurationSec','NumROIs', ...
        'LowFreqHz','HighFreqHz', ...
        'WindowSec','OverlapFraction', ...
        'MeanCoherenceInBand','MedianCoherenceInBand', ...
        'MeanPeakCoherenceInBand','MedianPeakCoherenceInBand', ...
        'MeanFreqOfPeakCoherenceHz','MedianFreqOfPeakCoherenceHz'});

    fileSummary = [fileSummary; fileRow];

end

writetable(pairwiseSummary, fullfile(outputFolder, 'vessel_pairwise_coherence_summary.xlsx'));
writetable(fileSummary, fullfile(outputFolder, 'vessel_file_coherence_summary.xlsx'));

disp('Finished. Outputs saved to:');
disp(outputFolder);