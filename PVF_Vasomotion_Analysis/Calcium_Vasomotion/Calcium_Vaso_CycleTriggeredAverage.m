% Cycle-triggered average analysis for cellular calcium and vessel diameter 
% signals
%
% This analysis aligns calcium and diameter traces to detected peaks in the
% filtered vessel diameter signal. For each recording, individual vasomotor
% cycles are extracted around diameter peaks and averaged to generate
% cycle-triggered traces. Calcium peak and trough timing is then quantified
% relative to the corresponding diameter peaks and troughs.
%
% Input files (.csv) should be in standardized spreadsheet format:
%   column E = time
%   column F = calcium signal
%   column G = vessel diameter signal
%   row 12 onward = numeric data

inputFolder = 'input\directory\here'; % Replace with path to user input data
outputFolder = fullfile(inputFolder, 'cycle_triggered_average_output');

% Standardized input file format
dataStartRow = 12;
numHeaderLines = dataStartRow - 1;

timeCol = 5;
calciumCol = 6;
diameterCol = 7;

Fs = 1.951; % Sampling frequency: 1.951 for 2Hz, 1.727 for 1.7Hz, 7.440476 for 7Hz
lowF = 0.025;  % Frequency band (lower and upper limits) used to isolate vasomotor-range oscillations
highF = 0.25;
filterOrder = 2;

% Minimum recording duration and PSD-based quality-control thresholding
minDurationSec = 50;
minBandFractionCalcium = 0.02;
minBandFractionDiameter = 0.02;

% Peak detection and cycle window settings
minPeakDistanceSec = 8;
windowSec = 10;
minPeakProminence = 0.2;

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

summaryResults = table();

% Recursively locate all input .csv files in folders and subfolders
files = dir(fullfile(inputFolder, '**', '*.csv')); 
if isempty(files)
    error('No CSV files found in: %s', inputFolder);
end

% Butterworth bandpass filter with filtfilt below for zero-phase filtering
[b, a] = butter(filterOrder, [lowF highF] / (Fs/2), 'bandpass');

for k = 1:numel(files)
    fname = files(k).name;
    fpath = fullfile(files(k).folder, fname);

    M = readmatrix(fpath, 'NumHeaderLines', numHeaderLines);

    if size(M,2) < max([timeCol calciumCol diameterCol])
        fprintf('Skipping %s: not enough columns.\n', fname);
        continue;
    end

    t = M(:, timeCol);
    calcium = M(:, calciumCol);
    diameter = M(:, diameterCol);

    valid = isfinite(t) & isfinite(calcium) & isfinite(diameter);
    t = t(valid);
    calcium = calcium(valid);
    diameter = diameter(valid);

    if numel(t) < 10
        fprintf('Skipping %s: too few valid samples.\n', fname);
        continue;
    end

    durationSec = t(end) - t(1);
    if durationSec < minDurationSec
        fprintf('Skipping %s: duration too short.\n', fname);
        continue;
    end

    % Input diameter data is percent of baseline, where 100 represents baseline
    % Subtracting 100 then converts it to deviation from baseline
    diameter_dev = diameter - 100;

    % Detrend to remove slow linear drift before frequency-domain analysis and filtering
    calcium_dt = detrend(calcium);
    diameter_dt = detrend(diameter_dev);

    % Zero-phase bandpass filtering avoids introducing phase shifts between signals
    calcium_f = filtfilt(b, a, calcium_dt);
    diameter_f = filtfilt(b, a, diameter_dt);

    % Welch PSD is used here as a quality control (qc) to estimate how much signal power
    % falls within the vasomotor frequency band
    nfft = max(256, 2^nextpow2(numel(calcium_dt)));
    windowPSD = min(numel(calcium_dt), round(60 * Fs));
    if mod(windowPSD,2) ~= 0
        windowPSD = windowPSD - 1;
    end
    if windowPSD < 8
        windowPSD = min(numel(calcium_dt), 8);
    end
    noverlap = floor(windowPSD * 0.5);

    [PxxCal, Fcal] = pwelch(calcium_dt, windowPSD, noverlap, nfft, Fs);
    [PxxDiam, Fdiam] = pwelch(diameter_dt, windowPSD, noverlap, nfft, Fs);

    totalBandMaskCal = Fcal > 0 & Fcal <= (Fs/2);
    totalBandMaskDiam = Fdiam > 0 & Fdiam <= (Fs/2);

    vasoMaskCal = Fcal >= lowF & Fcal <= highF;
    vasoMaskDiam = Fdiam >= lowF & Fdiam <= highF;

    totalPowerCal = trapz(Fcal(totalBandMaskCal), PxxCal(totalBandMaskCal));
    totalPowerDiam = trapz(Fdiam(totalBandMaskDiam), PxxDiam(totalBandMaskDiam));

    vasoPowerCal = trapz(Fcal(vasoMaskCal), PxxCal(vasoMaskCal));
    vasoPowerDiam = trapz(Fdiam(vasoMaskDiam), PxxDiam(vasoMaskDiam));

    if totalPowerCal > 0
        bandFractionCal = vasoPowerCal / totalPowerCal;
    else
        bandFractionCal = NaN;
    end

    if totalPowerDiam > 0
        bandFractionDiam = vasoPowerDiam / totalPowerDiam;
    else
        bandFractionDiam = NaN;
    end

    qcPassCal = bandFractionCal >= minBandFractionCalcium;
    qcPassDiam = bandFractionDiam >= minBandFractionDiameter;
    qcPassBand = qcPassCal && qcPassDiam;

    % Dominant diameter frequency is used to estimate the representative cycle period
    if any(vasoMaskDiam)
        [~, idxPeakDiamPSD] = max(PxxDiam(vasoMaskDiam));
        freqValsDiam = Fdiam(vasoMaskDiam);
        dominantFreqHz = freqValsDiam(idxPeakDiamPSD);
    else
        dominantFreqHz = NaN;
    end

    if ~isnan(dominantFreqHz) && dominantFreqHz > 0
        cyclePeriodSec = 1 / dominantFreqHz;
    else
        cyclePeriodSec = NaN;
    end

    minPeakDistanceSamples = round(minPeakDistanceSec * Fs);
    halfWin = round((windowSec / 2) * Fs);

    % Detect peaks in filtered diameter - used to define the alignment
    % points for the cycle-triggered average
    if isempty(minPeakProminence)
        [pks, locs] = findpeaks(diameter_f, 'MinPeakDistance', minPeakDistanceSamples);
    else
        [pks, locs] = findpeaks(diameter_f, 'MinPeakDistance', minPeakDistanceSamples, ...
            'MinPeakProminence', minPeakProminence);
    end

    % Exclude peaks too close to the beginning or end to contain a full window
    keep = locs - halfWin >= 1 & locs + halfWin <= numel(diameter_f);
    locs = locs(keep);
    pks = pks(keep);

    nCycles = numel(locs);

    if nCycles < 3
        fprintf('Skipping %s: not enough usable cycles.\n', fname);
        continue;
    end

    calciumCycles = NaN(nCycles, 2*halfWin + 1);
    diameterCycles = NaN(nCycles, 2*halfWin + 1);

    % Extract calcium and diameter windows centered around each diameter peak
    for i = 1:nCycles
        idx = locs(i);
        calciumCycles(i,:) = calcium_f(idx-halfWin : idx+halfWin);
        diameterCycles(i,:) = diameter_f(idx-halfWin : idx+halfWin);
    end

    % Average across detected cycles to estimate the representative waveform
    meanCalcium = mean(calciumCycles, 1, 'omitnan');
    meanDiameter = mean(diameterCycles, 1, 'omitnan');

    semCalcium = std(calciumCycles, 0, 1, 'omitnan') ./ sqrt(nCycles);
    semDiameter = std(diameterCycles, 0, 1, 'omitnan') ./ sqrt(nCycles);

    relTime = (-halfWin:halfWin) / Fs;

    % Identify global peak and trough positions in the averaged cycle
    [diamPeakVal, diamPeakIdx] = max(meanDiameter);
    [diamTroughVal, diamTroughIdx] = min(meanDiameter);

    [calPeakVal, calPeakIdx] = max(meanCalcium);
    [calTroughVal, calTroughIdx] = min(meanCalcium);

    diamPeakTimeSec = relTime(diamPeakIdx);
    diamTroughTimeSec = relTime(diamTroughIdx);
    calPeakTimeSec = relTime(calPeakIdx);
    calTroughTimeSec = relTime(calTroughIdx);

    % Compute timing offsets
    calPeakLagFromDiamPeakSec = calPeakTimeSec - diamPeakTimeSec;
    calPeakLagFromDiamTroughSec = calPeakTimeSec - diamTroughTimeSec;
    calTroughLagFromDiamPeakSec = calTroughTimeSec - diamPeakTimeSec;
    calTroughLagFromDiamTroughSec = calTroughTimeSec - diamTroughTimeSec;

    if abs(calPeakLagFromDiamPeakSec) <= abs(calPeakLagFromDiamTroughSec)
        calciumPeakCloserTo = "DiameterPeak";
    else
        calciumPeakCloserTo = "DiameterTrough";
    end

    if abs(calTroughLagFromDiamPeakSec) <= abs(calTroughLagFromDiamTroughSec)
        calciumTroughCloserTo = "DiameterPeak";
    else
        calciumTroughCloserTo = "DiameterTrough";
    end

    % Identify local (not global) calcium peaks and troughs on each side of 
    % time zero to capture that may be biologically relevant
    [calLocalPeaks, calPeakLocs] = findpeaks(meanCalcium);

    [calLocalTroughsNeg, calTroughLocs] = findpeaks(-meanCalcium);
    calLocalTroughs = -calLocalTroughsNeg;

    peakBeforeLocs = calPeakLocs(relTime(calPeakLocs) < 0);
    peakAfterLocs = calPeakLocs(relTime(calPeakLocs) > 0);

    if ~isempty(peakBeforeLocs)
        [~, idxPeakBefore] = min(abs(relTime(peakBeforeLocs)));
        closestCalciumPeakBeforeZeroTimeSec = relTime(peakBeforeLocs(idxPeakBefore));
        closestCalciumPeakBeforeZeroValue = meanCalcium(peakBeforeLocs(idxPeakBefore));
    else
        closestCalciumPeakBeforeZeroTimeSec = NaN;
        closestCalciumPeakBeforeZeroValue = NaN;
    end

    if ~isempty(peakAfterLocs)
        [~, idxPeakAfter] = min(abs(relTime(peakAfterLocs)));
        closestCalciumPeakAfterZeroTimeSec = relTime(peakAfterLocs(idxPeakAfter));
        closestCalciumPeakAfterZeroValue = meanCalcium(peakAfterLocs(idxPeakAfter));
    else
        closestCalciumPeakAfterZeroTimeSec = NaN;
        closestCalciumPeakAfterZeroValue = NaN;
    end

    troughBeforeLocs = calTroughLocs(relTime(calTroughLocs) < 0);
    troughAfterLocs = calTroughLocs(relTime(calTroughLocs) > 0);

    if ~isempty(troughBeforeLocs)
        [~, idxTroughBefore] = min(abs(relTime(troughBeforeLocs)));
        closestCalciumTroughBeforeZeroTimeSec = relTime(troughBeforeLocs(idxTroughBefore));
        closestCalciumTroughBeforeZeroValue = meanCalcium(troughBeforeLocs(idxTroughBefore));
    else
        closestCalciumTroughBeforeZeroTimeSec = NaN;
        closestCalciumTroughBeforeZeroValue = NaN;
    end

    if ~isempty(troughAfterLocs)
        [~, idxTroughAfter] = min(abs(relTime(troughAfterLocs)));
        closestCalciumTroughAfterZeroTimeSec = relTime(troughAfterLocs(idxTroughAfter));
        closestCalciumTroughAfterZeroValue = meanCalcium(troughAfterLocs(idxTroughAfter));
    else
        closestCalciumTroughAfterZeroTimeSec = NaN;
        closestCalciumTroughAfterZeroValue = NaN;
    end

    cycleTable = array2table([relTime(:), meanCalcium(:), semCalcium(:), meanDiameter(:), semDiameter(:)], ...
        'VariableNames', {'RelativeTimeSec','MeanCalcium','SEMCalcium','MeanDiameter','SEMDiameter'});

    individualCycleTable = table();
    for i = 1:nCycles
        Ti = table( ...
            repmat(i, numel(relTime), 1), ...
            relTime(:), ...
            calciumCycles(i,:)', ...
            diameterCycles(i,:)', ...
            'VariableNames', {'CycleNumber','RelativeTimeSec','Calcium','Diameter'});
        individualCycleTable = [individualCycleTable; Ti];
    end

    [~, nameNoExt, ~] = fileparts(fname);

    writetable(cycleTable, fullfile(outputFolder, [nameNoExt '_cycle_triggered_average.csv']));
    writetable(individualCycleTable, fullfile(outputFolder, [nameNoExt '_individual_cycles.csv']));

    % Plot filtered calcium and diameter traces with detected diameter peaks
    fig1 = figure('Visible', 'off');

    ax1 = subplot(2,1,1);
    plot(t, calcium_f, 'LineWidth', 1.2);
    ylabel('Filtered calcium');
    title(['Filtered Signals with Detected Diameter Peaks: ' nameNoExt], 'Interpreter', 'none');

    ax2 = subplot(2,1,2);
    plot(t, diameter_f, 'LineWidth', 1.2); hold on;
    scatter(t(locs), diameter_f(locs), 30, 'filled');
    xlabel('Time (s)');
    ylabel('Filtered diameter');

    linkaxes([ax1 ax2], 'x');

    saveas(fig1, fullfile(outputFolder, [nameNoExt '_detected_cycles.png']));
    print(fig1, fullfile(outputFolder, [nameNoExt '_detected_cycles.eps']), '-depsc');
    close(fig1);

    % Plot the cycle-triggered average with SEM around the mean waveform
    fig2 = figure('Visible', 'off');

    ax3 = subplot(2,1,1);
    plot(relTime, meanCalcium, 'LineWidth', 2); hold on;
    plot(relTime, meanCalcium + semCalcium, '--', 'LineWidth', 1);
    plot(relTime, meanCalcium - semCalcium, '--', 'LineWidth', 1);
    xline(0, '--');
    yline(calPeakVal, ':');
    yline(calTroughVal, ':');
    ylabel('Calcium');
    title(['Cycle-Triggered Average: ' nameNoExt], 'Interpreter', 'none');

    ax4 = subplot(2,1,2);
    plot(relTime, meanDiameter, 'LineWidth', 2); hold on;
    plot(relTime, meanDiameter + semDiameter, '--', 'LineWidth', 1);
    plot(relTime, meanDiameter - semDiameter, '--', 'LineWidth', 1);
    xline(0, '--');
    yline(diamPeakVal, ':');
    yline(diamTroughVal, ':');
    xlabel('Time relative to diameter peak (s)');
    ylabel('Diameter');

    linkaxes([ax3 ax4], 'x');

    saveas(fig2, fullfile(outputFolder, [nameNoExt '_cycle_triggered_average.png']));
    print(fig2, fullfile(outputFolder, [nameNoExt '_cycle_triggered_average.eps']), '-depsc');
    close(fig2);

    % Save PSD plots used for visual qc of vasomotor-band power
    fig3 = figure('Visible', 'off');

    ax5 = subplot(2,1,1);
    plot(Fcal, PxxCal, 'LineWidth', 1.2); hold on;
    xline(lowF, '--');
    xline(highF, '--');
    ylabel('Calcium PSD');
    title(['PSD and Vasomotion Band: ' nameNoExt], 'Interpreter', 'none');
    xlim([0 max(0.5, highF + 0.05)]);

    ax6 = subplot(2,1,2);
    plot(Fdiam, PxxDiam, 'LineWidth', 1.2); hold on;
    xline(lowF, '--');
    xline(highF, '--');
    xlabel('Frequency (Hz)');
    ylabel('Diameter PSD');
    xlim([0 max(0.5, highF + 0.05)]);

    linkaxes([ax5 ax6], 'x');

    saveas(fig3, fullfile(outputFolder, [nameNoExt '_psd_qc.png']));
    print(fig3, fullfile(outputFolder, [nameNoExt '_psd_qc.eps']), '-depsc');
    close(fig3);

    % Store one summary row per recording
    newRow = table( ...
        string(fname), durationSec, ...
        bandFractionCal, bandFractionDiam, qcPassCal, qcPassDiam, qcPassBand, ...
        dominantFreqHz, cyclePeriodSec, ...
        nCycles, ...
        diamPeakVal, diamPeakTimeSec, ...
        diamTroughVal, diamTroughTimeSec, ...
        calPeakVal, calPeakTimeSec, ...
        calTroughVal, calTroughTimeSec, ...
        calPeakLagFromDiamPeakSec, calPeakLagFromDiamTroughSec, ...
        calTroughLagFromDiamPeakSec, calTroughLagFromDiamTroughSec, ...
        string(calciumPeakCloserTo), string(calciumTroughCloserTo), ...
        closestCalciumPeakBeforeZeroValue, closestCalciumPeakBeforeZeroTimeSec, ...
        closestCalciumPeakAfterZeroValue, closestCalciumPeakAfterZeroTimeSec, ...
        closestCalciumTroughBeforeZeroValue, closestCalciumTroughBeforeZeroTimeSec, ...
        closestCalciumTroughAfterZeroValue, closestCalciumTroughAfterZeroTimeSec, ...
        'VariableNames', { ...
        'File','DurationSec', ...
        'BandFractionCalcium','BandFractionDiameter','QCPassCalcium','QCPassDiameter','QCPassBand', ...
        'DominantFreqHz','CyclePeriodSec', ...
        'NumCycles', ...
        'MeanDiameterPeakValue','MeanDiameterPeakTimeSec', ...
        'MeanDiameterTroughValue','MeanDiameterTroughTimeSec', ...
        'MeanCalciumPeakValue','MeanCalciumPeakTimeSec', ...
        'MeanCalciumTroughValue','MeanCalciumTroughTimeSec', ...
        'CalciumPeakLagFromDiameterPeakSec','CalciumPeakLagFromDiameterTroughSec', ...
        'CalciumTroughLagFromDiameterPeakSec','CalciumTroughLagFromDiameterTroughSec', ...
        'CalciumPeakCloserTo','CalciumTroughCloserTo', ...
        'ClosestCalciumPeakBeforeZeroValue','ClosestCalciumPeakBeforeZeroTimeSec', ...
        'ClosestCalciumPeakAfterZeroValue','ClosestCalciumPeakAfterZeroTimeSec', ...
        'ClosestCalciumTroughBeforeZeroValue','ClosestCalciumTroughBeforeZeroTimeSec', ...
        'ClosestCalciumTroughAfterZeroValue','ClosestCalciumTroughAfterZeroTimeSec'});

    summaryResults = [summaryResults; newRow];
end

writetable(summaryResults, fullfile(outputFolder, 'cycle_triggered_average_summary.xlsx'));

disp('Finished. Outputs saved to:');
disp(outputFolder);