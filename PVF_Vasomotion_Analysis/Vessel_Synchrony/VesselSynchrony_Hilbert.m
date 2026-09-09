% Hilbert phase synchrony analysis across vessel ROIs
%
% This analysis quantifies phase synchrony between multiple vessel diameter
% regions of interest (ROIs) within a recording. Each ROI signal is
% converted from percent of baseline to deviation from baseline, detrended,
% bandpass filtered in the vasomotor frequency range, and transformed using
% the Hilbert transform.
%
% Pairwise phase differences are calculated between all vessel ROIs, returning
% PLV, circular mean phase difference, and estimated lag for
% each ROI pair, along with full ROI-by-ROI matrices for each recording.
%
% Input files (.csv) should be in standardized spreadsheet format:
%   column E = time
%   column F onward = vessel ROI diameter signals
%   row 12 onward = numeric data

inputFolder = 'input\directory\here'; % Replace with path to user input data
outputFolder = fullfile(inputFolder, 'vessel_sync_output');

% Standardized input file format
dataStartRow = 12;
numHeaderLines = dataStartRow - 1;

timeCol = 5;
roiStartCol = 6;

Fs = 1.951; % Sampling frequency: 1.951 for 2Hz, 1.727 for 1.7Hz, 7.440476 for 7Hz
lowF = 0.025; % Frequency band (lower and upper limits) used to isolate vasomotor-range oscillations
highF = 0.25;
filterOrder = 2;
trimEdgeSec = 5; % Edge trimming removes filter/Hilbert edge artifacts

expectedCyclePeriodSec = 10; % Expected vasomotor cycle period used as a standardized lag estimate
minPeakDistanceSec = 3; % Minimum distance between detected diameter peaks for estimating each file-specific cycle period

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Butterworth bandpass filter with filtfilt below for zero-phase filtering
[b, a] = butter(filterOrder, [lowF highF] / (Fs/2), 'bandpass');

summaryResults = table();

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

    % Exclude set amount of data at the beginning and end of the trace where filtering and Hilbert
    % transform estimates are prone to edge artifacts
    trimMask = t >= (t(1) + trimEdgeSec) & t <= (t(end) - trimEdgeSec);

    if sum(trimMask) < 10
        trimMask = true(size(t));
    end

    roiFiltered = nan(size(roiData));
    roiPhaseRad = nan(size(roiData));

    % Preprocess each vessel ROI and extract instantaneous phase
    for r = 1:nROIs
        x = roiData(:, r);

        % Input diameter data is percent of baseline, where 100 represents baseline
        % Subtracting 100 then converts it to deviation from baseline
        x_dev = x - 100;
        x_dt = detrend(x_dev); % Remove slow linear drift before phase estimation
        x_f = filtfilt(b, a, x_dt); % Zero-phase bandpass filtering avoids introducing artificial phase shifts

        roiFiltered(:, r) = x_f;
        roiPhaseRad(:, r) = angle(hilbert(x_f));
    end

    % Use the mean filtered vessel signal to estimate the representative
    % oscillation period
    meanSignal = mean(roiFiltered, 2, 'omitnan');

    minPeakDistanceSamples = round(minPeakDistanceSec * Fs);

    [~, locs] = findpeaks(meanSignal, ...
        'MinPeakDistance', minPeakDistanceSamples);

    if numel(locs) >= 2
        peakTimes = t(locs);
        cyclePeriods = diff(peakTimes);
        fileCyclePeriodSec = median(cyclePeriods, 'omitnan');
    else
        fileCyclePeriodSec = NaN;
    end

    % Initialize ROI-by-ROI matrices
    PLVmatrix = nan(nROIs, nROIs);
    phaseDiffMatrixDeg = nan(nROIs, nROIs);
    lagSecMatrix = nan(nROIs, nROIs);

    % Calculate pairwise phase synchrony for every ROI combination,
    % including diagonal self-comparisons
    for i = 1:nROIs
        for j = 1:nROIs

            % Phase difference is wrapped to the circular range -pi to pi.
            % Positive values reflect ROI_i phase minus ROI_j phase
            phaseDiffRad = angle(exp(1i * (roiPhaseRad(:, i) - roiPhaseRad(:, j))));
            phaseDiffRad_use = phaseDiffRad(trimMask);

            % Circular mean gives the preferred phase offset between the two ROIs
            circMeanRad = atan2(mean(sin(phaseDiffRad_use)), mean(cos(phaseDiffRad_use)));
            circMeanDeg = rad2deg(circMeanRad);

            % PLV is the length of the mean resultant vector of phase differences.
            % Values near 1 indicate consistent phase locking between compared ROIs
            plv = abs(mean(exp(1i * phaseDiffRad_use)));

            % Convert phase offset into an approximate timing lag using the
            % file-specific cycle period. This time-based measure is
            % included for exploratory interpretation and is not a required or
            % standard part of phase analysis
            meanLagSec_FileCycle = (circMeanDeg / 360) * fileCyclePeriodSec;

            PLVmatrix(i, j) = plv;
            phaseDiffMatrixDeg(i, j) = circMeanDeg;
            lagSecMatrix(i, j) = meanLagSec_FileCycle;
        end
    end

    [~, nameNoExt, ~] = fileparts(fname);

    % Save ROI-by-ROI matrices for each recording
    PLVtable = array2table(PLVmatrix, ...
        'VariableNames', cellstr(roiNames), ...
        'RowNames', cellstr(roiNames));

    phaseTable = array2table(phaseDiffMatrixDeg, ...
        'VariableNames', cellstr(roiNames), ...
        'RowNames', cellstr(roiNames));

    lagTable = array2table(lagSecMatrix, ...
        'VariableNames', cellstr(roiNames), ...
        'RowNames', cellstr(roiNames));

    writetable(PLVtable, fullfile(outputFolder, [nameNoExt '_PLV_matrix.xlsx']), ...
        'WriteRowNames', true);

    writetable(phaseTable, fullfile(outputFolder, [nameNoExt '_phase_diff_deg_matrix.xlsx']), ...
        'WriteRowNames', true);

    writetable(lagTable, fullfile(outputFolder, [nameNoExt '_lag_sec_matrix.xlsx']), ...
        'WriteRowNames', true);

    % Save one summary row for each ROI pair
    for i = 1:nROIs
        for j = i+1:nROIs

            phaseDiffRad = angle(exp(1i * (roiPhaseRad(:, i) - roiPhaseRad(:, j))));
            phaseDiffRad_use = phaseDiffRad(trimMask);

            circMeanRad = atan2(mean(sin(phaseDiffRad_use)), mean(cos(phaseDiffRad_use)));
            circMeanDeg = rad2deg(circMeanRad);
            absPhaseDeg = abs(circMeanDeg);

            plv = abs(mean(exp(1i * phaseDiffRad_use)));

            if plv > 0
                circSpreadRad = sqrt(-2 * log(plv));
                circSpreadDeg = rad2deg(circSpreadRad);
            else
                circSpreadRad = NaN;
                circSpreadDeg = NaN;
            end

            meanLagSec_FileCycle = (circMeanDeg / 360) * fileCyclePeriodSec;
            meanLagSec_Expected10SecCycle = (circMeanDeg / 360) * expectedCyclePeriodSec;

            newRow = table( ...
                string(fname), durationSec, ...
                string(roiNames(i)), string(roiNames(j)), ...
                lowF, highF, filterOrder, trimEdgeSec, ...
                fileCyclePeriodSec, expectedCyclePeriodSec, ...
                circMeanRad, circMeanDeg, absPhaseDeg, ...
                circSpreadRad, circSpreadDeg, ...
                plv, ...
                meanLagSec_FileCycle, meanLagSec_Expected10SecCycle, ...
                'VariableNames', { ...
                'File','DurationSec', ...
                'ROI_1','ROI_2', ...
                'LowFreqHz','HighFreqHz','FilterOrder','TrimEdgeSec', ...
                'FileCyclePeriodSec','ExpectedCyclePeriodSec', ...
                'CircularMeanRad','CircularMeanDeg','AbsolutePhaseDeg', ...
                'CircularSpreadRad','CircularSpreadDeg', ...
                'PhaseLockingValue', ...
                'MeanLagSec_FileCycle','MeanLagSec_Expected10SecCycle'});

            summaryResults = [summaryResults; newRow];

        end
    end

end

writetable(summaryResults, fullfile(outputFolder, 'vessel_sync_summary.xlsx'));

disp('Finished. Outputs saved to:');
disp(outputFolder);