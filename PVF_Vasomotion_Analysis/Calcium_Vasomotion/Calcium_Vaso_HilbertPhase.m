% Hilbert phase analysis for cellular calcium and vessel diameter signals
%
% This analysis quantifies the phase relationship between cellular calcium and
% vessel diameter oscillations. Signals are detrended, bandpass filtered in
% the vasomotor frequency range, then transformed using the Hilbert
% transform to calculate instantaneous phase and the phase difference between
% calcium and vessel diameter oscillations.
%
% For each recording, the script calculates circular mean phase, phase-locking
% value (PLV), Rayleigh test statistics, and von Mises concentration.
%
% Input files (.csv) should be in standardized spreadsheet format:
%   column E = time
%   column F = calcium signal
%   column G = vessel diameter signal
%   row 12 onward = numeric data

inputFolder = 'input\directory\here'; % Replace with path to user input data
outputFolder = fullfile(inputFolder, 'hilbert_phase_output');

% Standardized input file format
dataStartRow = 12;
numHeaderLines = dataStartRow - 1;

timeCol = 5;
calciumCol = 6;
diameterCol = 7;

Fs = 1.951; % Sampling frequency: 1.951 for 2Hz, 1.727 for 1.7Hz, 7.440476 for 7Hz
lowF = 0.025; % Frequency band (lower and upper limits) used to isolate vasomotor-range oscillations
highF = 0.25;
filterOrder = 2;
trimEdgeSec = 5; % Edge trimming removes filter/Hilbert edge artifacts

% Low-amplitude Hilbert phase estimates can be unstable, so phase summaries
% are calculated from time points above this amplitude percentile
ampPercentile = 20;

% PLV threshold used to label phase estimates as stable enough for interpretation
minPLVForStablePhase = 0.10;
expectedCyclePeriodSec = 10; % Expected vasomotor cycle period used as a standardized lag estimate
minPeakDistanceSec = 3; % Minimum distance between detected diameter peaks for estimating each file-specific cycle period

histBinEdgesDeg = -180:20:180;

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

summaryResults = table();

% Recursively locate all input .csv files in folders and subfolders
files = dir(fullfile(inputFolder, '**','*.csv'));
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

    if numel(t) < 20
        fprintf('Skipping %s: too few valid samples.\n', fname);
        continue;
    end

    durationSec = t(end) - t(1);

    % Input diameter data is percent of baseline, where 100 represents baseline
    % Subtracting 100 then converts it to deviation from baseline
    diameter_dev = diameter - 100;

    % Detrend to remove slow linear drift before phase estimation
    calcium_dt = detrend(calcium);
    diameter_dt = detrend(diameter_dev);

    % Zero-phase bandpass filtering avoids introducing artificial phase lags
    calcium_f = filtfilt(b, a, calcium_dt);
    diameter_f = filtfilt(b, a, diameter_dt);

    % Estimate a file-specific diameter cycle period from filtered diameter peaks
    % This is used to convert phase angle into an approximate timing lag
    minPeakDistanceSamples = round(minPeakDistanceSec * Fs);

    [~, locs] = findpeaks(diameter_f, ...
        'MinPeakDistance', minPeakDistanceSamples);

    if numel(locs) >= 2
        peakTimes = t(locs);
        cyclePeriods = diff(peakTimes);
        fileCyclePeriodSec = median(cyclePeriods, 'omitnan');
    else
        fileCyclePeriodSec = NaN;
    end

    % Hilbert transform is applied to generate the analytic signal, from which the
    % amplitude envelope and instantaneous phase are derived
    calcium_h = hilbert(calcium_f);
    diameter_h = hilbert(diameter_f);

    calciumEnvelope = abs(calcium_h);
    diameterEnvelope = abs(diameter_h);

    % Define amplitude thresholds separately for calcium and diameter
    ampThreshCalcium = prctile(calciumEnvelope, ampPercentile);
    ampThreshDiameter = prctile(diameterEnvelope, ampPercentile);

    % Keep time points where both signals have sufficient oscillatory amplitude
    ampMask = calciumEnvelope > ampThreshCalcium & ...
              diameterEnvelope > ampThreshDiameter;

    phaseCalciumRad = angle(calcium_h);
    phaseDiameterRad = angle(diameter_h);

    phaseCalciumDeg = rad2deg(phaseCalciumRad);
    phaseDiameterDeg = rad2deg(phaseDiameterRad);

    % Phase difference is wrapped to the circular range -pi to pi, meaning
    % that positive values indicate calcium phase is ahead of diameter phase
    % according to the subtraction calcium - diameter
    phaseDiffRad = angle(exp(1i * (phaseCalciumRad - phaseDiameterRad)));
    phaseDiffDeg = rad2deg(phaseDiffRad);

    % Exclude set amount of data at the beginning and end of the trace where filtering and Hilbert
    % transform estimates are prone to edge artifacts
    trimMask = t >= (t(1) + trimEdgeSec) & t <= (t(end) - trimEdgeSec);

    if sum(trimMask) < 10
        trimMask = true(size(t));
    end

    summaryMask = trimMask & ampMask;

    if sum(summaryMask) < 10
        summaryMask = trimMask;
    end

    percentKept = 100 * sum(summaryMask) / numel(summaryMask);

    phaseDiffRad_trim = phaseDiffRad(summaryMask);
    phaseDiffDeg_trim = phaseDiffDeg(summaryMask);

    % Circular mean gives the preferred phase offset between calcium and diameter
    circMeanRad = atan2(mean(sin(phaseDiffRad_trim)), mean(cos(phaseDiffRad_trim)));
    circMeanDeg = rad2deg(circMeanRad);

    % PLV is the length of the mean resultant vector of phase differences.
    % Values near 1 indicate consistent phase locking while values near 0 indicate
    % broad distribution of phase differences
    plv = abs(mean(exp(1i * phaseDiffRad_trim)));

    % Circular statistics calculated from the same phase-difference vector.
    % Rayleigh statistics test for non-uniform phase clustering, while kappa
    % estimates concentration of the phase distribution
    [rayleighP, rayleighZ] = localRayleighTest(phaseDiffRad_trim);
    vonMisesKappa = localEstimateKappa(phaseDiffRad_trim);

    phaseStable = plv >= minPLVForStablePhase;

    % Convert mean phase offset into time using both the file-specific cycle
    % period and a standardized expected 10 s vasomotor cycle. This time-based measure is
    % included for exploratory interpretation and is not a required or
    % standard part of phase analysis
    meanLagSec_FileCycle = (circMeanDeg / 360) * fileCyclePeriodSec;
    meanLagSec_Expected10SecCycle = (circMeanDeg / 360) * expectedCyclePeriodSec;

    % Circular spread provides an interpretable angular dispersion estimate
    if plv > 0
        circSpreadRad = sqrt(-2 * log(plv));
        circSpreadDeg = rad2deg(circSpreadRad);
    else
        circSpreadRad = NaN;
        circSpreadDeg = NaN;
    end

    absPhaseDeg = abs(circMeanDeg);

    % Categorize the phase relationship relative to diameter peak/trough
    if absPhaseDeg <= 45
        phaseRelationship = "CalciumNearDiameterPeak";
    elseif absPhaseDeg >= 135
        phaseRelationship = "CalciumNearDiameterTrough";
    elseif circMeanDeg > 0
        phaseRelationship = "CalciumPositivePhaseOffset";
    else
        phaseRelationship = "CalciumNegativePhaseOffset";
    end

    % Welch PSD is calculated for quality control (qc) and visualization of frequency content
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

    [~, nameNoExt, ~] = fileparts(fname);

    % Save time-resolved intermediate values so phase calculations
    % can be inspected for individual recordings
    timeSeriesTable = table( ...
        t, calcium, diameter, diameter_dev, ...
        calcium_dt, diameter_dt, ...
        calcium_f, diameter_f, ...
        calciumEnvelope, diameterEnvelope, ...
        ampThreshCalcium * ones(size(t)), ampThreshDiameter * ones(size(t)), ...
        ampMask, trimMask, summaryMask, ...
        phaseCalciumRad, phaseCalciumDeg, ...
        phaseDiameterRad, phaseDiameterDeg, ...
        phaseDiffRad, phaseDiffDeg, ...
        'VariableNames', { ...
        'TimeSec', ...
        'RawCalcium','RawDiameter','DiameterDevFrom100', ...
        'DetrendedCalcium','DetrendedDiameter', ...
        'FilteredCalcium','FilteredDiameter', ...
        'CalciumAmplitudeEnvelope','DiameterAmplitudeEnvelope', ...
        'CalciumEnvelopeThreshold','DiameterEnvelopeThreshold', ...
        'AmplitudeMask','TrimMask','IncludedForSummary', ...
        'PhaseCalciumRad','PhaseCalciumDeg', ...
        'PhaseDiameterRad','PhaseDiameterDeg', ...
        'PhaseDiffRad','PhaseDiffDeg'});

    writetable(timeSeriesTable, fullfile(outputFolder, [nameNoExt '_hilbert_timeseries.xlsx']));

    % Save histogram counts for the included phase-difference samples
    [histCounts, histEdges] = histcounts(phaseDiffDeg_trim, histBinEdgesDeg);
    histCenters = histEdges(1:end-1) + diff(histEdges)/2;

    histTable = table(histCenters(:), histCounts(:), ...
        'VariableNames', {'PhaseDiffBinCenterDeg','Count'});

    writetable(histTable, fullfile(outputFolder, [nameNoExt '_phase_histogram.xlsx']));

    % Plot filtered calcium and diameter traces
    fig1 = figure('Visible','off');

    ax1 = subplot(2,1,1);
    plot(t, calcium_f, 'LineWidth', 1.2);
    ylabel('Filtered calcium');
    title(['Filtered Signals: ' nameNoExt], 'Interpreter','none');

    ax2 = subplot(2,1,2);
    plot(t, diameter_f, 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Filtered diameter');

    linkaxes([ax1 ax2], 'x');

    saveas(fig1, fullfile(outputFolder, [nameNoExt '_filtered_signals.png']));
    print(fig1, fullfile(outputFolder, [nameNoExt '_filtered_signals.eps']), '-depsc');
    close(fig1);

    % Plot phase difference across time with circular mean overlaid
    fig2 = figure('Visible','off');
    plot(t, phaseDiffDeg, 'LineWidth', 1.2); hold on;
    yline(circMeanDeg, '--');
    xlabel('Time (s)');
    ylabel('Phase difference (deg)');
    title(['Phase Difference Over Time: ' nameNoExt], 'Interpreter','none');
    ylim([-180 180]);

    saveas(fig2, fullfile(outputFolder, [nameNoExt '_phase_difference_timeseries.png']));
    print(fig2, fullfile(outputFolder, [nameNoExt '_phase_difference_timeseries.eps']), '-depsc');
    close(fig2);

    % Plot linear histogram of included phase differences
    fig3 = figure('Visible','off');
    histogram(phaseDiffDeg_trim, histBinEdgesDeg); hold on;
    xline(circMeanDeg, '--');
    xlabel('Phase difference (deg)');
    ylabel('Count');
    title(['Phase Difference Histogram: ' nameNoExt], 'Interpreter','none');
    xlim([-180 180]);

    saveas(fig3, fullfile(outputFolder, [nameNoExt '_phase_histogram.png']));
    print(fig3, fullfile(outputFolder, [nameNoExt '_phase_histogram.eps']), '-depsc');
    close(fig3);

    % Plot polar histogram to visualize the circular distribution of phase offsets
    figPolar = figure('Visible','off');
    polarhistogram(phaseDiffRad_trim, deg2rad(histBinEdgesDeg), 'Normalization','probability');
    title(['Polar Phase Histogram: ' nameNoExt], 'Interpreter','none');

    saveas(figPolar, fullfile(outputFolder, [nameNoExt '_polar_phase_histogram.png']));
    print(figPolar, fullfile(outputFolder, [nameNoExt '_polar_phase_histogram.eps']), '-depsc');
    close(figPolar);

    % Save PSD plots for visual inspection of frequency content
    fig4 = figure('Visible','off');

    ax3 = subplot(2,1,1);
    plot(Fcal, PxxCal, 'LineWidth', 1.2); hold on;
    xline(lowF, '--');
    xline(highF, '--');
    ylabel('Calcium PSD');
    title(['PSD Only: ' nameNoExt], 'Interpreter','none');
    xlim([0 max(0.5, highF + 0.05)]);

    ax4 = subplot(2,1,2);
    plot(Fdiam, PxxDiam, 'LineWidth', 1.2); hold on;
    xline(lowF, '--');
    xline(highF, '--');
    xlabel('Frequency (Hz)');
    ylabel('Diameter PSD');
    xlim([0 max(0.5, highF + 0.05)]);

    linkaxes([ax3 ax4], 'x');

    saveas(fig4, fullfile(outputFolder, [nameNoExt '_psd.png']));
    print(fig4, fullfile(outputFolder, [nameNoExt '_psd.eps']), '-depsc');
    close(fig4);

    % Store one summary row per recording
    newRow = table( ...
        string(fname), durationSec, ...
        lowF, highF, filterOrder, trimEdgeSec, ampPercentile, percentKept, ...
        minPLVForStablePhase, expectedCyclePeriodSec, minPeakDistanceSec, fileCyclePeriodSec, ...
        circMeanRad, circMeanDeg, absPhaseDeg, ...
        circSpreadRad, circSpreadDeg, ...
        plv, rayleighP, rayleighZ, vonMisesKappa, phaseStable, ...
        meanLagSec_FileCycle, meanLagSec_Expected10SecCycle, string(phaseRelationship), ...
        'VariableNames', { ...
        'File','DurationSec', ...
        'LowFreqHz','HighFreqHz','FilterOrder','TrimEdgeSec','AmplitudePercentileThreshold','PercentDataKept', ...
        'MinPLVForStablePhase','ExpectedCyclePeriodSec','MinPeakDistanceSec','FileCyclePeriodSec', ...
        'CircularMeanRad','CircularMeanDeg','AbsolutePhaseDeg', ...
        'CircularSpreadRad','CircularSpreadDeg', ...
        'PhaseLockingValue','RayleighP','RayleighZ','VonMisesKappa','PhaseStable', ...
        'MeanLagSec_FileCycle','MeanLagSec_Expected10SecCycle','PhaseRelationship'});

    summaryResults = [summaryResults; newRow];

end

writetable(summaryResults, fullfile(outputFolder, 'hilbert_summary.xlsx'));

disp('Finished. Outputs saved to:');
disp(outputFolder);

% Local helper functions used to calculate circular statistics from the phase
% difference distribution
%
% Rayleigh test for non-uniformity of circular data.
% The test evaluates whether phase angles are uniformly distributed
    % around the circle or show evidence of a preferred direction
function [p, z] = localRayleighTest(alpha)
    alpha = alpha(isfinite(alpha));
    n = numel(alpha);

    if n == 0
        p = NaN;
        z = NaN;
        return;
    end

    R = abs(sum(exp(1i*alpha)));
    Rbar = R / n;
    z = n * Rbar^2;

    p = exp(sqrt(1 + 4*n + 4*(n^2 - R^2)) - (1 + 2*n));
    p = min(max(p, 0), 1);
end

% Estimate von Mises kappa for circular concentration.
% Higher kappa values indicate a tighter phase distribution around the
% mean direction while lower values indicate broader phase dispersion
function kappa = localEstimateKappa(alphaOrRbar)
    if numel(alphaOrRbar) > 1
        alpha = alphaOrRbar(isfinite(alphaOrRbar));
        n = numel(alpha);

        if n == 0
            kappa = NaN;
            return;
        end

        Rbar = abs(mean(exp(1i*alpha)));
    else
        Rbar = alphaOrRbar;
        n = NaN;
    end

    if ~isfinite(Rbar) || Rbar < 0
        kappa = NaN;
    elseif Rbar < 0.53
        kappa = 2*Rbar + Rbar^3 + (5*Rbar^5)/6;
    elseif Rbar < 0.85
        kappa = -0.4 + 1.39*Rbar + 0.43/(1 - Rbar);
    else
        kappa = 1/(Rbar^3 - 4*Rbar^2 + 3*Rbar);
    end

    % Small-sample correction
    if isfinite(n) && n < 15 && n > 1
        if kappa < 2
            kappa = max(kappa - 2/(n*kappa), 0);
        else
            kappa = ((n - 1)^3 * kappa) / (n^3 + n);
        end
    end
    % Cap very large estimates to avoid unstable output
    kappa = min(max(kappa, 0), 100);
end
