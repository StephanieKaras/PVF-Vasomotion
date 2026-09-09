% Mean frequency and bandpower for cellular calcium and vessel diameter signals
% Welch PSD + Hilbert based
%
% For each paired signal of cellular calcium and vessel diameter oscillations, 
% this analysis estimates frequency content using Welch power
% spectral density (PSD) analysis and Hilbert-derived instantaneous frequency.
%
% Welch PSD provides an overall frequency-domain estimate across the recording,
% while the Hilbert transform estimates moment-to-moment instantaneous
% frequency after bandpass filtering.
%
% Input files (.csv) should be in standardized spreadsheet format:
%   column F = calcium signal
%   column G = vessel diameter signal
%   row 12 onward = numeric data

inputFolder = 'input\directory\here'; % Replace with path to user input data

Fs = 1.727; % Sampling frequency: 1.951 for 2Hz, 1.727 for 1.7Hz, 7.440476 for 7Hz
calciumColumn = 6;
diameterColumn = 7;
freqBand = [0.05 0.2]; % Frequency band (lower and upper limits) used to isolate vasomotor-range oscillations

trimSec = 5; % Edge trimming removes filter/Hilbert edge artifacts
ampPercentile = 20; % Low-amplitude Hilbert phase estimates can be unstable, so phase summaries
% are calculated from time points above this amplitude percentile

results = table();

% Recursively locate all input .csv files in folders and subfolders
files = dir(fullfile(inputFolder, '**', '*.csv'));

for k = 1:length(files)
    inputFile = fullfile(files(k).folder, files(k).name);
    [~, name, ~] = fileparts(files(k).name);

    try
        opts = detectImportOptions(inputFile);
        opts.DataLines = [12, Inf];
        data = readtable(inputFile, opts);
    catch
        warning('Failed to read %s. Skipping.', files(k).name);
        continue
    end

    if width(data) < max(calciumColumn, diameterColumn)
        warning('Skipping %s: not enough columns.', files(k).name);
        continue
    end

    try
        calcium = rmmissing(double(data{:, calciumColumn}));
        diameter = rmmissing(double(data{:, diameterColumn}));
    catch
        warning('Skipping %s: calcium or diameter column cannot be converted to numeric.', files(k).name);
        continue
    end

    if isempty(calcium) || isempty(diameter)
        warning('Skipping %s: missing data.', files(k).name);
        continue
    end

    % Detrend to remove slow linear drift from signals before PSD and Hilbert analyses
    calciumProc = detrend(calcium);

    % Convert diameter from percent baseline to deviation from baseline
    diameterProc = diameter - 100;
    % Detrend to remove slow linear drift from signals before PSD and Hilbert analyses
    diameterProc = detrend(diameterProc);

    % Estimate recording length for each signal
    calciumRecordingLengthSec = length(calciumProc) / Fs;
    diameterRecordingLengthSec = length(diameterProc) / Fs;

    % Choose an adaptive Welch window for calcium.
    % Shorter recordings use a shorter window to preserve enough segments
    % for PSD estimation while longer recordings use a 60 sec window
    if calciumRecordingLengthSec < 90
        calciumWindowLengthSec = 30;
    else
        calciumWindowLengthSec = 60;
    end

    calciumWindowLength = round(Fs * calciumWindowLengthSec);
    calciumWindowLength = min(calciumWindowLength, length(calciumProc));
    calciumWindowLengthSec = calciumWindowLength / Fs;
    calciumOverlap = round(calciumWindowLength * 0.5);
    calciumNFFT = max(256, 2^nextpow2(calciumWindowLength));

    % Choose an adaptive Welch window for diameter using the same rule as
    % above
    if diameterRecordingLengthSec < 90
        diameterWindowLengthSec = 30;
    else
        diameterWindowLengthSec = 60;
    end

    diameterWindowLength = round(Fs * diameterWindowLengthSec);
    diameterWindowLength = min(diameterWindowLength, length(diameterProc));
    diameterWindowLengthSec = diameterWindowLength / Fs;
    diameterOverlap = round(diameterWindowLength * 0.5);
    diameterNFFT = max(256, 2^nextpow2(diameterWindowLength));

    % Calculate Welch PSD for calcium
    [pxxCalcium, fCalcium] = pwelch(calciumProc, calciumWindowLength, calciumOverlap, calciumNFFT, Fs);

    calciumBandIdx = fCalcium >= freqBand(1) & fCalcium <= freqBand(2);
    fCalciumBand = fCalcium(calciumBandIdx);
    pxxCalciumBand = pxxCalcium(calciumBandIdx);

    meanFreqCalcium = sum(fCalciumBand .* pxxCalciumBand) / sum(pxxCalciumBand);
    bpCalcium = bandpower(pxxCalcium, fCalcium, freqBand, 'psd');

    % Calculate Welch PSD for diameter
    [pxxDiameter, fDiameter] = pwelch(diameterProc, diameterWindowLength, diameterOverlap, diameterNFFT, Fs);

    diameterBandIdx = fDiameter >= freqBand(1) & fDiameter <= freqBand(2);
    fDiameterBand = fDiameter(diameterBandIdx);
    pxxDiameterBand = pxxDiameter(diameterBandIdx);

    meanFreqDiameter = sum(fDiameterBand .* pxxDiameterBand) / sum(pxxDiameterBand);
    bpDiameter = bandpower(pxxDiameter, fDiameter, freqBand, 'psd');

    % Hilbert transform for calcium
    calciumFiltered = bandpass(calciumProc, freqBand, Fs);
    analyticCalcium = hilbert(calciumFiltered);
    calciumAmpEnvelope = abs(analyticCalcium);
    % Unwrap phase before differentiating so phase jumps at +/-pi don't
    % create artificial instantaneous frequency spikes
    calciumPhaseRad = unwrap(angle(analyticCalcium));
    % Instantaneous frequency is the time derivative of unwrapped phase
    calciumInstFreq = diff(calciumPhaseRad) * Fs / (2*pi);
    calciumAmpEnvelopeForFreq = calciumAmpEnvelope(1:end-1);

    % Hilbert transform for diameter
    diameterFiltered = bandpass(diameterProc, freqBand, Fs);
    analyticDiameter = hilbert(diameterFiltered);
    diameterAmpEnvelope = abs(analyticDiameter);
    diameterPhaseRad = unwrap(angle(analyticDiameter));
    diameterInstFreq = diff(diameterPhaseRad) * Fs / (2*pi);
    diameterAmpEnvelopeForFreq = diameterAmpEnvelope(1:end-1);

    % Apply edge trimming
    trimSamples = round(trimSec * Fs);

    if length(calciumInstFreq) > 2 * trimSamples
        calciumInstFreqTrimmed = calciumInstFreq(trimSamples+1:end-trimSamples);
        calciumAmpEnvelopeTrimmed = calciumAmpEnvelopeForFreq(trimSamples+1:end-trimSamples);
    else
        warning('Skipping calcium Hilbert metrics for %s: recording too short after trimming.', files(k).name);
        calciumInstFreqTrimmed = NaN;
        calciumAmpEnvelopeTrimmed = NaN;
    end

    if length(diameterInstFreq) > 2 * trimSamples
        diameterInstFreqTrimmed = diameterInstFreq(trimSamples+1:end-trimSamples);
        diameterAmpEnvelopeTrimmed = diameterAmpEnvelopeForFreq(trimSamples+1:end-trimSamples);
    else
        warning('Skipping diameter Hilbert metrics for %s: recording too short after trimming.', files(k).name);
        diameterInstFreqTrimmed = NaN;
        diameterAmpEnvelopeTrimmed = NaN;
    end

    % Apply amplitude masking for calcium
    if all(isnan(calciumInstFreqTrimmed))
        meanInstFreqCalciumHilbert = NaN;
        medianInstFreqCalciumHilbert = NaN;
        meanCalciumHilbertAmpEnvelope = NaN;
        percentKeptCalciumHilbert = NaN;
    else
        calciumAmpThresh = prctile(calciumAmpEnvelopeTrimmed, ampPercentile);
        calciumValidIdx = calciumAmpEnvelopeTrimmed > calciumAmpThresh;

        calciumInstFreqMasked = calciumInstFreqTrimmed(calciumValidIdx);
        calciumAmpEnvelopeMasked = calciumAmpEnvelopeTrimmed(calciumValidIdx);

        meanInstFreqCalciumHilbert = mean(calciumInstFreqMasked, 'omitnan');
        medianInstFreqCalciumHilbert = median(calciumInstFreqMasked, 'omitnan');
        meanCalciumHilbertAmpEnvelope = mean(calciumAmpEnvelopeMasked, 'omitnan');
        percentKeptCalciumHilbert = 100 * sum(calciumValidIdx) / length(calciumValidIdx);
    end

    % Apply amplitude masking for diameter
    if all(isnan(diameterInstFreqTrimmed))
        meanInstFreqDiameterHilbert = NaN;
        medianInstFreqDiameterHilbert = NaN;
        meanDiameterHilbertAmpEnvelope = NaN;
        percentKeptDiameterHilbert = NaN;
    else
        diameterAmpThresh = prctile(diameterAmpEnvelopeTrimmed, ampPercentile);
        diameterValidIdx = diameterAmpEnvelopeTrimmed > diameterAmpThresh;

        diameterInstFreqMasked = diameterInstFreqTrimmed(diameterValidIdx);
        diameterAmpEnvelopeMasked = diameterAmpEnvelopeTrimmed(diameterValidIdx);

        meanInstFreqDiameterHilbert = mean(diameterInstFreqMasked, 'omitnan');
        medianInstFreqDiameterHilbert = median(diameterInstFreqMasked, 'omitnan');
        meanDiameterHilbertAmpEnvelope = mean(diameterAmpEnvelopeMasked, 'omitnan');
        percentKeptDiameterHilbert = 100 * sum(diameterValidIdx) / length(diameterValidIdx);
    end

    results = [results; table(string(name), ...
        double(calciumRecordingLengthSec), double(diameterRecordingLengthSec), ...
        double(calciumWindowLengthSec), double(diameterWindowLengthSec), ...
        double(meanFreqCalcium), double(meanFreqDiameter), ...
        double(bpCalcium), double(bpDiameter), ...
        double(meanInstFreqCalciumHilbert), double(meanInstFreqDiameterHilbert), ...
        double(medianInstFreqCalciumHilbert), double(medianInstFreqDiameterHilbert), ...
        double(meanCalciumHilbertAmpEnvelope), double(meanDiameterHilbertAmpEnvelope), ...
        double(percentKeptCalciumHilbert), double(percentKeptDiameterHilbert))];
end

results.Properties.VariableNames = {'Filename', ...
    'CalciumRecordingLengthSec', 'DiameterRecordingLengthSec', ...
    'CalciumWelchWindowLengthSec', 'DiameterWelchWindowLengthSec', ...
    'MeanFreq_Calcium_WelchPSD_Hz', 'MeanFreq_Diameter_WelchPSD_Hz', ...
    'Bandpower_Calcium_WelchPSD', 'Bandpower_Diameter_WelchPSD', ...
    'MeanInstFreq_Calcium_Hilbert_Hz', 'MeanInstFreq_Diameter_Hilbert_Hz', ...
    'MedianInstFreq_Calcium_Hilbert_Hz', 'MedianInstFreq_Diameter_Hilbert_Hz', ...
    'MeanHilbertAmplitudeEnvelope_Calcium', 'MeanHilbertAmplitudeEnvelope_Diameter', ...
    'PercentKept_Calcium_HilbertAmpMask', 'PercentKept_Diameter_HilbertAmpMask'};

writetable(results, fullfile(inputFolder, 'calcium_vaso_frequency_output.csv'))