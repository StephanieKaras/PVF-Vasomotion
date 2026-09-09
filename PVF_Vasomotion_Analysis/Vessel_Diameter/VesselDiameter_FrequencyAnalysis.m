% Mean frequency and bandpower for vessel diameter signals only
% Welch PSD + Hilbert based
%
% For each diameter signal, this analysis estimates vasomotor frequency content using Welch power
% spectral density (PSD) analysis and Hilbert-derived instantaneous frequency.
%
% Welch PSD provides an overall frequency-domain estimate across the recording,
% while the Hilbert transform estimates moment-to-moment instantaneous
% frequency after bandpass filtering.
%
% Input files (.csv) should be in standardized spreadsheet format:
%   column F = vessel diameter signal
%   row 12 onward = numeric data

inputFolder = 'input\directory\here'; % Replace with path to user input data

Fs = 1.951; % Sampling frequency: 1.951 for 2Hz, 1.727 for 1.7Hz, 7.440476 for 7Hz
diameterColumn = 6;
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

   if width(data) < diameterColumn 
       warning('Skipping %s: not enough columns.', files(k).name);
       continue
   end

   diameterRaw = data{:, diameterColumn};
   try
       diameter = rmmissing(double(diameterRaw));
   catch
       warning('Skipping %s: diameter column cannot be converted to numeric.', files(k).name);
       continue
   end

   if isempty(diameter) 
       warning('Skipping %s: diameter data is empty.', files(k).name);
       continue
   end

   % Input diameter data is percent of baseline, where 100 represents baseline
   % Subtracting 100 then converts it to deviation from baseline
   diameterDev = diameter - 100;
   % Detrend to remove slow linear drift from signals before PSD and Hilbert analyses
   diameterDev = detrend(diameterDev);
   % Estimate recording length and choose an adaptive Welch window.
   % Shorter recordings use a shorter window to preserve enough segments
   % for PSD estimation while longer recordings use a 60 sec window
   recordingLengthSec = length(diameterDev) / Fs;

   if recordingLengthSec < 90
       windowLengthSec = 30;
   else
       windowLengthSec = 60;
   end

   windowLength = round(Fs * windowLengthSec);
   overlap = round(windowLength * 0.5);
   nfft = max(256, 2^nextpow2(windowLength));

   % Calculate Welch PSD
   [pxx, f] = pwelch(diameterDev, windowLength, overlap, nfft, Fs);
   % Select PSD values within the vasomotor frequency band
   bandIdx = f >= freqBand(1) & f <= freqBand(2);
   fBand = f(bandIdx);
   pxxBand = pxx(bandIdx);
   % Calculate mean frequency as the PSD-weighted average frequency within band
   meanFreqDiameter = sum(fBand .* pxxBand) / sum(pxxBand);
   % Calculate bandpower from Welch PSD within band
   bpDiameter = bandpower(pxx, f, freqBand, 'psd');
   % Bandpass filter signal for Hilbert transform
   diameterFiltered = bandpass(diameterDev, freqBand, Fs);
   % Calculate analytic signal using Hilbert transform
   analyticSignal = hilbert(diameterFiltered);
   % Amplitude envelope
   ampEnvelope = abs(analyticSignal);
   % Unwrapped phase
   phaseRad = unwrap(angle(analyticSignal));
   % Calculate instantaneous frequency (in Hz)
   instFreq = diff(phaseRad) * Fs / (2*pi);
   % Match envelope length to instantaneous frequency length
   ampEnvelopeForFreq = ampEnvelope(1:end-1);
   % Apply edge trimming
   trimSamples = round(trimSec * Fs);

   if length(instFreq) > 2 * trimSamples
       instFreqTrimmed = instFreq(trimSamples+1:end-trimSamples);
       ampEnvelopeTrimmed = ampEnvelopeForFreq(trimSamples+1:end-trimSamples);
   else
       warning('Skipping Hilbert metrics for %s: recording too short after trimming.', files(k).name);
       instFreqTrimmed = NaN;
       ampEnvelopeTrimmed = NaN;
   end

   % Apply amplitude masking
   if all(isnan(instFreqTrimmed))
       meanInstFreqHilbert = NaN;
       medianInstFreqHilbert = NaN;
       meanHilbertAmpEnvelope = NaN;
       percentKeptHilbert = NaN;
   else
       ampThresh = prctile(ampEnvelopeTrimmed, ampPercentile);
       validIdx = ampEnvelopeTrimmed > ampThresh;

       instFreqMasked = instFreqTrimmed(validIdx);
       ampEnvelopeMasked = ampEnvelopeTrimmed(validIdx);

       meanInstFreqHilbert = mean(instFreqMasked, 'omitnan');
       medianInstFreqHilbert = median(instFreqMasked, 'omitnan');
       meanHilbertAmpEnvelope = mean(ampEnvelopeMasked, 'omitnan');
       percentKeptHilbert = 100 * sum(validIdx) / length(validIdx);
   end

   results = [results; table(string(name), ...
       double(recordingLengthSec), double(windowLengthSec), ...
       double(meanFreqDiameter), double(bpDiameter), ...
       double(meanInstFreqHilbert), double(medianInstFreqHilbert), ...
       double(meanHilbertAmpEnvelope), double(percentKeptHilbert))];
end

results.Properties.VariableNames = {'Filename', ...
   'RecordingLengthSec', 'WelchWindowLengthSec', ...
   'MeanFreq_Diameter_WelchPSD_Hz', 'Bandpower_Diameter_WelchPSD', ...
   'MeanInstFreq_Diameter_Hilbert_Hz', 'MedianInstFreq_Diameter_Hilbert_Hz', ...
   'MeanHilbertAmplitudeEnvelope', 'PercentKept_HilbertAmpMask'};

writetable(results, fullfile(inputFolder, 'vessel_diameter_frequency_output.csv'))