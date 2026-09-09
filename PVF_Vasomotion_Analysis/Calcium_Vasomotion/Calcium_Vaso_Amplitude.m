% Maximum, minimum, and amplitude of cellular calcium and vessel diameter signals
%
% For both calcium and vessel diameter, this analysis identifies the maximum/minimum values
% and amplitude calculated as maximum minus minimum.
%
% Input files (.csv) should be in standardized spreadsheet format:
%   column F = calcium signal
%   column G = vessel diameter signal
%   row 12 onward = numeric data

inputFolder = 'input\directory\here'; % Replace with path to user input data
dataStartRow = 12;
numHeaderLines = dataStartRow - 1;

calciumCol = 6;
diameterCol = 7;

outXLSX = fullfile(inputFolder, 'amplitude_output.xlsx');

% Recursively locate all input .csv files in folders and subfolders
files = dir(fullfile(inputFolder, '**', '*.csv')); 
if isempty(files)
    error('No CSV files found in: %s', inputFolder);
end

% Preallocate an empty results table with one output row per recording
results = table( ...
    strings(0,1), ...
    nan(0,1), nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    nan(0,1), nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    'VariableNames', { ...
    'File', ...
    'CalciumMaxValue','InputRowOfCalciumMax','CalciumMinValue','InputRowOfCalciumMin','CalciumAmplitude', ...
    'DiameterMaxValue','InputRowOfDiameterMax','DiameterMinValue','InputRowOfDiameterMin','DiameterAmplitude'} ...
);

for k = 1:numel(files)
    fname = files(k).name;
    fpath = fullfile(files(k).folder, fname);

    M = readmatrix(fpath, 'NumHeaderLines', numHeaderLines);

    if size(M,2) < max(calciumCol, diameterCol)
        continue;
    end

    cal = M(:, calciumCol);
    diamRaw = M(:, diameterCol);
    % Input diameter data is percent of baseline, where 100 represents baseline
    % Subtracting 100 then converts it to deviation from baseline
    diam = diamRaw - 100;

    % Extract calcium maximum, minimum, and peak-to-trough amplitude
    if all(isnan(cal)) || isempty(cal)
        calMax = NaN;
        inputRowOfCalciumMax = NaN;
        calMin = NaN;
        inputRowOfCalciumMin = NaN;
        calAmp = NaN;
    else
        [calMax, idxCalMax] = max(cal);
        inputRowOfCalciumMax = idxCalMax + numHeaderLines;

        [calMin, idxCalMin] = min(cal);
        inputRowOfCalciumMin = idxCalMin + numHeaderLines;

        calAmp = calMax - calMin;
    end

    % Extract diameter maximum, minimum, and peak-to-trough amplitude
    if all(isnan(diam)) || isempty(diam)
        diamMax = NaN;
        inputRowOfDiameterMax = NaN;
        diamMin = NaN;
        inputRowOfDiameterMin = NaN;
        diamAmp = NaN;
    else
        [diamMax, idxDiamMax] = max(diam);
        inputRowOfDiameterMax = idxDiamMax + numHeaderLines;

        [diamMin, idxDiamMin] = min(diam);
        inputRowOfDiameterMin = idxDiamMin + numHeaderLines;

        diamAmp = diamMax - diamMin;
    end

    results = [results; table( ...
        string(fname), ...
        calMax, inputRowOfCalciumMax, calMin, inputRowOfCalciumMin, calAmp, ...
        diamMax, inputRowOfDiameterMax, diamMin, inputRowOfDiameterMin, diamAmp, ...
        'VariableNames', results.Properties.VariableNames)];
end

writetable(results, outXLSX);
fprintf('Finished. Outputs saved to: %s\n', outXLSX);