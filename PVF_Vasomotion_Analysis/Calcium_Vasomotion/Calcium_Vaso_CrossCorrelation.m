% Cross-correlation analysis for paired cellular calcium and vessel diameter signals
%
% This analysis quantifies temporal coupling between calcium and vessel diameter
% signals using cross-correlation. Signals are mean-centered and cross-correlated
% within a defined maximum lag, with results summarized across symmetric
% (full window) and calcium-leading/diameter-leading ("before" and "after") windows.
%
% Sign convention:
%   xcorr(calcium, diameter)
%   positive lag = diameter leads calcium
%   negative lag = calcium leads diameter
%
% Input files (.csv) should be in standardized spreadsheet format:
%   column F = calcium signal
%   column G = vessel diameter signal
%   row 12 onward = numeric data


inputFolder  = 'input\directory\here'; % Replace with path to user input data
outputFolder = 'crosscorr_output';
calciumCol   = 6;
diameterCol  = 7;
Fs = 1.727; % Sampling frequency: 1.951 for 2Hz, 1.727 for 1.7Hz, 7.440476 for 7Hz
maxLagSec    = 5;  % Sets maximum lag for cross correlation curve in seconds
plotOn       = true;  % Toggle true to create and save plots

winSec = 2;  % Peak-analysis window: computes peak metrics within +/- winSec

% Shuffle-based p-value settings
doPvalue    = true;
nShuffles   = 1000;
minShiftSec = 2*winSec;

% Classification thresholds for coupling direction and sign
ambigMargin = 0.10;
corrThresh  = 0.01;
syncTolSec  = (1/Fs)/2;

rng(0); % Seed rng for reproducible shuffle p-values


% Diameter gating removes traces that are too flat or lack a measurable
% diameter event. Gating is based on deviation-from-baseline units
doGating       = true;
minDiamStd     = 0.02;
minDiamRange   = 0.10;
eventAbsThresh = 5.0;
eventMinCount  = 1;

percentBaselineValue = 100;

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Recursively locate all input .csv files in folders and subfolders
files = dir(fullfile(inputFolder, '**', '*.csv')); 
summaryAll = table(); 

for k = 1:length(files)

    filePath = fullfile(files(k).folder, files(k).name); 
    [~, name, ~] = fileparts(filePath); 

    opts = detectImportOptions(filePath); 
    opts.DataLines = [12 Inf];
    data = readtable(filePath, opts); 

    if width(data) < max(calciumCol, diameterCol)  
        fprintf('Skipping %s: not enough columns\n', files(k).name); 
        continue 
    end

    % Keep only paired samples so calcium and diameter remain time-aligned
    validIdx = ~isnan(data{:, diameterCol}) & ~isnan(data{:, calciumCol});
    diameter_raw = data{validIdx, diameterCol}; 
    calcium_raw  = data{validIdx, calciumCol}; 

    minLen = min(length(calcium_raw), length(diameter_raw)); 
    if minLen == 0 
        fprintf('Skipping %s: no valid paired samples after NaN filtering\n', files(k).name);
        continue  
    end 
    diameter_raw = diameter_raw(1:minLen);
    calcium_raw  = calcium_raw(1:minLen);


    % Convert "percent-of-baseline (100=baseline)" to delta-from-baseline
    diameter_dev = diameter_raw - percentBaselineValue;  % example: 104 -> +4, 96 -> -4

    % Mean-center both signals before cross correlation
    diameter = diameter_dev - mean(diameter_dev,'omitnan');  
    calcium  = calcium_raw  - mean(calcium_raw,'omitnan');  

    % Gating
    GatingApplied = doGating;
    IsRejected = false; 
    RejectReason = "None"; 

    % Flatness metrics are calculated on diameter deviation values
    DiamStd = std(diameter_dev,'omitnan'); 
    DiamRangeRobust = prctile(diameter_dev,95) - prctile(diameter_dev,5);

    if doGating
        if DiamStd < minDiamStd && DiamRangeRobust < minDiamRange
            IsRejected = true;
            RejectReason = "FlatDiameter";
        end
    end

    % Event-based gating (if not rejected):
    % Require at least one diameter event crossing the specified threshold
    EventCount = 0;  
    HasDiameterEvent = true; 

    if doGating && ~IsRejected  
        EventCount = sum(abs(diameter_dev) >= eventAbsThresh); 
        HasDiameterEvent = (EventCount >= eventMinCount); 
        if ~HasDiameterEvent 
            IsRejected = true;  
            RejectReason = "NoDiameterEvent"; 
        end
    end

    ValidCouplingFlag = ~IsRejected; 


    % Compute cross-correlation between calcium and diameter
    maxLagSamples = round(maxLagSec * Fs);
    [xc, lags] = xcorr(calcium, diameter, maxLagSamples, 'coeff');
    lag_times = lags / Fs;

    Tcorr = table(lag_times(:), xc(:), ...  
        'VariableNames', {'Lag_Seconds','Correlation'});
    writetable(Tcorr, fullfile(outputFolder, [name '_CrossCorrelationData.xlsx']));

    % Initialize outputs (always, even rejected files still write summaries)
    peakCorr = NaN; peakLag = NaN; peakAbsCorr = NaN;
    maxCorr = NaN; maxLag = NaN; minCorr = NaN; minLag = NaN;
    PeakLagSide = "NA"; PeakClass = "NA";

    MeanCorr_PeakWin = NaN; MeanAbsCorr_PeakWin = NaN; 
    AUC_Corr_PeakWin = NaN; AUC_AbsCorr_PeakWin = NaN;

    CouplingSign_Peak  = "NA";
    CouplingLabel_Peak = "NA";

    % After-only (Diameter-leading window)
    peakCorr_after = NaN; peakLag_after = NaN;
    maxCorr_after = NaN; minCorr_after = NaN;
    maxLag_after = NaN;  minLag_after = NaN;
    CouplingSign_After  = "NA";
    CouplingLabel_After = "NA";
    MeanCorr_After = NaN;
    MeanAbsCorr_After = NaN;
    AUC_Corr_After = NaN;
    AUC_AbsCorr_After = NaN;
    PeakClass_After = "NA";
    PeakLagClass_After = "NA";

    % Before-only (Calcium-leading window)
    peakCorr_before = NaN; peakLag_before = NaN;  
    maxCorr_before = NaN; minCorr_before = NaN; 
    maxLag_before = NaN;  minLag_before = NaN; 
    CouplingSign_Before  = "NA";  
    CouplingLabel_Before = "NA"; 
    MeanCorr_Before = NaN;  
    MeanAbsCorr_Before = NaN; 
    AUC_Corr_Before = NaN;  
    AUC_AbsCorr_Before = NaN; 
    PeakClass_Before = "NA";  
    PeakLagClass_Before = "NA";

    % Direction-dominance summary 
    % Compares the absolute peak correlation magnitudes in the
    % calcium-leading vs diameter-leading windows and chooses dominant direction
    PeakAbsCorr_After  = NaN;
    PeakAbsCorr_Before = NaN;
    DirectionAbsDiff   = NaN;
    DirectionDominant  = "NA";  % "DiameterLeads" / "CalciumLeads" / "Ambiguous" / "NA"

    % Calculate correlation metrics only for recordings that pass gating
    if ~IsRejected

        % Symmetric peak window centered on zero lag
        inRange = lag_times >= -winSec & lag_times <= winSec;
        xc_in = xc(inRange);
        lags_in = lag_times(inRange);

        if ~isempty(xc_in)

            % Peak correlation is defined as the largest absolute correlation
            % within the analysis window, preserving its sign
            [~, idx] = max(abs(xc_in));
            peakCorr = xc_in(idx);
            peakLag  = lags_in(idx);
            peakAbsCorr = abs(peakCorr);

            % Also store the strongest positive and strongest negative correlations
            [maxCorr, iMax] = max(xc_in);
            [minCorr, iMin] = min(xc_in);
            maxLag = lags_in(iMax);
            minLag = lags_in(iMin);


            if abs(peakLag) <= syncTolSec
                PeakLagSide = "Synchronized";
            elseif peakLag > 0 
                PeakLagSide = "DiameterLeads"; 
            else 
                PeakLagSide = "CalciumLeads";
            end

            if peakCorr > corrThresh 
                PeakClass = "Positive";  
            elseif peakCorr < -corrThresh 
                PeakClass = "Negative"; 
            else 
                PeakClass = "Uncorrelated"; 
            end

            MeanCorr_PeakWin = mean(xc_in,'omitnan');
            MeanAbsCorr_PeakWin = mean(abs(xc_in),'omitnan');
            AUC_Corr_PeakWin = trapz(lags_in, xc_in);
            AUC_AbsCorr_PeakWin = trapz(lags_in, abs(xc_in));

            if peakCorr > 0 
                CouplingSign_Peak = "Pos"; 
            elseif peakCorr < 0 
                CouplingSign_Peak = "Neg"; 
            else  
                CouplingSign_Peak = "Zero"; 
            end

            if abs(abs(maxCorr) - abs(minCorr)) < ambigMargin
                CouplingLabel_Peak = "Ambiguous"; 
            elseif abs(maxCorr) > abs(minCorr)
                CouplingLabel_Peak = "PosDominant";  
            else 
                CouplingLabel_Peak = "NegDominant"; 
            end
        end


        % "After"
        % Positive-lag window: diameter-leading side of xcorr(calcium, diameter)
        inAfter = lag_times >= 0 & lag_times <= winSec;
        xc_after = xc(inAfter);
        lags_after = lag_times(inAfter);

        if ~isempty(xc_after)

            [~, idxA] = max(abs(xc_after));
            peakCorr_after = xc_after(idxA);
            peakLag_after  = lags_after(idxA);
            PeakAbsCorr_After = abs(peakCorr_after);

            [maxCorr_after, iMaxA] = max(xc_after); 
            [minCorr_after, iMinA] = min(xc_after); 
            maxLag_after = lags_after(iMaxA);
            minLag_after = lags_after(iMinA);

            if peakCorr_after > 0 
                CouplingSign_After = "Pos";
            elseif peakCorr_after < 0 
                CouplingSign_After = "Neg";
            else 
                CouplingSign_After = "Zero";
            end

            if abs(abs(maxCorr_after) - abs(minCorr_after)) < ambigMargin 
                CouplingLabel_After = "Ambiguous"; 
            elseif abs(maxCorr_after) > abs(minCorr_after) 
                CouplingLabel_After = "PosDominant"; 
            else 
                CouplingLabel_After = "NegDominant"; 
            end

            MeanCorr_After = mean(xc_after,'omitnan'); 
            MeanAbsCorr_After = mean(abs(xc_after),'omitnan');
            AUC_Corr_After = trapz(lags_after, xc_after); 
            AUC_AbsCorr_After = trapz(lags_after, abs(xc_after)); 

            if peakCorr_after > corrThresh 
                PeakClass_After = "Positive";
            elseif peakCorr_after < -corrThresh 
                PeakClass_After = "Negative"; 
            else 
                PeakClass_After = "Uncorrelated";
            end

            if abs(peakLag_after) <= syncTolSec 
                PeakLagClass_After = "Synchronized";
            else  
                PeakLagClass_After = "DiameterLeads"; 
            end
        end

        % "Before"
        % Negative-lag window: calcium-leading side of xcorr(calcium, diameter)
        inBefore = lag_times >= -winSec & lag_times <= 0;
        xc_before = xc(inBefore);
        lags_before = lag_times(inBefore); 

        if ~isempty(xc_before)

            [~, idxB] = max(abs(xc_before)); 
            peakCorr_before = xc_before(idxB);  
            peakLag_before  = lags_before(idxB);  
            PeakAbsCorr_Before = abs(peakCorr_before); 

            [maxCorr_before, iMaxB] = max(xc_before);  
            [minCorr_before, iMinB] = min(xc_before); 
            maxLag_before = lags_before(iMaxB); 
            minLag_before = lags_before(iMinB);

            if peakCorr_before > 0 
                CouplingSign_Before = "Pos";
            elseif peakCorr_before < 0  
                CouplingSign_Before = "Neg";
            else  
                CouplingSign_Before = "Zero";
            end

            if abs(abs(maxCorr_before) - abs(minCorr_before)) < ambigMargin  
                CouplingLabel_Before = "Ambiguous"; 
            elseif abs(maxCorr_before) > abs(minCorr_before)  
                CouplingLabel_Before = "PosDominant";
            else  
                CouplingLabel_Before = "NegDominant";
            end

            MeanCorr_Before = mean(xc_before,'omitnan');  
            MeanAbsCorr_Before = mean(abs(xc_before),'omitnan');  
            AUC_Corr_Before = trapz(lags_before, xc_before);  
            AUC_AbsCorr_Before = trapz(lags_before, abs(xc_before)); 

            if peakCorr_before > corrThresh  
                PeakClass_Before = "Positive"; 
            elseif peakCorr_before < -corrThresh  
                PeakClass_Before = "Negative"; 
            else  
                PeakClass_Before = "Uncorrelated";
            end

            if abs(peakLag_before) <= syncTolSec 
                PeakLagClass_Before = "Synchronized";
            else 
                PeakLagClass_Before = "CalciumLeads";
            end
        end

        % Direction dominance - compare best r "Before" vs "After"
        if ~isnan(PeakAbsCorr_After) || ~isnan(PeakAbsCorr_Before) 

            if isnan(PeakAbsCorr_Before) && ~isnan(PeakAbsCorr_After) 
                DirectionDominant = "DiameterLeads";  
                DirectionAbsDiff  = NaN;

            elseif isnan(PeakAbsCorr_After) && ~isnan(PeakAbsCorr_Before)  
                DirectionDominant = "CalciumLeads";  
                DirectionAbsDiff  = NaN;

            else 
                DirectionAbsDiff = PeakAbsCorr_After - PeakAbsCorr_Before;  

                if abs(DirectionAbsDiff) < ambigMargin  
                    DirectionDominant = "Ambiguous";
                elseif DirectionAbsDiff > 0  
                    DirectionDominant = "DiameterLeads";  
                else  
                    DirectionDominant = "CalciumLeads";  
                end
            end
        end
    end

    % Plots (always saved if plotOn, even if rejected)
    if plotOn 
        fig = figure('Visible','on');  % use 'off' to speed up batch processing
        hold on; 
        plot(lag_times, xc, 'b-', 'LineWidth', 2); 

        if ~isnan(peakLag) 
            plot(peakLag, peakCorr, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y'); 
        end
        if ~isnan(maxLag)  
            plot(maxLag, maxCorr, 'ro', 'MarkerSize', 8);  
        end
        if ~isnan(minLag) 
            plot(minLag, minCorr, 'go', 'MarkerSize', 8); 
        end

        xline(-winSec, '--k');  
        xline( winSec, '--k');  

        xlabel('Lag (seconds)'); 
        ylabel('Cross-correlation coefficient');
        title(['Cross-correlation: ' name], 'Interpreter', 'none');
        legend('Cross-correlation', 'Peak |r|', 'Max Corr', 'Min Corr', 'Location', 'best');
        grid on;
        hold off;
        saveas(fig, fullfile(outputFolder, [name '_CrossCorrelationPlot_PeakWin.png']));
        close(fig);

        fig = figure('Visible','on'); 
        hold on; 
        plot(lag_times, xc, 'b-', 'LineWidth', 2);

        xline(0, '--k'); 
        xline(winSec, '--k'); 

        if ~isnan(peakLag_after)
            plot(peakLag_after, peakCorr_after, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
            plot(maxLag_after, maxCorr_after, 'ro', 'MarkerSize', 8);  
            plot(minLag_after, minCorr_after, 'go', 'MarkerSize', 8);
        end

        xlabel('Lag (seconds)'); 
        ylabel('Cross-correlation coefficient');
        title(['Cross-correlation (After-only): ' name], 'Interpreter', 'none'); 
        legend('Cross-correlation', 'Peak |r|', 'Max Corr', 'Min Corr', 'Location', 'best');
        grid on; 
        hold off;  
        saveas(fig, fullfile(outputFolder, [name '_CrossCorrelationPlot_AfterOnly.png'])); 
        close(fig); 

        fig = figure('Visible','on'); 
        hold on; 
        plot(lag_times, xc, 'b-', 'LineWidth', 2);

        xline(-winSec, '--k'); 
        xline(0, '--k'); 

        if ~isnan(peakLag_before)
            plot(peakLag_before, peakCorr_before, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');  
            plot(maxLag_before, maxCorr_before, 'ro', 'MarkerSize', 8); 
            plot(minLag_before, minCorr_before, 'go', 'MarkerSize', 8); 
        end 

        xlabel('Lag (seconds)');  
        ylabel('Cross-correlation coefficient'); 
        title(['Cross-correlation (Before-only): ' name], 'Interpreter', 'none');
        legend('Cross-correlation', 'Peak |r|', 'Max Corr', 'Min Corr', 'Location', 'best'); 
        grid on;  
        hold off; 
        saveas(fig, fullfile(outputFolder, [name '_CrossCorrelationPlot_BeforeOnly.png']));
        close(fig); 
    end

    % Permutation test/shuffle p-values to test whether the observed peak correlation exceeds
    % correlations expected after disrupting temporal alignment between
    % signals (random signals)
    peakP = NaN; peakP_after = NaN; peakP_before = NaN; 

    if doPvalue && ~IsRejected && ~isnan(peakCorr)  
        N = length(calcium);  
        minShiftSamples = max(1, round(minShiftSec * Fs));  

        if N > 2*minShiftSamples 
            nullPeak   = zeros(nShuffles,1); 
            nullAfter  = zeros(nShuffles,1); 
            nullBefore = zeros(nShuffles,1);

            for s = 1:nShuffles  % Loop through shuffles
                shift = randi([minShiftSamples, N-minShiftSamples]);  % Random shift amount
                calcium_s = circshift(calcium, shift);  % Circularly shift calcium to break true coupling
                xc_s = xcorr(calcium_s, diameter, maxLagSamples, 'coeff'); 
                lags_s = (-maxLagSamples:maxLagSamples)/Fs; 

                nullPeak(s) = max(abs(xc_s(lags_s >= -winSec & lags_s <= winSec))); 

                if any(lags_s >= 0 & lags_s <= winSec)
                    nullAfter(s) = max(abs(xc_s(lags_s >= 0 & lags_s <= winSec))); 
                else
                    nullAfter(s) = NaN;
                end

                if any(lags_s >= -winSec & lags_s <= 0)
                    nullBefore(s) = max(abs(xc_s(lags_s >= -winSec & lags_s <= 0))); 
                else
                    nullBefore(s) = NaN;
                end
            end

            peakP = (1 + sum(nullPeak >= abs(peakCorr))) / (nShuffles + 1);  

            if ~isnan(peakCorr_after)  
                peakP_after = (1 + sum(nullAfter >= abs(peakCorr_after),'omitnan')) / (nShuffles + 1);
            end

            if ~isnan(peakCorr_before)  
                peakP_before = (1 + sum(nullBefore >= abs(peakCorr_before),'omitnan')) / (nShuffles + 1);
            end
        end
    end

    % Save summary
    summaryT = table( ... 
        peakAbsCorr, peakCorr, peakLag, PeakClass, PeakLagSide, peakP, ... 
        maxCorr, maxLag, minCorr, minLag, ...  
        CouplingSign_Peak, CouplingLabel_Peak, ...  
        MeanCorr_PeakWin, MeanAbsCorr_PeakWin, AUC_Corr_PeakWin, AUC_AbsCorr_PeakWin, ...  
        peakCorr_after, peakLag_after, PeakClass_After, PeakLagClass_After, peakP_after, ...  
        maxCorr_after, maxLag_after, minCorr_after, minLag_after, ...  
        CouplingSign_After, CouplingLabel_After, ... 
        MeanCorr_After, MeanAbsCorr_After, AUC_Corr_After, AUC_AbsCorr_After, ... 
        peakCorr_before, peakLag_before, PeakClass_Before, PeakLagClass_Before, peakP_before, ... 
        maxCorr_before, maxLag_before, minCorr_before, minLag_before, ...  
        CouplingSign_Before, CouplingLabel_Before, ...  
        MeanCorr_Before, MeanAbsCorr_Before, AUC_Corr_Before, AUC_AbsCorr_Before, ...  
        PeakAbsCorr_After, PeakAbsCorr_Before, DirectionAbsDiff, DirectionDominant, ...  
        GatingApplied, ValidCouplingFlag, IsRejected, RejectReason, DiamStd, DiamRangeRobust, EventCount, HasDiameterEvent, ...  
        'VariableNames', { ...  
            'PeakAbsCorrelation', 'PeakCorrelation', 'PeakLag_sec', 'PeakClass', 'PeakLagSide', 'PeakPvalue', ...
            'MaxCorrelation', 'MaxLag_sec', 'MinCorrelation', 'MinLag_sec', ...
            'CouplingSign_Peak', 'CouplingLabel_Peak', ...
            'MeanCorrelation_PeakWin', 'MeanAbsCorrelation_PeakWin', 'AUC_Corr_PeakWin', 'AUC_AbsCorr_PeakWin', ...
            'PeakCorrelation_After', 'PeakLag_sec_After', 'PeakClass_After', 'PeakLagClass_After', 'PeakPvalue_After', ...
            'MaxCorrelation_After', 'MaxLag_sec_After', 'MinCorrelation_After', 'MinLag_sec_After', ...
            'CouplingSign_After', 'CouplingLabel_After', ...
            'MeanCorrelation_After', 'MeanAbsCorrelation_After', 'AUC_Corr_After', 'AUC_AbsCorr_After', ...
            'PeakCorrelation_Before', 'PeakLag_sec_Before', 'PeakClass_Before', 'PeakLagClass_Before', 'PeakPvalue_Before', ...
            'MaxCorrelation_Before', 'MaxLag_sec_Before', 'MinCorrelation_Before', 'MinLag_sec_Before', ...
            'CouplingSign_Before', 'CouplingLabel_Before', ...
            'MeanCorrelation_Before', 'MeanAbsCorrelation_Before', 'AUC_Corr_Before', 'AUC_AbsCorr_Before', ...
            'PeakAbsCorr_After', 'PeakAbsCorr_Before', 'DirectionAbsDiff', 'DirectionDominant', ...
            'GatingApplied', 'ValidCouplingFlag', 'IsRejected', 'RejectReason', 'DiameterStd', 'DiameterRangeRobust', 'DiameterEventCount', 'HasDiameterEvent' ...
        } ... 
    );  

    summaryT.Filename = {name};  
    summaryAll = [summaryAll; summaryT]; 

    writetable(summaryT, fullfile(outputFolder, [name '_CorrelationSummary.xlsx'])); 
end

if ~isempty(summaryAll) && any(strcmp(summaryAll.Properties.VariableNames, 'Filename')) 
    summaryAll = movevars(summaryAll,'Filename','Before',1); 
end

if ~isempty(summaryAll) && width(summaryAll) > 0  
    writetable(summaryAll, fullfile(outputFolder,'Combined_CorrelationSummary.xlsx')); 
else  
    warning('No summaries were produced. Combined_CorrelationSummary.xlsx not written.'); 
end

% Merge correlation columns
ccFiles = dir(fullfile(outputFolder,'*_CrossCorrelationData.xlsx'));

if isempty(ccFiles)
    warning('No *_CrossCorrelationData.xlsx files found in outputFolder: %s. Combined_CorrelationColumns.xlsx not written.', outputFolder); 
else  
    mergedCorr = table(); 

    for i = 1:length(ccFiles)  
        T = readtable(fullfile(ccFiles(i).folder, ccFiles(i).name)); 

        if i == 1 
            mergedCorr = table(T.Lag_Seconds,'VariableNames',{'Lag_Seconds'}); 
        else 
            if height(T) ~= height(mergedCorr)  
                warning('Skipping %s: row count mismatch.', ccFiles(i).name);
                continue  
            end
        end

        [~, fname, ~] = fileparts(ccFiles(i).name);
        varName = matlab.lang.makeValidName(regexprep(fname,'_CrossCorrelationData','')); 
        mergedCorr.(varName) = T.Correlation;  
    end

    writetable(mergedCorr, fullfile(outputFolder,'Combined_CorrelationColumns.xlsx')); 
end
