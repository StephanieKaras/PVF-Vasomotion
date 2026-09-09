% Circular Statistical Tests: Watson-Williams and Mardia-Watson-Wheeler
%
% Watson-Williams: tests for differences in mean circular direction under
% parametric assumptions (such as sufficient angular concentration).
%
% Mardia-Watson-Wheeler: nonparametric comparison of overall circular
% distributions and can detect differences such as shape, 
% spread, and direction.
%
% Both overall and pairwise group comparisons are calculated. Pairwise
% p-values are corrected for multiple comparisons using the Holm method.
%
% Input files (.xlsx) should be in standardized spreadsheet format in which 
% each column represents one group/condition and each row contains one observation per column.

inputFolder = 'input\directory\here'; % Replace with path to user input data
outputFolder = fullfile(inputFolder, 'WW_MWW_Output');

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Locate input Excel files
files = dir(fullfile(inputFolder, '**', '*.xlsx'));
files = files(~startsWith({files.name}, '~$'));

if isempty(files)
    error('No Excel files found in: %s', inputFolder);
end

for f = 1:length(files)

    filename = files(f).name;
    filepath = fullfile(files(f).folder, filename);

    fprintf('\n============================================\n');
    fprintf('Processing: %s\n', filename);
    fprintf('============================================\n');

    T = readtable(filepath, ...
        'Sheet', 1, ...
        'VariableNamingRule', 'preserve');

    groupNames = T.Properties.VariableNames;
    numGroups = width(T);

    angles = cell(numGroups, 1);
    n = zeros(numGroups, 1);

    % Convert each group to numeric angles in radians
    for g = 1:numGroups

        x = T{:, g};

        if ~isnumeric(x)
            x = str2double(string(x));
        end

        x = x(~isnan(x));

        angles{g} = deg2rad(x(:));
        n(g) = length(x);
    end

    keep = n > 0;

    angles = angles(keep);
    groupNames = groupNames(keep);
    n = n(keep);

    numGroups = length(angles);

    if numGroups < 2
        warning('Skipping %s because fewer than two groups were found.', filename);
        continue
    end

    groupTable = table( ...
        string(groupNames(:)), ...
        n, ...
        'VariableNames', {'Group', 'N'});

    % Watson-Williams test for differences in mean circular direction
    % with parametric assumptions (kappa)
    [F_WW, df1_WW, df2_WW, p_WW, kappa_WW] = ...
        watsonWilliams(angles);

    % Low kappa (concentration) indicates that Watson-Williams assumptions may not be met
    if kappa_WW < 1
        WWwarning = ...
            "Kappa < 1. Watson-Williams assumptions may not be met.";
    else
        WWwarning = "";
    end

    overallWW = table( ...
        F_WW, ...
        df1_WW, ...
        df2_WW, ...
        p_WW, ...
        kappa_WW, ...
        WWwarning, ...
        'VariableNames', ...
        {'F', 'DF1', 'DF2', 'PValue', 'EstimatedKappa', 'Warning'});

    % Mardia-Watson-Wheeler for differences in overall circular distributions
    % without parametric assumptions
    [W_MWW, df_MWW, p_MWW] = ...
        mardiaWatsonWheeler(angles);

    % Flag small ample sizes where the chi-square approximation may be unreliable
    if any(n < 10)
        MWWwarning = ...
            "At least one group has fewer than 10 observations. Chi-square approximation may be unreliable.";
    else
        MWWwarning = "";
    end

    overallMWW = table( ...
        W_MWW, ...
        df_MWW, ...
        p_MWW, ...
        MWWwarning, ...
        'VariableNames', ...
        {'W', 'DF', 'PValue', 'Warning'});

    % Generate all pairwise group comparisons
    pairs = nchoosek(1:numGroups, 2);
    numPairs = size(pairs, 1);

    Group1 = strings(numPairs, 1);
    Group2 = strings(numPairs, 1);

    N1 = zeros(numPairs, 1);
    N2 = zeros(numPairs, 1);

    WW_F = zeros(numPairs, 1);
    WW_DF1 = zeros(numPairs, 1);
    WW_DF2 = zeros(numPairs, 1);
    WW_P = zeros(numPairs, 1);
    WW_Kappa = zeros(numPairs, 1);
    WW_Warning = strings(numPairs, 1);

    MWW_W = zeros(numPairs, 1);
    MWW_DF = zeros(numPairs, 1);
    MWW_P = zeros(numPairs, 1);
    MWW_Warning = strings(numPairs, 1);

    % Run both circular tests for each pair of groups
    for p = 1:numPairs

        g1 = pairs(p, 1);
        g2 = pairs(p, 2);

        Group1(p) = string(groupNames{g1});
        Group2(p) = string(groupNames{g2});

        N1(p) = n(g1);
        N2(p) = n(g2);

        pairAngles = {
            angles{g1}
            angles{g2}
        };

        [WW_F(p), ...
            WW_DF1(p), ...
            WW_DF2(p), ...
            WW_P(p), ...
            WW_Kappa(p)] = ...
            watsonWilliams(pairAngles);

        if WW_Kappa(p) < 1
            WW_Warning(p) = ...
                "Kappa < 1. Watson-Williams assumptions may not be met.";
        end

        [MWW_W(p), ...
            MWW_DF(p), ...
            MWW_P(p)] = ...
            mardiaWatsonWheeler(pairAngles);

        if N1(p) < 10 || N2(p) < 10
            MWW_Warning(p) = ...
                "At least one group has fewer than 10 observations. Chi-square approximation may be unreliable.";
        end
    end

    % Apply Holm correction separately to pairwise p-values from each test
    WW_P_Holm = holmCorrection(WW_P);
    MWW_P_Holm = holmCorrection(MWW_P);

    % Combine and store pairwise results
    pairwiseTable = table( ...
        Group1, ...
        Group2, ...
        N1, ...
        N2, ...
        WW_F, ...
        WW_DF1, ...
        WW_DF2, ...
        WW_P, ...
        WW_P_Holm, ...
        WW_Kappa, ...
        WW_Warning, ...
        MWW_W, ...
        MWW_DF, ...
        MWW_P, ...
        MWW_P_Holm, ...
        MWW_Warning, ...
        'VariableNames', { ...
        'Group1', ...
        'Group2', ...
        'N1', ...
        'N2', ...
        'WW_F', ...
        'WW_DF1', ...
        'WW_DF2', ...
        'WW_P_Raw', ...
        'WW_P_Holm', ...
        'WW_Kappa', ...
        'WW_Warning', ...
        'MWW_W', ...
        'MWW_DF', ...
        'MWW_P_Raw', ...
        'MWW_P_Holm', ...
        'MWW_Warning'});

    % Save group information, overall tests, and pairwise results
    [~, baseName, ~] = fileparts(filename);

    outputFile = fullfile( ...
        outputFolder, ...
        [baseName '_WW_MWW_results.xlsx']);

    if exist(outputFile, 'file')
        delete(outputFile);
    end

    writetable(groupTable, ...
        outputFile, ...
        'Sheet', 'Groups');

    writetable(overallWW, ...
        outputFile, ...
        'Sheet', 'Overall_WW');

    writetable(overallMWW, ...
        outputFile, ...
        'Sheet', 'Overall_MWW');

    writetable(pairwiseTable, ...
        outputFile, ...
        'Sheet', 'Pairwise');

    % Display results in the MATLAB command window
    fprintf('\nGroups:\n');

    for g = 1:numGroups
        fprintf('%s: n = %d\n', groupNames{g}, n(g));
    end

    fprintf('\nOverall Watson-Williams:\n');
    fprintf('F = %.6f\n', F_WW);
    fprintf('df = (%d, %d)\n', df1_WW, df2_WW);
    fprintf('p = %.6g\n', p_WW);
    fprintf('kappa = %.6f\n', kappa_WW);

    if strlength(WWwarning) > 0
        fprintf('Warning: %s\n', WWwarning);
    end

    fprintf('\nOverall Mardia-Watson-Wheeler:\n');
    fprintf('W = %.6f\n', W_MWW);
    fprintf('df = %d\n', df_MWW);
    fprintf('p = %.6g\n', p_MWW);

    if strlength(MWWwarning) > 0
        fprintf('Warning: %s\n', MWWwarning);
    end

    fprintf('\nPairwise comparisons:\n');

    for p = 1:numPairs

        fprintf('\n%s vs %s\n', Group1(p), Group2(p));

        fprintf('  WW:  F = %.6f, raw p = %.6g, Holm p = %.6g\n', ...
            WW_F(p), WW_P(p), WW_P_Holm(p));

        fprintf('  MWW: W = %.6f, raw p = %.6g, Holm p = %.6g\n', ...
            MWW_W(p), MWW_P(p), MWW_P_Holm(p));
    end

    fprintf('\nSaved: %s\n', outputFile);
end

fprintf('\n============================================\n');
fprintf('Finished.\n');
fprintf('============================================\n');


% Functions
%
%
% Watson-Williams test:
% Compares mean circular directions across groups and returns the F statistic,
% degrees of freedom, p-value, and estimated kappa

function [F, df1, df2, p, kappa] = watsonWilliams(angles)

    k = length(angles);
    n = cellfun(@length, angles);
    N = sum(n);

    Rgroup = zeros(k, 1);

    % Calculate the circular vector length for each group
    for g = 1:k

        theta = angles{g};

        C = sum(cos(theta));
        S = sum(sin(theta));

        Rgroup(g) = sqrt(C^2 + S^2);
    end

    % Calculate the circular vector length for all observations combined
    allAngles = vertcat(angles{:});

    Ctotal = sum(cos(allAngles));
    Stotal = sum(sin(allAngles));

    Rtotal = sqrt(Ctotal^2 + Stotal^2);

    sumRgroup = sum(Rgroup);

    % Estimate the concentration of angles across groups
    rbarWithin = sumRgroup / N;

    if rbarWithin <= 0

        kappa = 0;

    elseif rbarWithin >= 1

        kappa = Inf;

    else

        % Estimate kappa from the mean concentration
        A1 = @(x) ...
            besseli(1, x, 1) ./ besseli(0, x, 1) - rbarWithin;

        upperBound = 1;

        while A1(upperBound) < 0
            upperBound = upperBound * 2;
        end

        kappa = fzero(A1, [0 upperBound]);
    end

    % Apply the Watson-Williams correction based on kappa
    if kappa == 0

        correction = NaN;

    elseif isinf(kappa)

        correction = 1;

    else

        correction = 1 + 3 / (8 * kappa);
    end

    df1 = k - 1;
    df2 = N - k;

    if isnan(correction)

        F = NaN;
        p = NaN;

    else

        numerator = ...
            (N - k) * (sumRgroup - Rtotal);

        denominator = ...
            (k - 1) * (N - sumRgroup);

        F = correction * numerator / denominator;

        p = fcdf(F, df1, df2, 'upper');
    end
end


% Mardia-Watson-Wheeler test
%
% Compares circular distributions across groups using ranked angular data

function [W, df, p] = mardiaWatsonWheeler(angles)

    k = length(angles);

    n = cellfun(@length, angles);
    N = sum(n);

    allAngles = [];
    groupID = [];

    % Combine observations while retaining original group
    for g = 1:k

        allAngles = [allAngles; angles{g}];
        groupID = [groupID; repmat(g, n(g), 1)];
    end

    % Wrap angles to one full circle and assign circular ranks
    allAngles = mod(allAngles, 2*pi);

    ranks = tiedrank(allAngles);

    circularRanks = 2*pi*ranks/N;

    Cgroup = zeros(k, 1);
    Sgroup = zeros(k, 1);

    % Calculate rank-based circular components for each group
    for g = 1:k

        idx = groupID == g;

        Cgroup(g) = sum(cos(circularRanks(idx)));
        Sgroup(g) = sum(sin(circularRanks(idx)));
    end

    % Calculate the W test statistic for two or more groups
    if k == 2

        W = ...
            2 * (N - 1) * ...
            (Cgroup(1)^2 + Sgroup(1)^2) / ...
            (n(1) * n(2));

        df = 2;

    else

        W = 0;

        for g = 1:k

            W = W + ...
                (Cgroup(g)^2 + Sgroup(g)^2) / n(g);
        end

        W = 2 * W;

        df = 2 * (k - 1);
    end

    % Calculate the p-value from the chi-square distribution
    p = chi2cdf(W, df, 'upper');
end


% Holm multiple-comparisons correction
%
% Corrects the pairwise p-values for multiple comparisons

function adjustedP = holmCorrection(pValues)

    pValues = pValues(:);

    adjustedP = NaN(size(pValues));

    valid = ~isnan(pValues);

    p = pValues(valid);

    m = length(p);

    if m == 0
        return
    end

    % Sort p-values before applying the Holm correction
    [sortedP, order] = sort(p);

    adjustedSorted = zeros(m, 1);

    for i = 1:m

        adjustedSorted(i) = ...
            (m - i + 1) * sortedP(i);
    end

    % Keep adjusted p-values in increasing order
    for i = 2:m

        adjustedSorted(i) = ...
            max(adjustedSorted(i), adjustedSorted(i - 1));
    end

    adjustedSorted = min(adjustedSorted, 1);

    % Return corrected p-values to their original comparison order
    corrected = zeros(m, 1);
    corrected(order) = adjustedSorted;

    adjustedP(valid) = corrected;
end