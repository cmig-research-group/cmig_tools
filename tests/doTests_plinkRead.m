%% Read dummy data and perform a battery of read tests
txt = 'Starting tests';
disp([txt, repmat('.', 1, 80 - length(txt))]);

%% Get MATLAB version and suppress warning
matlab_ver      = version('-release');
matlab_ver(end) = '';
matlab_ver      = str2double(matlab_ver);
if matlab_ver   > 2016
    warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
else
    warning('OFF', 'MATLAB:table:ModifiedVarnames');
end

%% Identify location where the test files are located
dirTests = fullfile(fileparts(which('doTests_plinkRead.m')), 'plinkRead');

%% Wide format data - sample size < SNPs
% Settings
bFile_wide   = fullfile(dirTests, 'dummyBED_50samp-100snp-02miss');
nsubjs_wide  = 50;
nsnps_wide   = 100;

% Get ground truth information
truth_wide          = int8(ped2mat([bFile_wide, '_text.ped']));
truth_missing_wide  = readtable([bFile_wide, '_genoCounts.gcount'], 'FileType', 'text');

txt = 'Finished reading wide format data';
disp([txt, repmat('.', 1, 80 - length(txt))]);

% -------------------------------------------------------------------------
% Test 1 - get full data
% -------------------------------------------------------------------------

% Get data using PlinkRead_binary2
fullData_bin2 = PlinkRead_binary2(nsubjs_wide, 1:nsnps_wide, bFile_wide);

% Get data using PlinkRead_binary2_subj
fullData_bin2_subj = PlinkRead_binary2_subj(nsubjs_wide, 1:nsnps_wide, 1:nsubjs_wide, bFile_wide);

% Check if all output are the same
if any(truth_wide - fullData_bin2)
    error('Reading full data: Test failed');
else
    if any(truth_wide - fullData_bin2_subj)
        error('Reading full data: Test failed');
    else
        if any(fullData_bin2 - fullData_bin2_subj)
            error('Reading full data: Test failed');
        else
            txt = 'Reading full data for wide format';
            disp([txt, repmat('.', 1, 80 - length(txt) - 6), 'Passed']);
        end
    end
end

clear fullData*

% -------------------------------------------------------------------------
% Test 2 - get some 13 random SNPs for all subjects
% -------------------------------------------------------------------------

rng(0, 'twister');
whichSNPs = sort(randperm(nsnps_wide, 13));

% Get data using PlinkRead_binary2
SNPData_bin2 = PlinkRead_binary2(nsubjs_wide, whichSNPs, bFile_wide);

% Get data using PlinkRead_binary2_subj
SNPData_bin2_subj = PlinkRead_binary2_subj(nsubjs_wide, whichSNPs, 1:nsubjs_wide, bFile_wide);

% Check if all output are the same
if any(truth_wide(:, whichSNPs) - SNPData_bin2)
    error('Reading subset of SNP data: Test failed');
else
    if any(truth_wide(:, whichSNPs) - SNPData_bin2_subj)
        error('Reading subset of SNP  data: Test failed');
    else
        if any(SNPData_bin2 - SNPData_bin2_subj)
            error('Reading subset of SNP  data: Test failed');
        else
            txt = 'Reading subset of SNPs for wide format';
            disp([txt, repmat('.', 1, 80 - length(txt) - 6), 'Passed']);
        end
    end
end

clear SNPData* which*

% -------------------------------------------------------------------------
% Test 3 - get all SNPs for some 42 random subjects 
% -------------------------------------------------------------------------

rng(10, 'twister');
whichSubjs = randperm(nsubjs_wide, 42);

% Get data using PlinkRead_binary2
allData_subjs_bin2 = PlinkRead_binary2(nsubjs_wide, 1:nsnps_wide, bFile_wide);

% Get data using PlinkRead_binary2_subj
allData_bin2_subj = PlinkRead_binary2_subj(nsubjs_wide, 1:nsnps_wide, whichSubjs, bFile_wide);

% Check if all output are the same
if any(truth_wide(whichSubjs, :) - allData_subjs_bin2(whichSubjs, :))
    error('Reading full data on subset of subjects: Test failed');
else
    if any(truth_wide(whichSubjs, :) - allData_bin2_subj)
        error('Reading full data on subset of subjects: Test failed');
    else
        if any(allData_subjs_bin2(whichSubjs, :) - allData_bin2_subj)
            error('Reading full data on subset of subjects: Test failed');
        else
            txt = 'Reading subset of subjects for wide format';
            disp([txt, repmat('.', 1, 80 - length(txt) - 6), 'Passed']);
        end
    end
end

clear allData* which*

% -------------------------------------------------------------------------
% Test 4 - get 33 random SNPs for 27 random subjects
% -------------------------------------------------------------------------

rng(100, 'twister');
whichSubjs = randperm(nsubjs_wide, 27);
whichSNPs  = sort(randperm(nsnps_wide,  33));

% Get data using PlinkRead_binary2
SNPData_subjs_bin2 = PlinkRead_binary2(nsubjs_wide, whichSNPs, bFile_wide);

% Get data using PlinkRead_binary2_subj
SNPData_bin2_subj = PlinkRead_binary2_subj(nsubjs_wide, whichSNPs, whichSubjs, bFile_wide);

% Check if all output are the same
if any(truth_wide(whichSubjs, whichSNPs) - SNPData_subjs_bin2(whichSubjs, :))
    error('Reading subset of SNP for a subset of subjects: Test failed');
else
    if any(truth_wide(whichSubjs, whichSNPs) - SNPData_bin2_subj)
        error('Reading subset of SNP for a subset of subjects: Test failed');
    else
        if any(SNPData_subjs_bin2(whichSubjs, :) - SNPData_bin2_subj)
            error('Reading subset of SNP for a subset of subjects: Test failed');
        else
            txt = 'Reading subset of SNPs and subjects for wide format';
            disp([txt, repmat('.', 1, 80 - length(txt) - 6), 'Passed']);
        end
    end
end

clear SNPData* which*

% -------------------------------------------------------------------------
% Test 5 - missingness rate
% -------------------------------------------------------------------------

% Get data using PlinkRead_binary2
fullData_wide_bin2 = PlinkRead_binary2(nsubjs_wide, 1:nsnps_wide, bFile_wide);

% Get data using PlinkRead_binary2_subj
fullData_wide_bin2_subj = PlinkRead_binary2_subj(nsubjs_wide, 1:nsnps_wide, 1:nsubjs_wide, bFile_wide);

% SNP-wise missing information
% Sum should be across rows
missing_wide_bin2       = sum(fullData_wide_bin2      == -1, 1);
missing_wide_bin2_subj  = sum(fullData_wide_bin2_subj == -1, 1);

% Compare missing counts
% Transpose is necessary in this case
if any(truth_missing_wide.MISSING_CT - missing_wide_bin2') ||   ...
   any(truth_missing_wide.MISSING_CT - missing_wide_bin2_subj')
    error('Missing count test failed');
else
    txt = 'Comparing missingness counts for wide format';
    disp([txt, repmat('.', 1, 80 - length(txt)  - 6), 'Passed']);
end

% Finished tests on wide format data
txt = 'Reading running tests on wide format';
disp([txt, repmat('.', 1, 80 - length(txt))]);

%% Long format data: sample size > SNPs
% Settings
bFile_long   = fullfile(dirTests, 'dummyBED_100samp-50snp-02miss');
nsubjs_long  = 100;
nsnps_long   = 50;

% Get ground truth information
truth_long          = int8(ped2mat([bFile_long, '_text.ped']));
truth_missing_long  = readtable([bFile_long, '_genoCounts.gcount'], 'FileType', 'text');

txt = 'Finished reading long format data';
disp([txt, repmat('.', 1, 80 - length(txt))]);

% -------------------------------------------------------------------------
% Test 1 - get full data
% -------------------------------------------------------------------------

% Get data using PlinkRead_binary2
fullData_bin2 = PlinkRead_binary2(nsubjs_long, 1:nsnps_long, bFile_long);

% Get data using PlinkRead_binary2_subj
fullData_bin2_subj = PlinkRead_binary2_subj(nsubjs_long, 1:nsnps_long, 1:nsubjs_long, bFile_long);

% Check if all output are the same
if any(truth_long - fullData_bin2)
    error('Reading full data: Test failed');
else
    if any(truth_long - fullData_bin2_subj)
        error('Reading full data: Test failed');
    else
        if any(fullData_bin2 - fullData_bin2_subj)
            error('Reading full data: Test failed');
        else
            txt = 'Reading full data for long format';
            disp([txt, repmat('.', 1, 80 - length(txt) - 6), 'Passed']);
        end
    end
end

clear fullData*

% -------------------------------------------------------------------------
% Test 2 - get some 13 random SNPs for all subjects
% -------------------------------------------------------------------------

rng(0, 'twister');
whichSNPs = sort(randperm(nsnps_long, 13));

% Get data using PlinkRead_binary2
SNPData_bin2 = PlinkRead_binary2(nsubjs_long, whichSNPs, bFile_long);

% Get data using PlinkRead_binary2_subj
SNPData_bin2_subj = PlinkRead_binary2_subj(nsubjs_long, whichSNPs, 1:nsubjs_long, bFile_long);

% Check if all output are the same
if any(truth_long(:, whichSNPs) - SNPData_bin2)
    error('Reading subset of SNP data: Test failed');
else
    if any(truth_long(:, whichSNPs) - SNPData_bin2_subj)
        error('Reading subset of SNP  data: Test failed');
    else
        if any(SNPData_bin2 - SNPData_bin2_subj)
            error('Reading subset of SNP  data: Test failed');
        else
            txt = 'Reading subset of SNPs for wide format';
            disp([txt, repmat('.', 1, 80 - length(txt) - 6), 'Passed']);
        end
    end
end

clear SNPData* which*

% -------------------------------------------------------------------------
% Test 3 - get all SNPs for some 42 random subjects 
% -------------------------------------------------------------------------

rng(10, 'twister');
whichSubjs = randperm(nsubjs_long, 42);

% Get data using PlinkRead_binary2
allData_subjs_bin2 = PlinkRead_binary2(nsubjs_long, 1:nsnps_long, bFile_long);

% Get data using PlinkRead_binary2_subj
allData_bin2_subj = PlinkRead_binary2_subj(nsubjs_long, 1:nsnps_long, whichSubjs, bFile_long);

% Check if all output are the same
if any(truth_long(whichSubjs, :) - allData_subjs_bin2(whichSubjs, :))
    error('Reading full data on subset of subjects: Test failed');
else
    if any(truth_long(whichSubjs, :) - allData_bin2_subj)
        error('Reading full data on subset of subjects: Test failed');
    else
        if any(allData_subjs_bin2(whichSubjs, :) - allData_bin2_subj)
            error('Reading full data on subset of subjects: Test failed');
        else
            txt = 'Reading subset of subjects for long format';
            disp([txt, repmat('.', 1, 80 - length(txt) - 6), 'Passed']);
        end
    end
end

clear allData* which*

% -------------------------------------------------------------------------
% Test 4 - get 33 random SNPs for 27 random subjects
% -------------------------------------------------------------------------

rng(100, 'twister');
whichSubjs = randperm(nsubjs_long, 27);
whichSNPs  = sort(randperm(nsnps_long,  33));

% Get data using PlinkRead_binary2
SNPData_subjs_bin2 = PlinkRead_binary2(nsubjs_long, whichSNPs, bFile_long);

% Get data using PlinkRead_binary2_subj
SNPData_bin2_subj = PlinkRead_binary2_subj(nsubjs_long, whichSNPs, whichSubjs, bFile_long);

% Check if all output are the same
if any(truth_long(whichSubjs, whichSNPs) - SNPData_subjs_bin2(whichSubjs, :))
    error('Reading subset of SNP for a subset of subjects: Test failed');
else
    if any(truth_long(whichSubjs, whichSNPs) - SNPData_bin2_subj)
        error('Reading subset of SNP for a subset of subjects: Test failed');
    else
        if any(SNPData_subjs_bin2(whichSubjs, :) - SNPData_bin2_subj)
            error('Reading subset of SNP for a subset of subjects: Test failed');
        else
            txt = 'Reading subset of SNPs and subjects for long format';
            disp([txt, repmat('.', 1, 80 - length(txt) - 6), 'Passed']);
        end
    end
end

clear SNPData* which*

% -------------------------------------------------------------------------
% Test 5 - missingness rate
% -------------------------------------------------------------------------

% Get data using PlinkRead_binary2
fullData_long_bin2 = PlinkRead_binary2(nsubjs_long, 1:nsnps_long, bFile_long);

% Get data using PlinkRead_binary2_subj
fullData_long_bin2_subj = PlinkRead_binary2_subj(nsubjs_long, 1:nsnps_long, 1:nsubjs_long, bFile_long);

% SNP-wise missing information
% Sum should be across rows
missing_long_bin2       = sum(fullData_long_bin2      == -1, 1);
missing_long_bin2_subj  = sum(fullData_long_bin2_subj == -1, 1);

% Compare missing counts
% Transpose is necessary
if any(truth_missing_long.MISSING_CT - missing_long_bin2') ||   ...
   any(truth_missing_long.MISSING_CT - missing_long_bin2_subj')
    error('Missing count test failed');
else
    txt = 'Comparing missingness counts for long format';
    disp([txt, repmat('.', 1, 80 - length(txt) - 6), 'Passed']);
end

% Finished tests on wide format data
txt = 'Finished running tests on long format';
disp([txt, repmat('.', 1, 80 - length(txt))]);

%% Turn warnings back on
if matlab_ver > 2016
    warning('ON', 'MATLAB:table:ModifiedAndSavedVarnames');
else
    warning('ON', 'MATLAB:table:ModifiedVarnames');
end