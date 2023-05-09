function resMatrix = ped2mat(pedFile)
% Function that reads a ped file and returns a -1, 0, 1, and 2 coded matrix
%% Input(s):
% pedFile:      full path to a ped file, including extension
%
%% Output(s):
% resMatrix:    a subj x SNP matrix where -1 indicates missing data, 0
%               indicates homozygous recessive, 1 indicates heterozygous,
%               and 2 indicates homozygous dominant

%% Check inputs
if ~exist('pedFile', 'var') || isempty(pedFile)
    error('Please provide full path to a valid PED file');
else
    if ~exist(pedFile, 'file')
        error(['Unable to find: ', pedFile]);
    end
end

%% Read text file
data = readtable(pedFile, 'FileType', 'text');

% First six columns can be ignored
data(:, 1:6) = [];

% Convert to matrix
data = data{:,:};

% Convert to characters
gtruth_txt = cellfun(@(x) strrep(x, ' ', ''), cellstr(num2str(data)), 'UniformOutput', false);

% Every two characters form one SNP value
beginPositions = 1:2:size(data,2);

%% Resulting matrix
resMatrix        = NaN(size(data, 1), length(beginPositions));
tmpR             = cell(size(data, 1), length(beginPositions));
for locs         = 1:length(beginPositions)
    tmpR(:,locs) = cellfun(@(x) x(beginPositions(locs):beginPositions(locs)+1), gtruth_txt, 'UniformOutput', false);
end

%% Replace with results
% 00        is -1
% 11        is  0
% 01 and 10 are 1
% 02 and 20 are 1
% 12 and 21 are 1
% 22        is  2
resMatrix(ismember(tmpR, '00'))                         = -1;
resMatrix(ismember(tmpR, '11'))                         =  0;
resMatrix(ismember(tmpR, '01') | ismember(tmpR, '10'))  =  1;
resMatrix(ismember(tmpR, '02') | ismember(tmpR, '20'))  =  1;
resMatrix(ismember(tmpR, '12') | ismember(tmpR, '21'))  =  1;
resMatrix(ismember(tmpR, '22'))                         =  2;