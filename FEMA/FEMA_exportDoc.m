function FEMA_exportDoc(inDir, outDir)
% Function that exports the header of all .m files in a specified inDir to
% a plain text file in outDir, removing all occurrances of % and %%
% 
% Assumes that documentation of an input file ends when we encounter
% discontinuity in lines begining with %
% 
%% Inputs:
% inDir:        full path to an input directory having .m files
% 
% outDir:       full path to where the .txt files should be saved

%% Make a list of .m files
listFiles = dir(fullfile(inDir, '*.m'));

%% Loop over files, parse them, save them
for files = 1:length(listFiles)

    % Read the file
    % https://www.mathworks.com/matlabcentral/answers/455574#answer_370050
    fid     = fopen(fullfile(inDir, listFiles(files).name), 'rt');
    content = fread(fid, '*char')';
    fclose(fid);

    % Split the content by lines
    content = regexp(content,'\n','split')';
    % content = strsplit(content, '\n')';

    % Find all lines starting with '%'
    temp     = ~cellfun(@isempty, regexpi(content, '^%'));
    allLocs  = find(temp);
    % content  = content(temp);

    if isempty(allLocs)
        % No documentation found in this file; write out an empty file
        content = {''};
        fid = fopen(fullfile(outDir, strrep(listFiles(files).name, '.m', '.txt')), 'w');
        fprintf(fid, '%s\n', content{:});
        fclose(fid);
        continue;
    else
        locStart = allLocs(1);
        locEnd   = find(~temp); % allLocs(find(diff(allLocs) > 1, 1));
        locEnd   = locEnd(find(locEnd > locStart, 1)) - 1;

        if isempty(locEnd)
            locEnd = allLocs(end);
        end

        % Subset content
        content  = content(locStart:locEnd);

        % Replace all occurrences of %% and %
        content = strrep(strrep(content, '%%', ''), '%', '');
    
        % Write out
        fid = fopen(fullfile(outDir, strrep(listFiles(files).name, '.m', '.txt')), 'w');
        fprintf(fid, '%s\n', content{:});
        fclose(fid);
    end

    % Display status
    disp(['Completed: ', listFiles(files).name]);
end