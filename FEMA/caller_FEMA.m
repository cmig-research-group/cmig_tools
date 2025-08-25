function caller_FEMA(moduleName, varargin)
% Caller function for all FEMA modules (compiled version)
%% Inputs:
% --------
% moduleName:       should be one of the following (case in-sensitive):
%                   values separated by OR (|) are equivalent
%                       * help
%                       * FEMA_fit | fit
%                       * createBasisFunctions | createBF | makeBF
%                       * FEMA_fit_GWAS | GWAS | fitGWAS
% 
%% Modules to be added/worked on
% * makeContrasts | parseContrasts
% * evaluateContrasts

%% Make decisions based on moduleName
if ~exist('moduleName', 'var') || isempty(moduleName)
    error('No FEMA module specified');
else
    validModules = {'help', 'h', ...
                    'FEMA_fit', 'fit', ...
                    'createBasisFunctions', 'createBF', 'makeBF', ...
                    'FEMA_fit_GWAS', 'GWAS', 'fitGWAS', ...
                    'makeContrasts', 'parseContrasts', ...
                    'evaluateContrasts'};

    if ~ismember(lower(moduleName), lower(validModules))
        error(['Unknown FEMA module specified: ', moduleName]);
    else
        p                 = inputParser;
        p.KeepUnmatched   = true;
        validationFcn_str = @(s) isstring(s) | ischar(s);

        switch lower(moduleName)
            case {'help', 'h'}
                if isempty(varargin)
                    disp('Please use the syntax: help <name of the module you need help for>');
                else
                    dirHelp = fullfile(fileparts(fileparts(which('caller_FEMA'))), 'docs', 'help');

                    if ~exist(dirHelp, 'dir')
                        disp(dirHelp);
                        warning('Could not find the help folder');
                    else
                        forHelp = varargin{1};

                        switch lower(forHelp)
                            case {'fema_fit', 'fit'}
                                system(['cat ', fullfile(dirHelp, 'caller_FEMA_fit.txt')]);

                            case {'createbasisfunctions', 'createbf', 'makebf'}
                                system(['cat ', fullfile(dirHelp, 'caller_createBasisFunctions.txt')]);

                            case {'fema_fit_gwas', 'gwas', 'fitgwas'}
                                system(['cat ', fullfile(dirHelp, 'caller_FEMA_fit_GWAS.txt')]);

                            otherwise
                                disp('Invoking help of non-compiled versions');
                                toDisp = fullfile(dirHelp, [strrep(forHelp, '.m', ''), '.txt']);
                                if exist(toDisp, 'file')
                                    system(['cat ', toDisp]);
                                else
                                    warning(['No help found for: ', forHelp]);
                                end
                        end
                    end
                end

            case {'fema_fit', 'fit'}
                def_outPrefix = ['FEMA_fit-', char(datetime('now', 'Format', 'yyyyMMMdd-HHmmSS'))];

                % Add positional and optional required arguments
                addRequired(p, 'file_X',    validationFcn_str);
                addRequired(p, 'file_ymat', validationFcn_str);
                addRequired(p, 'dirOutput', validationFcn_str);
                addOptional(p, 'outPrefix', def_outPrefix, validationFcn_str);

                % Parse inputs
                parse(p, varargin{:});

                % Call FEMA
                caller_FEMA_fit(p.Results.file_X,    p.Results.file_ymat, ...
                                p.Results.dirOutput, p.Results.outPrefix, p.Unmatched);

            case {'createbasisfunctions', 'createbf', 'makebf'}
                def_outPrefix = ['FEMA_createBasisFunctions-', ...
                                 char(datetime('now', 'Format', 'yyyyMMMdd-HHmmSS'))];

                % Add positional and optional required arguments
                addRequired(p, 'file_valvec', validationFcn_str);
                addRequired(p, 'dirOutput',   validationFcn_str);
                addRequired(p, 'outType',     validationFcn_str);
                addOptional(p, 'outPrefix', def_outPrefix, validationFcn_str);

                % Parse inputs
                parse(p, varargin{:});

                % Call createBasisFunctions
                caller_createBasisFunctions(p.Results.file_valvec, p.Results.dirOutput, ...
                                            p.Results.outType, p.Results.outPrefix, p.Unmatched);

            case {'fema_fit_gwas', 'gwas', 'fitgwas'}
                def_outPrefix = '';

                % Add positional and optional required arguments
                addRequired(p, 'file_PLINK',    validationFcn_str);
                addRequired(p, 'file_FEMA_fit', validationFcn_str);
                addRequired(p, 'GWASType',      validationFcn_str);
                addRequired(p, 'dirOutput',     validationFcn_str);
                addOptional(p, 'outPrefix', def_outPrefix, validationFcn_str);

                % Parse inputs
                parse(p, varargin{:});

                % Call FEMA_fit_GWAS
                caller_FEMA_fit_GWAS(p.Results.file_PLINK, p.Results.file_FEMA_fit, ...
                                     p.Results.GWASType, p.Results.dirOutput, p.Results.outPrefix);
        end
    end
end