%% Local compilation code 
% To be integrated into the make file
DEAPDir = '/space/ceph/1/ABCD/users/parekh/2026-01-02_DEAP_Parallel/code/DEAP';
FEMADir = '/space/ceph/1/ABCD/users/parekh/2026-01-02_DEAP_Parallel/code/FEMA';
utilDir = '/space/ceph/1/ABCD/users/parekh/2026-01-02_DEAP_Parallel/code/cmig_tools_utils/matlab';
docsDir = '/space/ceph/1/ABCD/users/parekh/2026-01-02_DEAP_Parallel/code/docs';
recipes = '/space/ceph/1/ABCD/users/parekh/2026-01-02_DEAP_Parallel/code/recipes';
surfDir = '/space/ceph/1/ABCD/users/parekh/2026-01-02_DEAP_Parallel/code/showSurf';
volDir  = '/space/ceph/1/ABCD/users/parekh/2026-01-02_DEAP_Parallel/code/showVol';
gifti   = '/space/ceph/1/ABCD/users/parekh/2026-01-02_DEAP_Parallel/code/gifti_toolbox';
reqdDir = [{FEMADir}, {utilDir}, {docsDir}, {recipes}, {surfDir}, {volDir}, {gifti}];
outDir  = '/space/ceph/1/ABCD/users/parekh/2026-01-02_DEAP_Parallel/compiledCode';

%% FEMA_DEAP_wrapper
doCompile(fullfile(DEAPDir, 'FEMA_DEAP_wrapper.m'), reqdDir, outDir);

%% FEMA_DEAP_worker
doCompile(fullfile(DEAPDir, 'FEMA_DEAP_worker.m'), reqdDir, outDir);

%% FEMA_DEAP_gencache
doCompile(fullfile(DEAPDir, 'FEMA_DEAP_gencache.m'), reqdDir, outDir);

function doCompile(inName, reqdDir, outDir)
% Compile app options
opts = compiler.build.StandaloneApplicationOptions(inName, 'TreatInputsAsNumeric', 'Off', ...
                                                   'AdditionalFiles', reqdDir, 'OutputDir', outDir);

% Actual compilation
compiler.build.standaloneApplication(opts);

end