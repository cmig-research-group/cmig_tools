% Script for running ABCD_DEMO tests before merging branches
% Diana Smith, Nov 2022
%
% This script is designed to write out a "diary" file with the outputs of FEMA for the five modalities run in the 
% FEMA_wrapper_demo.m script. if folder==today, then the script will run FEMA_wrapper_demo.m. If folder==old, then 
% this script will print the outputs from a previous date the script was run.
%
% As of now, in order for the script to work properly you need to have all of your FEMA test results stored in
% a directory ABCD_DEMOS_results, with each set of results stored in a subdirectory corresponding to that day's date.
%
% Because this script calls "FEMA_wrapper_demo.m" it will still require the user to specify inputs from within 
% the demo script. This is something to change in the future.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up diary

% add path to the cmig_tools_internal directory
git_dir=input('Please specify the path to the cmig_tools directory you wish to test (make sure you have switched to the correct branch). \nIf this is already in your path leave blank and press ENTER: ','s');
addpath(genpath(git_dir))

% set "old" to whichever date you want to use as your comparison
today = datestr(now,'dd-mmm-yyyy','local');
old = '2022-11-09';

folder = today;

if folder == today
    run(fullfile(git_dir, 'ABCD_DEMOS', 'FEMA_wrapper_demo.m'))
end

diary(['diary_',datestr(now,'dd-mm-yy','local'),'_',datestr(now,'hh-MM-ss','local'),'.txt'])
logging(sprintf('folder: %s', folder));

fnames_results = {'FEMA_wrapper_output_corrmat_rsfmri_fd0.20.mat';...
'FEMA_wrapper_output_external_RNI.mat'; 'FEMA_wrapper_output_vertex_thickness_ic5_sm256_synth.mat'; ...
'FEMA_wrapper_output_vertex_thickness_ic5_sm256.mat'; 'FEMA_wrapper_output_voxel_RNI.mat'};

% TODO: generate list of indexes rather than just 1:10

for i = 1:length(fnames_results)
    disp(fnames_results{i});
    load(fullfile('ABCD_DEMOS_results',folder,'FEMA_demo_rel4.0_designmatrix_age',fnames_results{i}))

    disp('whos');
    whos
    disp('X');
    X(1:10)
    disp('beta_hat');
    beta_hat(1:10)
    disp('beta_se');
    beta_se(1:10)
    disp('eid');
    eid(1:10)
    disp('iid');
    iid(1:10)
    disp('logpmat');
    logpmat(1:10)
    disp('sig2mat');
    sig2mat(:,1:10)
    disp('zmat');
    zmat(1:10)

    clear X beta_hat beta_se colnames_imaging colnames_model contrasts datatype eid iid logpmat mask sig2mat sig2tvec zmat;
end

logging(sprintf('Diary saved to %s',get(0,'DiaryFile')));

diary off