%% set up paths
atlas_id = 'ABCD3_cor10';
dirname_atlas = '';
fname_atlas = fullfile(dirname_atlas, sprintf('atlas_dspace_%s.mat', atlas_id));
dirname_out = strrep(fileparts(which('FEMA_fit')), 'FEMA', 'support_files');

%% volumetric parcellations from ABCD3
fname_out = fullfile(dirname_out, 'roi2ABCDAtlasMaps.mat');
roi2ABCD3 = roi2ABCDAtlas(fname_atlas, fname_out);

% save as luts
aseg = table(roi2ABCD3.aseg.roinames, roi2ABCD3.aseg.roicodes, roi2ABCD3.aseg.roirgb, 'VariableNames', {'roinames', 'roicodes', 'roirgb'});
aparcaseg = table(roi2ABCD3.aparcaseg.roinames, roi2ABCD3.aparcaseg.roicodes, roi2ABCD3.aparcaseg.roirgb, 'VariableNames', {'roinames', 'roicodes', 'roirgb'});
fiber = table(roi2ABCD3.fiber.roinames, roi2ABCD3.fiber.roicodes, roi2ABCD3.fiber.roirgb, 'VariableNames', {'roinames', 'roicodes', 'roirgb'});

fname_aseg = fullfile(dirname_out, 'aseg_lut.csv');
fname_aparcaseg = fullfile(dirname_out, 'aparcaseg_lut.csv');
fname_fiber = fullfile(dirname_out, 'fiber_lut.csv');
writetable(aseg, fname_aseg);
writetable(aparcaseg, fname_aparcaseg);
writetable(fiber, fname_fiber);

%% surface parcellations 
icoNum = 5;
splitLR = true;
fname_out = fullfile(dirname_out, 'roi2SurfaceAtlasMaps.mat');
roi2Surf = roi2FreesurferAtlas(icoNum, fname_out, 'splitLR', splitLR);

% save as luts
% aparc
aparc_lh = table(roi2Surf.aparc.roinames_lh, roi2Surf.aparc.roicodes_lh, ...
              'VariableNames', {'roinames', 'roicodes'});
aparc_rh = table(roi2Surf.aparc.roinames_rh, roi2Surf.aparc.roicodes_rh, ...
              'VariableNames', {'roinames', 'roicodes'});
aparc = [aparc_lh; aparc_rh];

% aparc.a2009s
aparc_a2009s_lh = table(roi2Surf.aparc_a2009s.roinames_lh, roi2Surf.aparc_a2009s.roicodes_lh, ...
                     'VariableNames', {'roinames', 'roicodes'});    
aparc_a2009s_rh = table(roi2Surf.aparc_a2009s.roinames_rh, roi2Surf.aparc_a2009s.roicodes_rh, ...
                     'VariableNames', {'roinames', 'roicodes'});    
aparc_a2009s = [aparc_a2009s_lh; aparc_a2009s_rh];

fname_aparc = fullfile(dirname_out, 'aparc_lut.csv');
fname_aparc_a2009s = fullfile(dirname_out, 'aparc_a2009s_lut.csv');
writetable(aparc, fname_aparc);
writetable(aparc_a2009s, fname_aparc_a2009s);