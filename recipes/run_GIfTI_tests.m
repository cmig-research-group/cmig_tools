% This function showcases the use of writeGIfTI_underlay and writeGIfTI to
% generate surface underlay and ovelay files

%% Test underlays
% Since no names are provided, output goes to pwd
dirFreesurfer = '/Applications/freesurfer/7.2.0';
icoNum        = 5:7;
underlays     = {'pial', 'inflated', 'white', 'aparc', 'aparc.a2009s', 'yeo17', 'yeo7'};

% Test without splitting LR
splitLR = false;
for ic = 1:length(icoNum)
    for un = 1:length(underlays)
        writeGIfTI_underlay(icoNum(ic), underlays{un}, dirFreesurfer, [], splitLR);
    end
end

% Test with splitting LR
splitLR = true;
for ic = 1:length(icoNum)
    for un = 1:length(underlays)
        writeGIfTI_underlay(icoNum(ic), underlays{un}, dirFreesurfer, [], splitLR);
    end
end

%% Test overlays
% Generate random data for five statistical coefficients (ico = 5)
icoNum = 5;
data   = randn(5, 10242*2);

% Do not split left and right
writeGIfTI(data, icoNum, [], false);

% Split left and right
writeGIfTI(data, icoNum, [], true);