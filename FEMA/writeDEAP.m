function writeDEAP(results, dirname_out, fstem_imaging, ivnames, colnames_model)
% writeDEAP write voxelwise results to binary format mosaics for DEAP, called from SSE_wrapper
%
%   writeDEAP(results, dirname_out, fstem_imaging, ivnames, colnames_model)
%
% writes a volume for each variable type in results and each IV in ivnames and in each of three planes

randomFields = {'sig2tvec', 'sig2mat'};
planes = {'Axial','Sagittal','Coronal'};

% parse IVs
if isempty(ivnames)
  excludeCol = strmatch('mri_info_',colnames_model);
  nCol = length(colnames_model);
  ivCol = setdiff(1:nCol, excludeCol);
else
  [~,ivCol,~] = intersect(colnames_model,ivnames);
end
if length(ivCol) < 1, error('No IVs found! Not writing DEAP.'), end

% =========================================================================
% write out main effects (for IVs of interest)
fieldnamelist = setdiff(fieldnames(results),randomFields);
for fi = 1:length(fieldnamelist)
  fieldname = fieldnamelist{fi};
  vol = results.(fieldname);  %all IVs
  name = [fstem_imaging '_' fieldname];
  
  %write one file per IV
  for iv = ivCol(:)'
    %colname = colnames_model{iv};
    data = vol(:,:,:,iv); %2mm voxel
    data = upsample_volume(data); %1mm voxel
    d = size(data);
    for iP = 1:length(planes)
      fname = sprintf('%s/%s_%s_%02d_%ssingle.dat',dirname_out, name, planes{iP}, iv, num2str(d,"%d_"));
      writeDeapMosaic(data,planes{iP},fname)
    end
  end %loop over IVs
end %loop on data type

% =========================================================================
% write out the random effects
fieldnamelist = randomFields;
for fi = 1:length(fieldnamelist)
  fieldname = fieldnamelist{fi};
  vol = results.(fieldname);  %all IVs
  name = [fstem_imaging '_' fieldname];
  
  %write one file per volume
  for iv = 1:size(vol,4)
    data = vol(:,:,:,iv); %2mm voxel
    data = upsample_volume(data); %1mm voxel
    d = size(data);
    for iP = 1:length(planes)
      fname = sprintf('%s/%s_%s_%02d_%ssingle.dat',dirname_out, name, planes{iP}, iv, num2str(d,"%d_"));
      writeDeapMosaic(data,planes{iP},fname)
    end
  end %loop over volumes
end %loop over random effects

% =========================================================================
%write column names to json for DEAP
%colnames_model = colnames_model(ivCol); %for the iv index to be used to lookup, need to write all IV names
fname_col = sprintf('%s/SSE_results_colnames.json',dirname_out);
out = struct('colnames_model',{colnames_model});
jsonStr = jsonencode(out);
fid = fopen(fname_col,'w');
fprintf(fid,'%s\n',jsonStr);
fclose(fid);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper to write out mosaic in a given plane

function writeDeapMosaic(data, plane, fname)
% saveVolToMosaic Save volume as a mosaic of axial slices, in single precision binary, for DEAP
%
%  writeDeapMosaic(data, plane, name)
%
%   data   imaging volume (3d)
%   plane  slice plane to make mosaic: Axial, Sagittal, Coronal
%   fname  filename for save
%
% JRI (jiversen@ucsd.edu)
%
% original code from Hauke


planes = size(data,4); %bit planes in image - will be one for this usage
d = size(data);

switch lower(plane)
  case 'axial'
    numSlices = d(1);
    imgWidth = d(3);
    imgHeight = d(2);
    numImages = ceil(sqrt(numSlices));
    d2d = [numImages * imgWidth, numImages * imgHeight]; %NB, because of rot90 w,h are swapped for axial
    data2d = single(zeros([d2d planes]));
    h = 1; w = 1;
    for i=1:numSlices
      im = squeeze(data(i,:,:,:));
      data2d(w:(w+imgWidth-1), h:(h+imgHeight-1),:) = rot90(im,1);
      if (w+imgWidth-1) >= d2d(1)
        w = 1;
        h = h + imgHeight;
      else
        w = w + imgWidth;
      end
    end
    
  case 'sagittal'
    numSlices = d(2);
    imgWidth = d(3);
    imgHeight = d(1);
    numImages = ceil(sqrt(numSlices));
    d2d = [numImages * imgHeight, numImages * imgWidth];
    data2d = single(zeros([d2d planes]));
    h = 1; w = 1;
    for i=1:numSlices
      im = squeeze(data(:,i,:,:));
      data2d(h:(h+imgHeight-1),w:(w+imgWidth-1),:) = im;
      if (h+imgHeight-1) >= d2d(1)
        h = 1;
        w = w + imgWidth;
      else
        h = h + imgHeight;
      end
    end
    
  case 'coronal'
    numSlices = d(3);
    imgWidth = d(2);
    imgHeight = d(1);
    numImages = ceil(sqrt(numSlices));
    d2d = [numImages * imgHeight, numImages * imgWidth];
    data2d = single(zeros([d2d planes]));
    h = 1; w = 1;
    for i=1:numSlices
      im = squeeze(data(:,:,i,:));
      data2d(h:(h+imgHeight-1),w:(w+imgWidth-1),:) = im;
      if (h+imgHeight-1) >= d2d(1)
        h = 1;
        w = w + imgWidth;
      else
        h = h + imgHeight;
      end
    end
    
  otherwise
    error('invalid plane')
    
end

%% Save image
disp(['Saving to ' fname])

fid = fopen(fname,'w');
fwrite(fid, data2d', 'single');
fclose(fid);