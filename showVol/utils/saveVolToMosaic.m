function saveVolToMosaic(vol, plane, name)
% saveVolMosaicPng Save volume as a mosaic of axial slices, in png, for DEAP
%
%  saveVolMosaic(fname,plane, vol)
%
%   vol   imaging volume. vol.imgs is ABCD brain volume
%   plane slice plane to make mosaic: Axial, Sagittal, Coronal
%   name  path/name describing volume, saved as '<name>_<plane>.png' e.g. name='/path/T1_ABCD2_cor10'
%
% JRI (jiversen@ucsd.edu)
%
% original code from Hauke

doScale16 = true;
doPlot = true;

data = vol.imgs;

planes = size(data,4); %bit planes in image
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

%data2d8bit = uint8(double(data2d)/double(max(max(data2d))) * (2^8-1));

if ~doScale16
  data2d16bit = uint16(data2d);
else
  data2d16bit = uint16(double(data2d)/double(max(data2d(:))) * (2^16-1));
end

opts={'CreationTime',datestr(now)};
if isfield(vol,'name')
  opts = [opts 'Description', vol.name];
end

fname = strcat(name,'_', plane, '.png');
disp(['Saving to ' fname])

imwrite(data2d16bit,fname,'PNG','BitDepth',16,opts{:}) 

if doPlot
  figure
  imshow(data2d16bit)
  title(fname,'interpreter','none')
end

  