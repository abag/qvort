function varargout = showYSlice(img, sliceIndex)
%SHOWYSLICE Show XZ slice of a 3D image
%
%   showYSlice(IMG, INDEX)
%   Display the given slice as a 3D planar image. INDEX is the slice index,
%   between 1 and stackSize(img, 1).
%
%   Example
%   % Display orthoslices of a humain head
%   img = analyze75read(analyze75info('brainMRI.hdr'));
%   figure(1); clf; hold on;
%   showZSlice(img, 13);
%   showXSlice(img, 60);
%   showYSlice(img, 80);
%   axis equal
%   xlabel('x'); ylabel('y'); zlabel('z');
%
%   See also
%   showXSlice, showZSlice, getSlice
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-06-30,    using Matlab 7.9.0.529 (R2009b)
% http://www.pfl-cepia.inra.fr/index.php?page=slicer
% Copyright 2010 INRA - Cepia Software Platform.


%% Extract image info

dim = stackSize(img);

% compute voxel positions
ly = 1:dim(2);

% position vectors of voxel corners
vx = ((0:dim(1))+.5);
vz = ((0:dim(3))+.5);

% global parameters for surface display
params = {'facecolor','texturemap', 'edgecolor', 'none'};

% compute position of voxel vertices in 3D space
[xz_x xz_z] = meshgrid(vx, vz);
xz_y = ones(size(xz_x))*ly(sliceIndex);

% extract slice in y direction
slice = stackSlice(img, 'y', sliceIndex);

% eventually converts to uint8, rescaling data between 0 and max value
if ~strcmp(class(slice), 'uint8')
    slice = double(slice);
    slice = uint8(slice*255/max(slice(:)));
end

% convert grayscale to rgb (needed by 'surface' function)
if length(size(slice))==2
    slice = repmat(slice, [1 1 3]);
end

% for ZX slice, we permute
slice = permute(slice, [2 1 3]);

% repeat slice three times to manage a color image
hs = surface(xz_x, xz_y, xz_z, slice, params{:});


%% process output arguments

if nargout>0
    varargout{1} = hs;
end
