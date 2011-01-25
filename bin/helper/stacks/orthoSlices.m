function orthoSlices(img, varargin)
%ORTHOSLICES  Show three orthogonal slices of a 3D image
%
%   orthoSlices(IMG, POS)
%   POS is 1*3 row vector containing position of slices intersection point,
%   in image index coordinate between 1 and image size, in order [XPOS,
%   YPOS ZPOS].
%
%   Example
%     img = analyze75read(analyze75info('brainMRI.hdr'));
%     figure(1); clf; hold on;
%     orthoSlices(img, [60 80 13]);
%     axis equal;                          % to have equal sizes
%
%   See also
%   stackSlice, showXSlice, showYSlice, showZSlice
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-06-30,    using Matlab 7.9.0.529 (R2009b)
% http://www.pfl-cepia.inra.fr/index.php?page=slicer
% Copyright 2010 INRA - Cepia Software Platform.

% get stack size (in x, y, z order)
siz = stackSize(img);

% use a default position if not specified
if isempty(varargin)
    pos = ceil(siz/2);
else
    pos = varargin{1};
end

% display three orthogonal slices
hold on;
showXSlice(img, pos(1));
showYSlice(img, pos(2));
showZSlice(img, pos(3));

% compute display extent (add a 0.5 limit around each voxel)
extent = [zeros(1, 3) ; siz] + .5;
extent = extent(:)';

% setup display
axis equal;
axis(extent);
