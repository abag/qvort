function res = stackSlice(img, dir, slice)
%STACKSLICE Extract a planar slice from a 3D image
%
%   SLICE = stackSlice(IMG, DIR, INDEX)
%   IMG is either a 3D or a 4D (3D+color) image.
%   DIR is 1, 2 or 3, or 'x', 'y', 'z'
%   INDEX is the slice index, between 1 and the number of voxels in the DIR
%   direction.
%
%   The functions returns a planar image the same type as the original.
%
%
%   Example
%   % Display 3 slices of a MRI head
%     img = analyze75read(analyze75info('brainMRI.hdr'));
%     figure(1); clf; hold on;
%     for i=1:3
%         subplot(3, 1, i);
%         imshow(stackSlice(img, 'x', 20+20*i)');
%         set(gca, 'ydir', 'normal')
%     end
%
%   See also
%
%
%   ------
%   author: David Legland, david.legland(at)grignon.inra.fr
%   INRA - Cepia Software Platform
%   Created: 2007-08-14,    using Matlab 7.4.0.287 (R2007a)
%   http://www.pfl-cepia.inra.fr/index.php?page=slicer
%   Licensed under the terms of the new BSD license, see file license.txt


% image size
dim = size(img);

% convert to an index between 1 and 3
dir = parseAxisIndex_ijk(dir);


if length(dim)==3
    % gray-scale image
    switch dir
        case 1, res = squeeze(img(slice, :, :));
        case 2, res = squeeze(img(:, slice, :));
        case 3, res = img(:, :, slice);
        otherwise, error('Wrong direction');
    end

elseif length(dim)==4
    % case of color image
    switch dir
        case 1, res = squeeze(permute(img(slice, :, :,:), [2 4 3 1]));
        case 2, res = squeeze(permute(img(:, slice, :, :), [1 4 3 2]));
        case 3, res = img(:, :, :, slice);
        otherwise, error('Wrong direction');
    end
else
    error('Input array should have dimension 3 or 4');
end
