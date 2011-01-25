% STACKS Manipulation and display of 3D images
% Version 0.7 02-Jul-2010 .
%
% Visualisation
%   slicer         - Interactive visualization of 3D images
%   orthoSlices    - Show three orthogonal slices of a 3D image
%   showXSlice     - Show YZ slice of a 3D image
%   showYSlice     - Show XZ slice of a 3D image
%   showZSlice     - Show XY slice of a 3D image
%
% Read/Write 3D images
%   readstack      - Read either a list of 2D images (slices), or a 3D image
%   savebinstack   - Save an binary stack to a file, as RGB Image.
%   savestack      - Save an image stack to a file or a serie of files
%
% Read/Write images in MetaImage format (used by ITK)
%   metaImageInfo  - Read information header of meta image data
%   metaImageRead  - Read an image in MetaImage format
%   metaImageWrite - Write header and data files of an image in MetaImage format
%
% Get information on 3D images
%   stackSize      - Compute the size of a 3D stack in [x y z] form
%   isColorStack   - Check if a 3D stack is color or gray-scale
%
% Manipulation of 3D images
%   stackSlice     - Extract a planar slice from a 3D image
%   rotateStack90  - Rotate a 3D image by 90 degrees around one image axis
%   stackRotate90  - Rotate a 3D image by 90 degrees around one image axis
%
% -----
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% created the  07/11/2005.
% Copyright INRA - Cepia Software Platform.
% http://www.pfl-cepia.inra.fr/index.php?page=slicer
% Licensed under the terms of the BSD License, see the file license.txt

%
