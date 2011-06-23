function [] = writeVTK(vol,vtkfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: writeVTK(vol,vtkfile)
%
%   vol:     The 3D matrix to be saved to file
%   vtkfile: The output filename (string)
%   notes:   Only writes binary STRUCTURED_POINTS
%  
% Erik Vidholm 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dimensions
volinfo = whos('vol');

sz = volinfo.size;
X = sz(1); Y = sz(2); Z = 1;
if( length(sz) == 3 )
  Z = sz(3);
end

% open file (OBS! big endian format)
fid = fopen(vtkfile,'w','b');

% write header
fprintf(fid, '%s\n', '# vtk DataFile Version 3.0');
fprintf(fid, '%s\n', 'created by writeVTK ');
fprintf(fid, '%s\n', 'BINARY');  
fprintf(fid, '%s\n', 'DATASET STRUCTURED_POINTS');  
fprintf(fid, '%s%d%c%d%c%d\n', 'DIMENSIONS ', X, ' ', Y, ' ', Z);
fprintf(fid, '%s%f%c%f%c%f\n', 'ORIGIN ', 0.0, ' ', 0.0, ' ', 0.0); 
fprintf(fid, '%s%f%c%f%c%f\n', 'SPACING ', 1.0, ' ', 1.0, ' ', 1.0); 
fprintf(fid, '%s%d\n', 'POINT_DATA ', X*Y*Z);

tp = volinfo.class;
if( strcmp(tp, 'uint8') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data unsigned_char');
elseif( strcmp(tp, 'int8') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data char');
elseif( strcmp(tp, 'uint16') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data unsigned_short');
elseif( strcmp(tp, 'int16') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data short');
elseif( strcmp(tp, 'uint32') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data unsigned_int');
elseif( strcmp(tp, 'int32') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data int');
elseif( strcmp(tp, 'single') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data float');
elseif( strcmp(tp, 'double') > 0 )
  fprintf(fid, '%s\n', 'SCALARS image_data double');
end

fprintf(fid, '%s\n', 'LOOKUP_TABLE default');

% write data as binary
fwrite(fid,vol,tp);

% close file
fclose(fid);
