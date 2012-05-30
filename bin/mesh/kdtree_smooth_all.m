function kdtree_smooth_all(start,finish,skip,mesh_size,rms)
global wmesh
if nargin<3
  disp('I do not have enough arguements to run, please supply start, final, skip')
  return
elseif nargin==3
  mesh_size=32; 
  rms=2.;
elseif nargin==4
  rms=2.
end
counter=1;
for i=start:skip:finish
  [L_inside L_outside L_total]=kdtree_smooth(i,mesh_size,0,rms)
end

