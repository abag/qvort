function kdtree_smooth_all(start,finish,skip,mesh_size,rms)
global LL_inside LL_outside Wrms
if nargin<3
  disp('I do not have enough arguements to run, please supply start, final, skip')
  return
elseif nargin==3
  mesh_size=32; 
  rms=1.5;
elseif nargin==4
  rms=2.
end
counter=1;
for i=start:skip:finish
  i
  %[L_inside L_outside L_total rms_w]=kdtree_smooth(i,mesh_size,0,rms);
  [L_inside L_outside L_total rms_w]=kdtree_vort_smooth(i,0,rms);
  LL_inside(counter,:)=L_inside;
  LL_outside(counter,:)=L_outside;
  LL_total(counter)=L_total;
  Wrms(counter)=rms_w;
  counter=counter+1;
end

