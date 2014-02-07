function mesh_plot(varargin)
global ux uy uz msize unormx unormy unormz unorm_mrms u_mrms x t
close all
optargin = size(varargin,2);
%set options based on varargin
slice=0 ; iso=0 ; spect=0 ; struct=0 ; print=0 ; para=0 ; 
vort_slice=0 ; vort_iso=0 ; smoothme=0 ; get_phi=0 ; 
double_spect=0 ;
for i=1:optargin
  switch cell2str(varargin(i))
    case 'slice'
      slice=1;
    case 'vort_slice'
      vort_slice=1;
    case 'iso'
      iso=1;
    case 'vort_iso'
      vort_iso=1;
    case 'spect'
      spect=1;
    case 'double_spect'
      double_spect=1;
    case 'para'
      para=1;
    case 'struct'
      struct=1;
    case 'smooth'
      smoothme=1;
    case 'get_phi'
      get_phi=1;
    case 'print'
      print=1;
      disp('printing to file')
  end
end
if smoothme==1
  ux=smooth3(ux,'box',3);
  uy=smooth3(uy,'box',3);
  uz=smooth3(uz,'box',3);
end
if get_phi==1
  mesh_get_phi(ux,uy,uz,msize,0) %don't plot
  mesh_get_phi(ux,uy,uz,msize,1) %do plot
end
%plot slices of field
if slice==1
  if u_mrms>0.
    mesh_slices(x,ux,uy,uz,msize,'super')
  end
  if unorm_mrms>0.
    mesh_slices(x,unormx,unormy,unormz,msize,'normal')
  end
end
if vort_slice==1
  if u_mrms>0.
    [curlx,curly,curlz,cav] = curl(ux,uy,uz) ;
    mesh_slices(x,curlx,curly,curlz,msize,'super-vorticity')
  end
  if unorm_mrms>0.
    [curlx,curly,curlz,cav] = curl(unormx,unormy,unormz) ;
    mesh_slices(x,curlx,curly,curlz,msize,'normal-vorticity')
  end
end
%output to paraview
if para==1
  if u_mrms>0.
    disp('printing to vtk file para_sup for paraview')
    writevtk(sqrt(ux.^2+uy.^2+uz.^2),'para_sup.vtk')
    savevtkvector(ux,uy,uz,'para_sup_vector.vtk');
  end
  if unorm_mrms>0.
    disp('printing to vtk file para_norm for paraview')
    writevtk(sqrt(unormx.^2+unormy.^2+unormz.^2),'para_norm.vtk')
  end
end
%plot isosurface
if iso==1
  if u_mrms>0.
    mesh_iso(x,ux,uy,uz,msize,'super')
  end
  if unorm_mrms>0.
    mesh_iso(x,unormx,unormy,unormz,msize,'normal')
  end
end
if vort_iso==1
  if u_mrms>0.
    [curlx,curly,curlz,cav] = curl(ux,uy,uz) ;
    mesh_iso(x,curlx,curly,curlz,msize,'super-vorticity')
  end
  if unorm_mrms>0.
    [curlx,curly,curlz,cav] = curl(unormx,unormy,unormz) ;
    mesh_iso(x,curlx,curly,curlz,msize,'normal-vorticity')
  end
end
%spectrum
if u_mrms>0.
  if spect==1
    mesh_spectrum(ux,uy,uz,msize,'super',0)
  end
  if struct==1
    %mesh_structure_func(x,ux,uy,uz,msize,'super')
  end
end
if unorm_mrms>0.
  if spect==1
    mesh_spectrum(unormx,unormy,unormz,msize,'normal',0)
  end
  if struct==1
    mesh_structure_func(x,unormx,unormy,unormz,msize,'normal')
  end
end
if double_spect==1
  mesh_double_spectrum(ux,uy,uz,unormx,unormy,unormz,msize)
end
  
