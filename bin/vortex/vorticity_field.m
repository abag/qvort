function vorticity_field(filenumber,mesh)
filename=sprintf('data/var%04d.log',filenumber);
%get the dimensions information from dims.log
dims=load('./data/dims.log');
if dims(4)==1
  fid=fopen(filename);
  if fid<0
      disp('var file does not exist, exiting script')
      return
  end
  time=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  x=fread(fid,number_of_particles,'float64');
  y=fread(fid,number_of_particles,'float64');
  z=fread(fid,number_of_particles,'float64');
  f=fread(fid,number_of_particles,'int');
  u=fread(fid,number_of_particles,'float64');
else 
  fid=fopen(filename);
  if fid<0
      disp('var file does not exist, exiting script')
      return
  end
  %read the time
  tline=fgetl(fid);
  dummy=textscan(tline, '%f');
  time=dummy{:};
  %how many particles
  tline=fgetl(fid);
  dummy=textscan(tline, '%d');  
  number_of_particles=dummy{:};
  %get the particles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j=1:number_of_particles
    tline=fgetl(fid);
    dummy=textscan(tline, '%f');
    dummy_vect=dummy{:};
    x(j)=dummy_vect(1);
    y(j)=dummy_vect(2);
    z(j)=dummy_vect(3);
    f(j)=dummy_vect(4);
    u(j)=dummy_vect(5);
  end
  f=uint16(f);
end
%what is the box size?
bsize=dims(2)
%divide this up into sections
for i=1:mesh+1
  bx(i)=-bsize/2+bsize*(i-1)/mesh;
  by(i)=-bsize/2+bsize*(i-1)/mesh;
  bz(i)=-bsize/2+bsize*(i-1)/mesh;
end
%now do a loop over the box#
vorticity_mesh(1:mesh,1:mesh,1:mesh,1:3)=0.;
for k=1:mesh  ; for j=1:mesh  ; for i=1:mesh 
  %must loop over particles
  for vi=1:number_of_particles
    if f(vi)==0
      continue
    end 
    %does particle lie within meshpoint?
    if (x(vi)>bx(i)) && (x(vi)<=bx(i+1))
      if (y(vi)>by(j)) && (y(vi)<=by(j+1))
        if (z(vi)>bz(k)) && (z(vi)<=bz(k+1))
          infront=f(vi);
          dist=sqrt((x(vi)-x(infront))^2+(y(vi)-y(infront))^2+(z(vi)-z(infront))^2);
          if (dist<.4*bsize)
            vorticity_mesh(k,j,i,1)=vorticity_mesh(k,j,i,1)+(x(infront)-x(vi));
            vorticity_mesh(k,j,i,2)=vorticity_mesh(k,j,i,2)+(y(infront)-y(vi));
            vorticity_mesh(k,j,i,3)=vorticity_mesh(k,j,i,3)+(z(infront)-z(vi));
          end
        end
       end
    end
  end
end ; end ; end 
vorticity_mesh=vorticity_mesh ; %*(mesh/bsize);
modv=sqrt(vorticity_mesh(:,:,:,1).^2+vorticity_mesh(:,:,:,2).^2+vorticity_mesh(:,:,:,3).^2) ;
if (0==1)
  p=patch(isosurface(modv));
  isonormals(modv, p)
  set(p, 'FaceColor', 'm', 'EdgeColor', 'none');
  daspect([1 1 1]); axis tight;
  camup([0 0 1 ]); campos([0.7686    0.1432    0.3043])
  camlight; lighting phong
  axis([0 mesh 0 mesh 0 mesh])
end
sum(sum(sum(vorticity_mesh)))
pcolor(squeeze(vorticity_mesh(:,mesh/2,:,3)));
%modv=sqrt(vorticity_mesh(:,:,:,1).^2+vorticity_mesh(:,:,:,2).^2+vorticity_mesh(:,:,:,3).^2) ;
%SliceBrowser(modv) ; colormap('jet')
%slicer(modv)
