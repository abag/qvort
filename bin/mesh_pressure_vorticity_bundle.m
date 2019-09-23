function [ux uy uz wx wy wz]=mesh_pressure_vorticity_jason(filenumbers, downsample,display)


sigma= 1;
skip = 8;
mesh_out= 1;
iso_vort = 0.3;
iso_press = 0.3;




if nargin<2
    downsample=0;
    display=0;figure
end
if nargin<3
    display=0;
else
    display=1;
    disp('producing iso surfaces before and after + paraview output')
end
load ./data/dims.log;
msize=dims(3);
if (msize==0) 
  disp('mesh size is zero exiting script')
  return
end
disp(sprintf('mesh size is: %04d',msize))
counter=0;
for ii=filenumbers
  filename=sprintf('data/mesh%03d.dat',ii);
  fid=fopen(filename);
  if fid<0
    disp('mesh file does not exist, exiting script')
    return
  end
  t=fread(fid,1,'float64');
  x=fread(fid,msize,'float64');
  
  
  unormx=fread(fid,msize^3,'float64');
  unormy=fread(fid,msize^3,'float64');
  unormz=fread(fid,msize^3,'float64');
  ux=fread(fid,msize^3,'float64');
  uy=fread(fid,msize^3,'float64');
  uz=fread(fid,msize^3,'float64');
  ux=reshape(ux,msize,msize,msize);
  uy=reshape(uy,msize,msize,msize);
  uz=reshape(uz,msize,msize,msize);

  ux=permute(ux, [2 3 1]);
  uy=permute(uy, [2 3 1]);
  uz=permute(uz, [2 3 1]);
 
  
 ux = smooth_gauss(ux,sigma);
 uy = smooth_gauss(uy,sigma);
 uz = smooth_gauss(uz,sigma);
  
  
  u_rms= sqrt(sum(ux(:).^2 + uy(:).^2 + uz(:).^2)/(msize^3.0))


  ux_hat=fftn(ux);
  uy_hat=fftn(uy);
  uz_hat=fftn(uz);

  k_vec=(2.0*pi/range(x))*[0:(msize/2)-1 0 (-(msize/2)+1):-1];


  
  for i=1:msize
  	ux_hat_x(:,i,:) = 1i*k_vec(i)*ux_hat(:,i,:);
  	ux_hat_y(i,:,:) = 1i*k_vec(i)*ux_hat(i,:,:);
  	ux_hat_z(:,:,i) = 1i*k_vec(i)*ux_hat(:,:,i);
	
  	uy_hat_x(:,i,:) = 1i*k_vec(i)*uy_hat(:,i,:);
  	uy_hat_y(i,:,:) = 1i*k_vec(i)*uy_hat(i,:,:);
  	uy_hat_z(:,:,i) = 1i*k_vec(i)*uy_hat(:,:,i);
	
  	uz_hat_x(:,i,:) = 1i*k_vec(i)*uz_hat(:,i,:);
  	uz_hat_y(i,:,:) = 1i*k_vec(i)*uz_hat(i,:,:);
  	uz_hat_z(:,:,i) = 1i*k_vec(i)*uz_hat(:,:,i);
  end

  ux_x = ifftn(ux_hat_x);
  ux_y = ifftn(ux_hat_y);
  ux_z = ifftn(ux_hat_z);

  uy_x = ifftn(uy_hat_x);
  uy_y = ifftn(uy_hat_y);
  uy_z = ifftn(uy_hat_z);

  uz_x = ifftn(uz_hat_x);
  uz_y = ifftn(uz_hat_y);
  uz_z = ifftn(uz_hat_z);
  
   
  press = ux_x.^2 + uy_y.^2 + uz_z.^2 + 2.0*(ux_y.*uy_x + ux_z.*uz_x + uy_z.*uz_y);
  
  wx=uz_y-uy_z;
  wy=ux_z-uz_x;
  wz=uy_x-ux_y;


  vort= sqrt((wx).^2 + (wy).^2 + (wz).^2);

  press_hat = fftn(press);

 
  for i=1:msize
  	for j=1:msize
  		for k=1:msize
  			k2 = k_vec(i)^2 + k_vec(j)^2 + k_vec(k)^2;
  			if k2 ==0
  				k2=1;
  				press_hat(i,j,k)=0;
			end
  			press_hat(i,j,k) = press_hat(i,j,k) / k2 ;
		end
  	end
  end
  press = ifftn(press_hat);
  
  
  segment_x = linspace(-0.05,0.05,msize)

   for i=1:msize
  	segment_press(i) = press(i,msize/2 -10, msize/2);
  	segment_vort(i) = vort(i,msize/2 -10,msize/2) ;
  end

end


plot(segment_x, segment_press,'b','LineWidth',3)
hold on
plot(segment_x, segment_vort,'r','LineWidth',3)
  %set(h_legend, 'interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast')
xlab = xlabel('Distance', 'interpreter', 'latex')
set(xlab, 'FontSize', 14);
ylab=ylabel('Pressure', 'interpreter', 'latex');
set(ylab, 'FontSize', 14);

print(gcf,'press_profile','-depsc','-r600')



x=linspace(-0.05,0.05,128);
y=linspace(-0.05,0.05,128);
z=linspace(-0.05,0.05,128);


figure
v = vort;
v = v/max((abs(v(:))));
hold on
vortex_plot(mesh_out*filenumbers)
patch(isocaps(x,y,z,v,iso_vort),'FaceColor','interp','EdgeColor','none');
p1 = patch(isosurface(x,y,z,v,iso_vort),'FaceColor','blue','EdgeColor','none');
isonormals(x,y,z,v,p1);
view(3);
alpha('0.5')
axis vis3d tight
camlight left
colormap('jet');
lighting gouraud
%xlim([ ]);
%ylim([ ]);
%zlim([ ]);
hold off
print(gcf,'vorticity_iso','-depsc','-r600')


figure
v = -press;
v = v/max((abs(v(:))));

hold on
patch(isocaps(x,y,z,v,iso_press),'FaceColor','interp','EdgeColor','none');
p1 = patch(isosurface(x,y,z,v,iso_press),'FaceColor','blue','EdgeColor','none');
isonormals(x,y,z,v,p1);
vortex_plot(mesh_out*filenumbers)
hold off
view(3);
alpha('0.5')
axis vis3d tight
camlight left
colormap('jet');
lighting gouraud
%xlim([ ]);
%ylim([ ]);
%zlim([ ]);
print(gcf,'pressure_iso', '-depsc','-r600')
end

function q = smooth_gauss(q,h)

N=size(q);
q_hat = fftn(q);
k=[0:(N/2 -1) 0 (1-N/2):-1];
[kx,ky,kz]= meshgrid(k,k,k);
filter = exp(-(0.5/h^2)*(kx.^2+ky.^2+kz.^2).^2);
q_hat = q_hat.*filter;
q=ifftn(q_hat);


end




