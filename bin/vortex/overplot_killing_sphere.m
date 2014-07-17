function overplot_killing_sphere(var_file)
hold on
A=load('./data/adaptive_killing_radius.log');
centre=A(var_file,1:3);
rad=A(var_file,4);
[x y z] = sphere(128);
c(1:129,1:129)=100;
x=(x+centre(1))*rad;
y=(y+centre(2))*rad;
z=(z+centre(3))*rad;
h = surf(x, y, z,c); 
set(h, 'FaceAlpha', 0.3)
shading interp
camlight; lighting phong
colormap(fireprint)
