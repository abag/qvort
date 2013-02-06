function povray_out
global dims
global x y z
global f u u2
global number_of_particles
delete('test.log')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:number_of_particles
  if round(f(j))==0
  else
    dummy_x(1,1)=x(j);
    dummy_x(2,1)=x(round(f(j)));
    dummy_x(1,2)=y(j);
    dummy_x(2,2)=y(round(f(j)));
    dummy_x(1,3)=z(j);
    dummy_x(2,3)=z(round(f(j)));
    dist=sqrt((dummy_x(1,1)-dummy_x(2,1))^2+(dummy_x(1,2)-dummy_x(2,2))^2+(dummy_x(1,3)-dummy_x(2,3))^2);
    if dist<dims(2)/2.
        dummy_x=dummy_x;
      stra=strcat('<',num2str(dummy_x(1,1)),',',num2str(dummy_x(1,2)),',',num2str(dummy_x(1,3)),'>');
      strb=strcat('<',num2str(dummy_x(2,1)),',',num2str(dummy_x(2,2)),',',num2str(dummy_x(2,3)),'>');
      A=strcat(stra,',','     ',strb,',','\n');
      fid=fopen('test.log','a');
      fprintf(fid,A);
      fclose(fid);
    end
   end
end
system('sed -ie ''$s/.$//'' test.log')
