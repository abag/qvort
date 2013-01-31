function small_loop_tracker(lnumber)
A=load('./data/small_loops.log');
var_snap=ceil(A(lnumber,1));
figure('visible','off');
for i=-15:2
    if var_snap+i<1 ; continue ;end
    i
    vortex_plot(var_snap+i,'thin_line');
    dims=load('./data/dims.log');
    zoom_fac=50;
    axis([A(lnumber,4)-dims(2)/zoom_fac A(lnumber,4)+dims(2)/zoom_fac A(lnumber,5)-dims(2)/zoom_fac A(lnumber,5)+dims(2)/zoom_fac A(lnumber,6)-dims(2)/zoom_fac A(lnumber,6)+dims(2)/zoom_fac])
    box off
    axis off
    fout=sprintf('./data/loop_form%02d.jpeg',i+16);
    print('-djpeg',fout)
end
figure('visible','on');