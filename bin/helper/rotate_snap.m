function rotate_snap(snap_count)
for i=1:snap_count
    [Az, El] = view;
    El       = ceil(El);
    Az_new = Az+360*(i-1)/snap_count;
    view([Az_new, El]);
    drawnow
    fout=sprintf('data/rot_snap%03d',i);
    print ('-dpng',fout)
end
