function rotate_snap(snap_count)
for i=1:snap_count
    campan(360/snap_count,1)
    drawnow
    fout=sprintf('data/rot_snap%03d',i);
    print ('-dpng',fout)
end