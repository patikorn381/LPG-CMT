function write_csv(fname, header, M)
% WRITE_CSV  Write a numeric matrix with a header line (MATLAB and Octave).
% Replaces dlmwrite, which is not recommended in recent MATLAB releases.
fid = fopen(fname, 'w');
if fid < 0, error('write_csv: cannot open %s', fname); end
fprintf(fid, '%s\n', header);
fmt = [repmat('%.10g,', 1, size(M,2)-1) '%.10g\n'];
fprintf(fid, fmt, M.');
fclose(fid);
end
