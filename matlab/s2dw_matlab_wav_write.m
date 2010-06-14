function s2dw_matlab_wav_write(J, N, B, bl_scoeff, alpha, scoeff, ...
                              wcoeff, filename, comment)

% function s2dw_matlab_wav_write(J, N, B, bl_scoeff, alpha, scoeff, ..
%                               wcoeff, filename, comment)
% 
%  S2DW Matlab code
%  Written by Jason McEwen (mcewen@mrao.cam.ac.uk)
%
%  Write S2DW data manipulated in Matlab to a S2DW .m file and
%  corresponding data .dat files that can be read by Matlab and
%  the S2DW Fortran code.
%
%  Variables:
%   J: Maximum analysis scale depth.
%   B: Harmonic band limit.
%   N: Azimuthal band limit.
%   bl_scoeff: Upper band limit for scaling coefficients.
%   alpha: Basis dilation factor.
%   scoeff(1:bl_scoeff,1:bl_scoeff): Scaling coefficients.
%   wcoeff{1:J+1}(aa,bb,gg): Wavelet coefficients (ranges of aa, bb 
%     and gg depend on scale).

% Open file.
fileid = fopen(filename,'w');

% Write header.
fprintf(fileid, '%s%s\n', 'function [J, N, B, bl_scoeff, alpha, scoeff, wcoeff] = ', ...
        filename(1:end-2) );
fprintf(fileid, '%s%s\n', '% function [J, N, B, bl_scoeff, alpha, scoeff, wcoeff] = ', ...
        filename(1:end-2) );
fprintf(fileid, '%s\n', '%');

fprintf(fileid, '%s\n', '%  Scale discretised wavelet coefficients created by S2DW code');
fprintf(fileid, '%s\n', '%  S2DW code written by Jason McEwen (mcewen@mrao.cam.ac.uk)');
fprintf(fileid, '%s\n', '%');
fprintf(fileid, '%s\n', '%  This file written by: Matlab');
if(nargin > 8) then
  fprintf(fileid, '%s%s\n', '%  Comment: ', comment);
else
  fprintf(fileid, '%s\n', '%  Comment: ');
end
fprintf(fileid, '%s\n', '%');
fprintf(fileid, '%s\n', '%  Variables:');
fprintf(fileid, '%s\n', '%   J: Maximum analysis scale depth.');
fprintf(fileid, '%s\n', '%   B: Harmonic band limit.');
fprintf(fileid, '%s\n', '%   N: Azimuthal band limit.');
fprintf(fileid, '%s\n', '%   bl_scoeff: Upper band limit for scaling coefficients.');
fprintf(fileid, '%s\n', '%   alpha: Basis dilation factor.');
fprintf(fileid, '%s\n', '%   scoeff(1:bl_scoeff,1:bl_scoeff): Scaling coefficients.');
fprintf(fileid, '%s%s\n', '%   wcoeff{1:J+1}(aa,bb,gg): Wavelet coefficients ', ...        
        '(ranges of aa, bb and gg, depend on scale).');
fprintf(fileid, '\n');

% Write parameters.
fprintf(fileid, '%s\n', '% Size parameters');
fprintf(fileid, '%s%27d%s\n', 'J         = ', J, ';');
fprintf(fileid, '%s%27d%s\n', 'N         = ', N, ';');
fprintf(fileid, '%s%27d%s\n', 'B         = ', B, ';');
fprintf(fileid, '%s%27d%s\n', 'bl_scoeff = ', bl_scoeff, ';');
fprintf(fileid, '%s%27.20e%s\n', 'alpha     = ', alpha, ';');
fprintf(fileid, '\n');

% Write wavelet sizes.
fprintf(fileid, '%s\n', '% Wavelet coefficient array sizes');
for jj = [0:J]
  fprintf(fileid,'%s%5i%s%15i%s\n', 'wcoeff_size_aa(', jj+1, ') = ', ...
          size(wcoeff{jj+1},1), ';');
  fprintf(fileid,'%s%5i%s%15i%s\n', 'wcoeff_size_bb(', jj+1, ') = ', ...
          size(wcoeff{jj+1},2), ';');
  fprintf(fileid,'%s%5i%s%15i%s\n', 'wcoeff_size_gg(', jj+1, ') = ', ...
          size(wcoeff{jj+1},3), ';');
end
fprintf(fileid, '\n'); 

% Fast way to write scaling coefficients to data file.
filename_scoeff = sprintf('%s%s', filename(1:end-2), '_scoeff.dat');
scoeff_array = [];
for el = [0: bl_scoeff-1]
  for m = [0: el]
    scoeff_array = [scoeff_array; real(scoeff(el+1, m+1))];
  end
end
for el = [0: bl_scoeff-1]
  for m = [0: el]
    scoeff_array = [scoeff_array; imag(scoeff(el+1, m+1))];
  end
end
save(filename_scoeff, 'scoeff_array', '-ASCII', '-DOUBLE');

% Write matlab code to read scaling coefficients.
fprintf(fileid, '%s\n', '% Scaling coefficients');
fprintf(fileid, '%s%s%s\n', 'scoeff_array = load(''', filename_scoeff, ''');');
fprintf(fileid, '%s\n', 'n = length(scoeff_array) / 2;');
fprintf(fileid, '%s\n', 'scoeff_array = scoeff_array(1:n) + i*scoeff_array(n+1:end);');
fprintf(fileid, '%s\n', 'scoeff = zeros(bl_scoeff,bl_scoeff);');
fprintf(fileid, '%s\n', 'ind_start = 1;');
fprintf(fileid, '%s\n', 'for el = 0:bl_scoeff-1');
fprintf(fileid, '%s\n', '    ind_end = ind_start - 1 + el + 1;');
fprintf(fileid, '%s\n', '    scoeff(el+1, 1:el+1) = scoeff_array(ind_start:ind_end).'';');
fprintf(fileid, '%s\n', '    ind_start = ind_end + 1;');
fprintf(fileid, '%s\n', 'end');
fprintf(fileid, '\n'); 

% Fast way to write wavelet coefficients to data file.
filename_wav = sprintf('%s%s', filename(1:end-2), '_wcoeff.dat');
wcoeff_array = [];
for jj=1:J+1
    wcoeff_array = [wcoeff_array; wcoeff{jj}(:)];
end
save(filename_wav, 'wcoeff_array', '-ASCII', '-DOUBLE');

% Write matlab code to read wavelet coefficients.
fprintf(fileid, '%s\n', '% Wavelet coefficients');
fprintf(fileid, '%s%s%s\n', 'wcoeff_array = load(''', filename_wav, ''');');
fprintf(fileid, '%s\n', 'ind_start = 1;');
fprintf(fileid, '%s\n', 'for jj = 1:J+1');
fprintf(fileid, '%s%s\n', '    ind_end = ind_start - 1 +', ...
        'wcoeff_size_aa(jj) * wcoeff_size_bb(jj) * wcoeff_size_gg(jj);');
fprintf(fileid, '%s%s\n', '    wcoeff{jj} = reshape(wcoeff_array(ind_start:ind_end),', ...
    'wcoeff_size_aa(jj), wcoeff_size_bb(jj), wcoeff_size_gg(jj));');
fprintf(fileid, '%s\n', '    ind_start = ind_end+1;');
fprintf(fileid, '%s\n', 'end');
fprintf(fileid, '\n'); 

% Close file.
fclose(fileid);

