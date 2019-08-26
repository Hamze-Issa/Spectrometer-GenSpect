function [lins,truncation_info]=hitloadr(param_filename,gasname,iso,start_wavnum,stop_wavnum,disp_log, gasNames)

% HITLOADR(param_filename,gasname,iso,start_wavnum,stop_wavnum,disp_log) 
% loads the -mat hitran format data file for gas GASNAME.  These files will have been
% saved in -mat form, with names such as genhit96_05.lin and must be
% included in the MATLAB path. They contain the structure variable LINS.
%
% Inputs:
% PARAM_FILENAME is the HITRAN format line database filename
% GASNAME is the HITRAN gasname (string or cell array)
% ISO is an array of isotopes being considered (column vector)
% START_WAVNUM is the lower boundary of the calculation in wavenumbers (scalar) 
% STOP_WAVNUM is the upper boundary of the calculation in wavenumbers (scalar) 
% DISP_LOG is set to 'off' to supress all non-error output (for www version)
%
% Outputs:
% LINS is a structure with fields gas_file, gasname, run_time, iso, wavenum, intens, 
%			air_broad, self_broad, E_lower, temp_coeff, press_shift, transition, grid_index
%
% (C) Ben Quine, 29-Feb-2000 Rev 1.01.

if exist('disp_log')==1
   if strcmp(lower(disp_log),'off')
      disp_log=0;
   elseif strcmp(lower(disp_log),'www')
      disp_log=2;
   else
      disp_log=1;
   end
else
   disp_log=1;
end
if disp_log==1
   disp(['GENSPECT HITLOADR:loading parameter data for ' gasname])
end
% load line data stored in lins
lins = getGasFile(param_filename, gasNames);
% load(param_filename, '-mat');
% Record Gas File
lins.gas_file = param_filename;
lins.gasname = gasname;
lins.run_time = datestr(now);
%truncate the lins array to the desired range [Wavenumber and Isotope]
if iso==0		% All Isotopes
   trunc_vec=find(lins.wavenum>start_wavnum & lins.wavenum<stop_wavnum);
else				% Specific Isotope
   trunc_vec=find(lins.wavenum>start_wavnum & lins.wavenum<stop_wavnum & lins.iso==iso);
end
% truncate line database if www calculation is to be performed
max_line_limit = 5e4;
if disp_log==2 & length(trunc_vec)>max_line_limit
      trunc_vec=trunc_vec(1:max_line_limit);
      truncation_info = [lins.gasname ' line database truncated at ' num2str(lins.wavenum(trunc_vec(end))) ' cm<SUP>-1</SUP>. '];
else
      truncation_info='';
end
lins.iso=lins.iso(trunc_vec);
lins.wavenum=lins.wavenum(trunc_vec);
lins.intens=lins.intens(trunc_vec);
lins.air_broad=lins.air_broad(trunc_vec);
lins.self_broad=lins.self_broad(trunc_vec);
lins.E_lower=lins.E_lower(trunc_vec);
lins.temp_coeff=lins.temp_coeff(trunc_vec);
lins.press_shift=lins.press_shift(trunc_vec);
lins.transition=lins.transition(trunc_vec);
lins.local_q_low=lins.local_q_low(trunc_vec);
lins.local_q_up=lins.local_q_up(trunc_vec);
lins.n_lines = length(trunc_vec);