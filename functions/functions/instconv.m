function k_inst=instconv(k,n_wavnum,calc_int,instfn)
%
% Function to convolve a spectral calculation with an instrument function.
%
% function k_inst=instconv(k,n_wavnum,calc_int,instfn)
%
% k is the input intensity function vector
% k_inst is the output intensity function vector (dimensions as k)
% n_wavnum is the scalar 1 sigma point (width of the slit function) in wavenumbers
% calc_int is the scalar calculation grid interval (wavenumbers cm-1)
% instfn is an optional parameter specifying the type of slit function
%
% Currently a 2 sigma truncated gaussian is implemented called 'gaus'.
% Experimentally, a sine^2 instrument function is also implemented 'sin2'.
%
% (C) Ben Quine 07-SEP-1999, University of Toronto.

% Revision: Ben Quine 11-NOV-1998
% Notes:		Made function use truncated gaussian function.
%				Revised code to convolve a matrix of k vectors.
%				Variables passed to function are now in wavenumber space.
% Revision: Ben Quine 19-JAN-1999
% Notes:		Improved help to include dimensions and input units

if ~exist('instfn')
   instfn = 'gaus';
else
   instfn = [instfn '    '];
end
n = (n_wavnum./calc_int);
% force n odd
n_points = floor(n/2)*2 + 1;

% choose function
if lower(instfn(1:4)) == 'sin2'
	% calc instrument function
	fn=(sin(linspace(0,pi,2*n_points)).^2);
else   
	% calc truncated gaussian instrument function
	x = linspace(-3*n_points,3*n_points,6*n_points);
	fn = exp(-0.5*(x/n).^2)./(n*(2*pi)^0.5);
end
% normalise (setting boundary values to zero)
fn = (fn-fn(1)) ./ sum(fn-fn(1));

k_inst = zeros(size(k,1)+length(fn)-1,size(k,2));
for i=1:size(k,2)
	% convolve the slit function with the spectral calc
	k_inst(:,i) = conv(k(:,i),fn);
end
% truncate the convolution
%k_inst = k_inst(n_points+1:length(k)+n_points);
k_inst = k_inst(3*n_points+1:length(k)+3*n_points,:);
