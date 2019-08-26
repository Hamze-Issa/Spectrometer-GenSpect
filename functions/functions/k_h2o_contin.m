function Kc=k_h20_contin(cells,cells_to_calc,gas_index)

% Kc=k_h20_contin(cells,cells_to_calc,h20_index)
%
% This function computes the values of the absorption 
% coefficent k for the water vapour continuum for the
% 700-1200cm-1 window. Functional form of continuum taken
% from Roberts et al. Applied Optics Vol.15 No 9, 1976
%
%
% The function forms part of the GENSPEC suite which
% defines gencalc and cells. Only an emprical fit for the
% self-broadening term is considered for the moment.
%
%
% Ben Quine, U of T, 24-FEB-1999, Revised Jan 2001

%Empirical results are taken from
% D. E. Burch and R. L. Alt: Continuum Absorption by H20
% in the 700-1200cm-1 and 2400-2800cm-1 windows, AFGL, 1984.
% Raw self-broading data at 284 and 296, 700:20:1100
% note 1080 to 1160 are extrapolated
%C.wavnum=[700:20:1160];
%C.self(1,:) = [775 710 650 600 550 500 460 420 385 355 319 292 268 246 224 205 187 170 153 153-17 153-2*17 153-3*17 153-4*17 153-5*17]*1e-24;
%C.self(2,:) = [620 550 485 435 390 345 312 284 260 240 219 199 183 168 157 145 133 122 112 112-10 112-2*10 112-3*10 112-4*10 112-5*10]*1e-24;
%C.self_temp = [284; 296];
%% Raw foreign-broading data at 296K, 700:20:1100
%% 1020:20:1160 are extrapolated [using LOWTRAN6 data as a guide]
%C.foreign_temp = [296];
%C.foreign(1,:) = [285 242 200 162 127 100 79 63 51 43 38 33 32 30 29 28 28-1 28-1 28-0 28+1 29+1 30 30 30]*1e-26;
% find wavnumber points that fall within data limits
%ind = find(gencalc.calc_grid>min(C.wavnum) & gencalc.calc_grid<max(C.wavnum));
% interpolate data onto frequency grid
%C.interp_self = interp1(C.wavnum,C.self,gencalc.calc_grid(ind));
%C.interp_foreign = interp1(C.wavnum,C.foreign,gencalc.calc_grid(ind));

a = 1.25e-22;			% mol-1cm2atm-1
b = 2.34e-19;			% mol-1cm2atm-1
beta = 8.3e-3;			% cm
C_self = a + b*exp(-beta*cells.wavnum);

T0 = 1800;				% K
temp_dep = exp(T0*(1./cells.tmean(cells_to_calc)-1/296));
p0=1e5;
gamma = 0.005;
% pressure dependence [self only]
%press_dep = cells.ppartial(:,gas_index)./p0;
% pressure dependence [self + foriegn]
press_dep = (1-0.001)*cells.ppartial(cells_to_calc,gas_index)./p0 + 0.0008*cells.pmean(cells_to_calc)/p0;

Kc = C_self * (press_dep'.*temp_dep');
