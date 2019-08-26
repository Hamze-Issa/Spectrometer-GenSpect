function qtweight=qtweight(mol,Tout)

% qt=QTWEIGHT(Mol,Tout) looks up the molecular weight and Q value for molecule
% number Mol (in the hitran database) for all included isotopes at temperature Tout.
%
% MOL is the number of the gas as referenced in the HITRAN database (scalar)
% TOUT is the temperature of the gas (scalar)
%
% Internal partition function is calculated as specified in HITRAN96. This
% function uses a four parameter estimate over two temperature ranges [70-2000K].
% The function takes Molicule, Isotope number and temperature inputs and estimates
% qt.Q_T (Q(T)) and qt.Q296 (Q(296K)) for each line. qt.ratio is Q(296)/Q(T)
%
% This function read coeficient data from QT_DAT.MAT located in the matlab search
% path. The user can modify their own TIP parameters directly or use the utility
% QTREAD to read in an ascii script such as QT.DAT
%
% (C) Ben Quine 1999.
%

% keep data global to this function
global qt
% load data is not already there
if isempty(qt)
   load qt_data.mat
   disp('GENSPECT QTWEIGHT: Loading QT data')
end
% checkout dataset
if max(mol)>qt.nmol
   error(['GENSPECT QTWEIGHT: Molicule number of ' num2str(max(mol)) ' exceeds limits'])
end
% form index for lookup
index = qt.mol_index(mol)+[1:qt.ison(mol)]';
% determine temperature range
if Tout>qt.tmax | Tout<qt.tmin
   error(['Temperature of ' num2str(Tout) ' exceeds range limits']);
end
if length(Tout)~=1
   error('GENSPECT QTWEIGHT: Temperature must be scalar')
end
trange=qt.n_coeff;
for i=qt.n_coeff-1:-1:1
   if Tout<qt.t_limits(i)
      trange=i;
   end
end
%return Q296
qtweight.Q296 = qt.Q296(index)';
% return Qt_T
n_terms = length(qt.coeff(1,1,:));
qtweight.Q_T = reshape(qt.coeff(trange,index,1:end),length(index),n_terms);
for i=2:n_terms
  qtweight.Q_T(:,i)=qtweight.Q_T(:,i)*(Tout.^(i-1));
end
% sum power series terms
qtweight.Q_T=sum(qtweight.Q_T,2);

% set entries with no supplied temperature variation to Q296
ind=find(qtweight.Q_T==-1)';
qtweight.Q_T(ind)=qtweight.Q296(ind);

% return ratio
qtweight.ratio = qtweight.Q296./qtweight.Q_T;
