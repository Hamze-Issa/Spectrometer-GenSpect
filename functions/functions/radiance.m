function paths = radiance(paths,sub_layers,calc_type)

% paths = RADIANCE(paths,sub_layers,calc_type)
%
% function to calculate radiance through a path
% Use path_gas, path_source and path_reflect to form a path 
% This function generates the radiance after each path segment
%
% PATHS     is a structure specifying the radiation path.
% SUBLAYERS is the number of sub-layers to include during the calculation 
%           of radiance (settings other than 1 should be used with caution)
% CALC_TYPE is set to 'by segment' (default) to calculate radiance after
%           each path segment or to 'by path' to calculate only the radiance
%           for the total path.
%
% Radiances are in units [Wm^{-2}sr^{-1}(cm^{-1})^{-1}].
%
% (C) Ben Quine 02-OCT-2001.
%

n_k = length(paths.wavnum);					% number of k-vector elements

% interpolate combination onto finer grid
if sub_layers>1
   warning('GENSPECT RADIANCE: Sub-layer settings other than 1 should be used with extreme caution')
   % form distance grid, interpolated distance grid and weight factor
   j=1;
   for i=1:paths.n_elements
      if paths.type_index==1 & i~=1		% gas cell case (excluding first)
         dist(i)=paths.thickness(i);
         dist_interp(j:j+sub_layers-1) = ones(sub_layers,1)*paths.thickness(i)/sub_layers;
         weight_interp(j:j+sub_layers-1) = ones(sub_layers,1)./sub_layers;
         j=j+sub_layers;
      elseif i~=1							% case where not a gas cell
            dist(i)=paths.thickness(i);
            dist_interp(j) = paths.thickness(i);
            weight_interp(j) = 1;
         j=j+1;
      else										% start case (first layer not interpolated)
         dist(i)=0;
   		weight_interp(j)=1;
         dist_interp(j)=0;
         j=j+1;
      end
   end
   dist=cumsum(dist);											% form monotonic dist axis
   dist_interp=min(cumsum(dist_interp),max(dist));		% form monotonic interpolated distance axis
   n_layers = length(dist_interp);
   % iterpolate log_transmission onto finer grid
   log_trans = interp1(dist,paths.log_transmission',dist_interp)'.*(ones(n_k,1)*weight_interp);
   % interpolate temperatures onto same latitude grid
   temps = interp1(dist,paths.tmean,dist_interp);
   % iterpolate emissivities onto finer grid
   emissivity = interp1(dist,paths.emissivity',dist_interp)'.*(ones(n_k,1)*weight_interp);
   % calculate emission for each sublayer
   paths.emission = emissivity.*planck(paths.wavnum,temps);
	% compute the cumalative top-down transmissions (starting from top)
   cum_log_trans = cumsum(log_trans,2);
else
   % number of layers
   n_layers = paths.n_elements;
    % calculate emission for all layers
    paths.emission = paths.emissivity.*planck(paths.wavnum,paths.tmean);
	% compute the cumalative top-down transmissions (starting from top)
	cum_log_trans = cumsum(paths.log_transmission,2);
end
% If there are radiant sources then add them here
ind=find(paths.type_index==3);
if ~isempty(ind)
    paths.emission(:,ind)=paths.source_radiance(:,ind);
end
% include the angular size of source terms
paths.emission=paths.emission.*(ones(length(paths.wavnum),1)*paths.ang_size);

% determine the type of calculation required and restrict accordingly
if exist('calc_type')~=1
   calc_type='by segment';
end
if strcmp(lower(calc_type),'by path')==1
   start_layer = n_layers;
   radiance = zeros(n_k,1);
elseif strcmp(lower(calc_type),'by segment')==1
   start_layer = 1;
	radiance = zeros(n_k,n_layers);
else 
   warning(['GENSPECT RADIANCE: Unrecognised calculation type - ' upper(setstr(calc_type)) ' Default to ''by segment'''])
    start_layer = 1;
    radiance = zeros(n_k,n_layers);
end

% calculation index for numbering radiance calculations
j=0;
% compute radiances for layer i [nth layer backwards]
for i=start_layer:n_layers
   j=j+1;
   % compute downward transmission for all the layers above [assert 0 as lowest transmission to avoid LSB errors]
   transmission = max(exp(cum_log_trans(:,i)*ones(1,i)-cum_log_trans(:,1:i)),0);
   % sum emission contributions transmitted from layer above
   radiance(:,j) = radiance(:,j)+sum(transmission.*paths.emission(:,1:i),2);
end
% extract layer radiance from sublayer calculation
if sub_layers>1 & start_layer~=n_layers
	paths.radiance = interp1(dist_interp,radiance',dist)';
else
   paths.radiance = radiance;
end
% calculate transmission
if start_layer==n_layers
	paths.transmission=max(exp(cum_log_trans(:,end)),0);
else		% show all layer transmissions
   paths.transmission=max(exp(cum_log_trans),0);
end
