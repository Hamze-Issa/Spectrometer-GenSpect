function g=voigtwells(x,y,p)

% g=VOIGTWELLS(x,y,p) generates a voigt profile using the complex 
%   error function approximation computed by region:
%
%  Function is from R. J. Wells, Rapid Approximation of the Voigt Function
%								and its derivatives, J.Quant Spec & Rad Trans 62 (1999) 29-48.
%
%
%  g(x,y) = y/pi * INT[exp(-t^2)/((x-y)^2 + y^2),-inf,inf,dt]
%
% x may be a row vector
% y must be scalar
% p is an accuracy parameter. Calculation error less than 10^(-p) [2<=p=<6]
%
% VOIGTWELLS calls no other functions
% Vector calcalutions must be column formatted
%
% (c) Ben Quine, 2000.

if ~exist('p')
   p=4;
end
if size(x,1)~=1 | length(y)~=1
   disp('GENSPECT VOIGTWELLS - (x,y) parameters are incorrectly dimensioned.')
   return
end
yq = y.*y;
abx=abs(x);
xq = abx.*abx;

% configure region boundary 0,1,3
if p==4
	xlim0 = 15100.0 + y*(40.0 + 3.6*y);
   xlim1 = 164.0 - y*(4.3 + y*1.8);
	xlim3 = 5.76*yq;
else
	R0 = 1.51*exp(1.144*p);
	R0 = min(max(14.88,R0),460.4);
	R1 = 1.6*exp(0.554*p);
	R1 = min(max(4.85,R1),25.5);
   xlim0 = R0-y;
   xlim1 = R1-y;
   xlim3 = 3.097*y -0.45;
end   
% set calculation boundaries for non-wing calcs
xlim2 = 6.8 - y;
xlim4 = 18.1*y + 1.65;
% avoid W4 algorithm for small y
if y<=1e-6
	xlim1 = xlim0;
	xlim2 = xlim0;
end
% look at region 0
region0 = find(abx > xlim0);

if length(region0)==length(x)			% for wing calcs use region 0 only
   g = (0.56418958*y)./(xq + yq);	
else											% other regions must be computed
	% set flagged array to indicate region calculated
	region = ones(size(x));
	% configure array for results
	g = zeros(size(x));
   % Region 0
   g(region0) = (0.56418958*y)./(xq(region0) + yq);
   region(region0) = 0;					% mark as comuputed
   
   % Region I - Humlicek W4 region I
   region1 = find(abx > xlim1 & region);
   if length(region1)~=0
%      disp('Region 1')
	   A0 = yq + 0.5;
      D2 = 2*yq - 1.0;
      xq1 = xq(region1);
	   g(region1) = (y*0.56418958) ./ (A0*A0 + xq1.*(D2+xq1)).*(A0 + xq(region1));
	   region(region1) = 0;					% mark as computed
	end   
   % Region II - Humlicek W4 region II
   region2 = find(abx > xlim2 & region);
   if length(region2)~=0
%      disp('Region 2')
	   xq2 = xq(region2);
	   H0 = 0.5625 + yq .* (4.5 + yq .*(10.5 + yq .*(6.0 + yq)));
	   H2 = -4.5 + yq .* (9.0 + yq.*(6.0 + yq .* 4.0));
	   H4 = 10.5 - yq .* (6.0 - yq.*6.0);
	   H6 = -6.0 + yq * 4.0;
	   E0 = 1.875 + yq.*(8.25 + yq.*(5.5 + yq));
      E2 = 5.25 + yq .* (1.0 + yq * 3.0);
      
      g(region2) = (0.56418958*y) ./ (H0 + xq2 .* (H2 + xq2 .*(H4 + xq2 .*(H6 + xq2)))) ...
         		.*(E0 + xq2.*(E2 + xq2.*(0.75 * H6 + xq2)));
	   region(region2) = 0;
	end	   
   % Region III - Humlicek W4 region III
   region3 = find(abx < xlim3 & region);
   if length(region3)~=0
%      disp('Region 3')
	   Z0 = 272.1014 + y.*(1280.829 + y.*(2802.870 + y.*(3764.944 ...
	      + y.*(3447.629 + y.*(2256.981 + y.*(1074.409 + y.*(369.1989 ...
	      + y.*(88.26741 + y.*(13.39880 + y)))))))));
	   Z2 = 211.678 + y.*(902.3066 + y.*(1758.336 + y.*(2037.310 ...
	      + y.*(1549.675 + y.*(793.4273 + y.*(266.2987 ...
	      + y.*(53.59518 + y.*5.0)))))));
	   Z4 = 78.86585 + y.*(308.1852 + y.*(497.3014 + y.*(479.2576 ...
	      + y.*(269.2916 + y.*(80.39278 + y*10.0)))));
	   Z6 = 22.03523 + y.*(55.02933 + y.*(92.75679 + y.*(53.59518 ...
	      + y*10)));
	   Z8 = 1.496460 + y.*(13.39880 + y*5.0);
	   
	   P0 = 153.5168 + y.*(549.3954 + y.*(919.4955 + y.*(946.8970 ...
	      + y.*(662.8097 + y.*(328.2151 + y.*(115.3772 + y.*(27.93941 ...
	      + y.*(4.264678 + y.*0.3183291))))))));
	   P2 = -34.16955 + y.*(-1.322256 + y.*(124.5975 + y.*(189.7730 ...
	      + y.*(139.4665 + y.*(56.81652 + y.*(12.79458 ...
	      + y*1.2733163))))));
	   P4 = 2.584042 + y.*(10.46332 + y.*(24.01655 + y.*(29.81482 ...
	      + y.*(12.79568 + y*1.9099744))));
	   P6 = -0.07272979 + y.*(0.9377051 + y.*(4.266322 + y*1.273316));
	   P8 = 0.0005480304 + y*0.3183291;
	   
	   xq3 = xq(region3);   
      g(region3) = 1.7724538 ./ (Z0 + xq3.*(Z2 + xq3.*(Z4 + xq3.*(Z6 + xq3.*(Z8 + xq3))))) ...
         				.*(P0 + xq3.*(P2 + xq3.*(P4 + xq3.*(P6 + xq3.*P8))));
	   region(region3)=0;
	end   
   
   % Region 4 - Humlicek CPF12 algorithm
   region4 = find(region);
   if length(region4)~=0
%      disp('Region 4')
      L = length(region4);
      YO=1.5;
      YOQ = YO*YO;
      YPYO = y + YO;
      YPYOQ = YPYO*YPYO;
	   % required data
	   C = [1.0117281, -0.75197147, 0.012557727, 0.010022008, -0.00024206814, 0.00000050084806];
	   S = [1.393237, 0.23115241, -0.15535147, 0.0062183662, 0.000091908299, -0.00000062752596];
	   T = [0.31424038, 0.94778839, 1.5976826, 2.2795071, 3.0206370, 3.8897249];
      
      D = zeros(6,L);
      for j=1:6
         D(j,:) = x(region4)-T(j);
      end
     	MQ = D.*D;
      MF = 1.0 ./(MQ + YPYOQ);
      XM = MF.*D;
      YM = MF.* YPYO;
      for j=1:6
         D(j,:) = x(region4)+T(j);
      end
     	PQ = D.*D;
      PF = 1.0 ./ (PQ + YPYOQ);
      XP = PF .* D;
      YP = PF .* YPYO;
   
      inda=find(abx(region4) <= xlim4);
      if ~isempty(inda)
         ga=zeros(1,length(inda));
	      for j = 1:6
	         ga = ga + C(j).*(YM(j,inda)+YP(j,inda)) - S(j)*(XM(j,inda)-XP(j,inda));
         end
	      g(region4(inda))=ga;
      end
      indb=find(abx(region4) > xlim4);
      if ~isempty(indb)
         gb=zeros(1,length(indb));
         YF = y + YO + YO;
	      for j = 1:6
            gb = gb + ...
               (C(j).*(MQ(j,indb).*MF(j,indb)-YO*YM(j,indb)) + S(j)*YF.*XM(j,indb))./(MQ(j,indb)+YOQ) + ...
               (C(j).*(PQ(j,indb).*PF(j,indb)-YO*YP(j,indb)) - S(j)*YF.*XP(j,indb))./(PQ(j,indb)+YOQ);
         end
         gb = y*gb + exp(-xq(region4(indb)));
	      g(region4(indb))=gb;
      end 
    end
end   