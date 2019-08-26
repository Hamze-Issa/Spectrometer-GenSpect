function g=voigthum(x,y)

% g=VOIGTHUM(x,y) generates a voigt profile using the complex 
%   error function approximation computed by region:
%
%  Region I/II: from J. Humlicek, "optimised Computation of the Voigt and 
%					 complex Probability Functions", JQSRT Vol.27, No.4 pp437-444,1982. 
%  Region III:  from A.K. Hui et al, "Rapid Computation of the Voigt and 
%  				 Complex Error Functions", JQSRT 1978.
%
%	This function uses the p=6 approximation and has error less than
%	10^(-6)
%
% VOIGTHUM calls no other genasis functions
%
%
% Modified for speed 23-09-98 by Ben Quine.

T = y - i*x;
S = abs(x)+y;

region1 = find(S >= 15);

if length(region1)==length(S)			% for wing calcs use region 1 only
   g = real(T.*0.5641896./(0.5+T.*T));	
else
	g = zeros(size(S));
	% Region I
	T1 = T(region1);
	g(region1) = real(T1.*0.5641896./(0.5+T1.*T1));	
   
   region2 = find(S >= 5.5 & S < 15);
	% Region II
	T1 = T(region2);
	U1 = T1.*T1;
	g(region2) = real(T1.*(1.410474+U1*0.5641896)./(0.75+U1.*(3+U1)));
   
   if (length(region1)+length(region2))<length(S)			% use region 1 and 2 only?
		region3 = find(S < 5.5);
	
		% region III
		%define the coefficients for the rational approximation 
		a= [122.607931777104326 214.382388695706425 181.92853309218549 93.155580458138441 30.180142196210589 5.9126209773153 0.564189583562615];
		b= [122.607931773875350 352.730625110963558 457.334478783897737 348.703917719495792 170.354001821091472 53.992906912940207 10.479857114260399];
		T1=T(region3);
		g(region3)=real(((((((a(7)*T1 + a(6) ).*T1 + a(5)).*T1 + a(4)).*T1 + a(3)).*...
	   T1 + a(2)).*T1+a(1))./(((((((T1+b(7)).*T1+b(6)).*T1+b(5)).*...
	   T1+b(4)).*T1+b(3)).*T1+b(2)).*T1+b(1)));
   end
end