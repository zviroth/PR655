function lut = createLUT(LumValues)

%create look up table
for i = 1:size(LumValues.red,1)
	gunVal(i) = LumValues.red(i).gunValue;
	R_Y(1,i) = LumValues.red(i).xyYJudd(3); 
	G_Y(1,i) = LumValues.green(i).xyYJudd(3); 
	B_Y(1,i) = LumValues.blue(i).xyYJudd(3); 
end

lut.x = [0:65535];
lut.R = makeLUTperGun(gunVal, R_Y);
lut.G = makeLUTperGun(gunVal, G_Y);
lut.B = makeLUTperGun(gunVal, B_Y);

end 


function xlut = makeLUTperGun(xdata, ydata)
%MAKELUT - Create lookup table for gamma correction.
%
%MAKELUT(XDATA, YDATA) creates a lookup table for
%gamma correction based on the pixel values XDATA and the measured
%luminances YDATA for these pixel values. XDATA and YDATA must be of the
%same size. 
%
% Example
% xdata = [0 4 8 16 24 32 64  96 128 140 150 160 170 180 190 200 210 220 ...
%           230 240 250 255];
% ydata = [0.15 0.15 0.16 0.18 0.19 0.22 0.38 0.67 1.11 1.24 1.39 1.71 ...
%          1.90 2.09 2.29 2.75 3.00 3.25 3.52 4.11 4.42 4.43];
% output = makelut(xdata, ydata);
%

% (c) 2005-10-19  Thorsten Hansen

ydata = (ydata-min(ydata))/(max(ydata)-min(ydata));
xdata = xdata /65535; % check the initial values of xdata RB


% initial estimates of values to be fitted
alpha_start = min(ydata); % Is 0 with normalization RB
beta_start =  max(ydata)/max(xdata).^gamma_start; % Is 1 anyway with normalization RB


% now do the fitting
p0 = [alpha_start beta_start gamma_start]; % define start values 
                                          % of the parameters to be fitted
fun = inline('p(1) + p(2)*xdata.^p(3)', 'p', 'xdata');
pfit = lsqcurvefit(fun, p0, xdata, ydata, min(ydata), max(ydata));

% create lut
% by evaluating the inverse function of the gamma fit 
% y = alpha+beta*x.^gamma => x = ((y-alpha)/beta).^(1/gamma);

alpha = pfit(1); beta = pfit(2); gamma = pfit(3);
N = 65536; % size of lut
% create lut based on measured data (probably Karl's lino approach)
ylut = linspace(min(ydata), max(ydata), N);

% creata lut based on yvalues of FITTED!!! function
% pro:  real(xlut) below is not necessary
% pro:  smooth function even if gamma_start is not that good
% con:  but some deviation from lino fit % not really, bad fit can come from noisy measurements
%
% ylut = linspace(feval(fun, pfit, xdata(1)), feval(fun, pfit, xdata(end)), N);

% inverse function 
xlut = (ylut).^(1/gamma);
xlut = real(xlut); % if alpha > min(ylut), the gamma-sqrt of a neg. value
                   % is taken, resulting in complex values

xlut = round(xlut*65535);

end 
