function makelut
	
[file,path,indx] = uigetfile; 
NAME  = who('-file', [path, file]); 
DAT   = load([path, file]); 

for i = 1:size(DAT.(NAME{1}).red,1)
	gunVal(i) = DAT.(NAME{1}).red(i).gunValue;
	R_Y(1,i) = DAT.(NAME{1}).red(i).xyYJudd(3); 
	G_Y(1,i) = DAT.(NAME{1}).green(i).xyYJudd(3); 
	B_Y(1,i) = DAT.(NAME{1}).blue(i).xyYJudd(3); 
end

[file,path] = uiputfile; 
makeLUTperGun(gunVal, R_Y, [path, file, '_Rlum_LUT_', date], 2.2);
makeLUTperGun(gunVal, G_Y, [path, file, '_Glum_LUT_', date], 2.2);
makeLUTperGun(gunVal, B_Y, [path, file, '_Blum_LUT_', date], 2.2);

end 


function makeLUTperGun(xdata, ydata, outputfile, gamma_start)
%MAKELUT - Create lookup table for gamma correction.
%
%MAKELUT(XDATA, YDATA, OUTPUTFILE, [GAMMA_START]) creates a lookup table for
%gamma correction based on the pixel values XDATA and the measured
%luminances YDATA for these pixel values. XDATA and YDATA must be of the
%same size. The resulting lookup table is stored in OUTPUTFILE.
%Optinal parameter GAMMA_START is the starting value of gamma for the
%fitting routine, default is 2.5. 
%
% Example
% xdata = [0 4 8 16 24 32 64  96 128 140 150 160 170 180 190 200 210 220 ...
%           230 240 250 255];
% ydata = [0.15 0.15 0.16 0.18 0.19 0.22 0.38 0.67 1.11 1.24 1.39 1.71 ...
%          1.90 2.09 2.29 2.75 3.00 3.25 3.52 4.11 4.42 4.43];
% makelut(xdata, ydata, 'out.lut');
%
% % To Define estimate of gamma, which also generates some plots
% 
% makelut(xdata, ydata, 'out.lut', 3.0);
%
% (c) 2005-10-19  Thorsten Hansen
  
%cd C:\Research\MGH_MRI\stimuli\MGH_lowercontrast\
ydataold = ydata;
xdataold = xdata;
ydata = (ydata-min(ydata))/(max(ydata)-min(ydata));
xdata = xdata /65535; % check the initial values of xdata RB

if nargin < 4
  gamma_start = 2.2; % 2.5 was default; this should match the data fairly good, otherwise the fit
                     % may be bad
end

% initial estimates of values to be fitted
alpha_start = min(ydata); % Is 0 with normalization RB
beta_start =  max(ydata)/max(xdata).^gamma_start; % Is 1 anyway with normalization RB


% now do the fitting
p0 = [alpha_start beta_start gamma_start]; % define start values 
                                          % of the parameters to be fitted
fun = inline('p(1) + p(2)*xdata.^p(3)', 'p', 'xdata');
pfit = lsqcurvefit(fun, p0, xdata, ydata, min(ydata), max(ydata));

% plot fit for testing
% 
%plot(xdata, ydata)
%hold on
%h = plot(xdata, feval(fun, pfit, xdata), 'r-');

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
% save LUT x xlut to outputfile
x = [0:65535];
%x = [0:255];
fid = fopen(outputfile, 'w');
for i=1:N
  fprintf(fid, '%d %d\n', x(i), round(xlut(i)));
end
fclose(fid);




%
% some plot & stats if gamma value is given
%
if nargin > 3 
  
  subplot(2,1,1)
  % plot data
  
  
  plot(xdata*65535, (ydata*(max(ydataold)-min(ydataold))+min(ydataold)), '-o');

  % and function based on starting values
  hold on
  plot(xdata*65535, (xdata*65535).^(gamma_start), '--b')   

   
  % plot fitted function in red
  fun = inline('p(1) + p(2)*(xdata).^p(3)', 'p', 'xdata');
  h = plot(xdata*65535, feval(fun, pfit, xdata), 'r-');


  subplot(2,1,2)
  % plot the lut
  plot(x, xlut)
  axis([0 N-1 0 N-1])
  title('LUT')

 % compare to Karl's lino data
 % xlino = textread(['/home/hansen/bib/data/lut/sonygdmf520.' phosphor]);
 % hold on
 % plot(x, xlino, 'k')

 
 % some statistics
 yfit = feval(fun, pfit, xdata);
 
 residuals = yfit - ydata
 
 % root mean squared error
 %   RMS = sqrt(mean(residuals.^2))
 % same as standard deviation
 std(residuals)
 
end

end 
