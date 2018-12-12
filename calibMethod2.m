%	- 9/25/11 -
%
%	This script computes from the three R, G and B guns (characterized in
%   the CIE31 color space) the L, M, and S corresponding responses, and
%   the transformation matrix allowing color space conversion.
%	The data can be used to determine the cone contrast elicited by a 
%   stimulus, etc.  The values from this program are also used to construct
%   the ldrgyv2rgb matrix. This program should be run AFTER a monitor 
%   calibration, and so with the most recent XYZ values.
%   The solver is based on Zaidi and Halevy (1993).
%   Romain Bachy

%   04/28/13
%
%   Equations have be rewritten to make them readable for human researchers.
%   Comments and questions have been added for human researchers.
%   Romain Bachy

%  01/17/14
%
%  Comment improvements
%  Computed value easier to read ( and so to be double checked)
%  Romain Bachy

%% Main Function
function [ldrgyvMAT, lmrgyvMAT] = calibMethod2
%clear al; close all; clc;

[file,path] = uigetfile; 
NAME  = who('-file', [path, file]); 
DAT   = load([path, file]); 
RED   = DAT.(NAME{1}).red; 
GREEN = DAT.(NAME{1}).green; 
BLUE  = DAT.(NAME{1}).blue; 

rx = RED(end).xyYcie(1); 
gx = GREEN(end).xyYcie(1); 
bx = BLUE(end).xyYcie(1); 
ry = RED(end).xyYcie(2); 
gy = GREEN(end).xyYcie(2);
by = BLUE(end).xyYcie(2); 
white(1) = RED(end).xyYcie(3);
white(2) = GREEN(end).xyYcie(3);
white(3) = BLUE(end).xyYcie(3);

%We enter x,y relative values of R, G and B guns (experimental measurements)

%        x     y     z
%   R = 0.632, 0.343, 0.0249
%   G = 0.292, 0.611, 0.0971
%   B = 0.148, 0.0749, 0.777 


%R, G and B into L, M, and S answers in a cone contrast space (/L+M so
%is there a von kries adaptation here?)
  
 disp([char(10), 'R gun LMS values']);
  R    = cie2lms(rx, ry);
 disp([char(10), 'G gun LMS values']); 
  G    = cie2lms(gx, gy);
 disp([char(10), 'B gun LMS values']);
  B    = cie2lms(bx, by);

  % Let it makes easier to read for everybody:

  L(1) = R(1);
  L(2) = G(1);
  L(3) = B(1);
  M(1) = R(2);
  M(2) = G(2);
  M(3) = B(2);
  S(1) = R(3);
  S(2) = G(3);
  S(3) = B(3);
  
  % index 1 for LMS values of R
  % index 2 for LMS values of G
  % index 3 for LMS values of B  
  
  Lr = R(1);
  Lg = G(1);
  Lb = B(1);
  Mr = R(2);
  Mg = G(2);
  Mb = B(2);
  Sr = R(3);
  Sg = G(3);
  Sb = B(3);  

    % Numerators are luminances of Rmax, Gmax, and Bmax  [cd/m^2] 
  % (experimental measurements or computed from absolute spectral distribution)

    
  Wr = white(1);
  Wg = white(2);
  Wb = white(3);

  disp([char(10), 'Middle Grey Point: LMS responses']);
  ConeContrastMG = lumchrm (white(1), white(2), white(3), L, M, S);

  disp([char(10), 'R/G Axis: G and B modulation']); 
  MyNewdBrg2 = newsolvex((Lg+Mg), Sr, (Lr+Mr), Sg, Sb, (Lb+Mb));
  MyNewdGrg2 = newsolvey((Lg+Mg), Sr, (Lr+Mr), Sg, Sb, (Lb+Mb));
  
  %The solver computes the relative variation of 2 guns with regard to
  %the third one. We fixed the R gun value at first, so we have 
  %to set up everything relatively to the respective luminance values 
  %(from min to max luminance value) for the matrix construction:
  
  %Here, I believe that luminance min has been put equal to 0 (is it good?)
  
  dGrg = diff(MyNewdGrg2)*Wr/Wg;
  dBrg = diff(MyNewdBrg2)*Wr/Wb;
  
  disp(['dGrg = ',num2str(dGrg),' and dBrg =', num2str(dBrg), char(10)]);

  disp([char(10), 'Gamut extrema on R/G axis (when R varies)']); 
  disp('RG axis: LMS responses for R = 0:');
  ConeContrastML = lumchrm(0.0, Wg-diff(MyNewdGrg2)*Wr, Wb-diff(MyNewdBrg2)*Wr, L, M, S);

  disp([char(10), 'RG axis: LMS responses for R = 1:']);
  ConeContrastLM = lumchrm(2*Wr, Wg+diff(MyNewdGrg2)*Wr, Wb+diff(MyNewdBrg2)*Wr, L, M, S);

  %lumchrm(0.0, Wg-diff(MyNewdGrg2)*Wr, Wb-diff(MyNewdBrg2)*Wr, L, M, S)
  %lumchrm(Wr*2.0, Wb+diff(MyNewdGrg2)*Wr, Wg+diff(MyNewdBrg2)*Wr, L, M, S);

  disp([char(10), 'B/Y Axis: R and G modulation']);
  
  MyNewdGyv2 = newsolvex(Mr, Lb, Mb,Lr, Lg, Mg);
  MyNewdRyv2 = newsolvey(Mr, Lb, Mb, Lr, Lg, Mg);
  
  dRyv    = diff(MyNewdRyv2)*Wb/Wr;
  dGyv    = diff(MyNewdGyv2)*Wb/Wg;

  disp(['dRyv = ',num2str(dRyv),' and dGyv =', num2str(dGyv), char(10)]);
  
  disp([char(10), 'Gamut extrema on B/Y axis (when B varies)']); 

  disp('YV axis: LMS responses for B = 0:');
  ConeContrastY = lumchrm(Wr-diff(MyNewdRyv2)*Wb, Wg-diff(MyNewdGyv2)*Wb, 0.0, L, M, S);
  disp([char(10), 'YV axis: LMS responses for B = 1:']); 
  ConeContrastB = lumchrm(Wr+diff(MyNewdRyv2)*Wb, Wg+diff(MyNewdGyv2)*Wb, Wb*2.0, L, M, S);

%     disp([char(10), 'Cone contrast responses for each pole']); 
%   % Display the L, M and S cone contrast responses on each axis at contrast maximum:
% ConeContrastMG   % response for mid grey level
% ConeContrastLD = 2* ConeContrastMG % response to white
% ConeContrastLM   % response to L-M
% ConeContrastML   % response to M-L
% ConeContrastY    % response to Y
% ConeContrastB    % response to V


    disp([char(10), 'Cone contrast excursion for each 1/2 axis']); 
  % Compute and display the L, M and S cone contrast excursions on each axis btw contrast
  % maximum and the mid grey level:
DeltaBMax   = ConeContrastB - ConeContrastMG 
% gives the displayable maximum cone contrast on YV axis
%DeltaYMax   = ConeContrastY - ConeContrastMG % only here to check symmetry
DeltaLMMax  = ConeContrastLM - ConeContrastMG 
% gives the displayable maximum cone contrast on RG axis
DeltaMLMax  = ConeContrastML - ConeContrastMG % only here to check symmetry

%DeltaLDMax = ConeContrastWhite - ConeContrastMG is equal to ConeContrastW
% by definition
DeltaLDMax  = ConeContrastMG
% gives the displayable maximum cone contrast on LD axis axis

  disp([char(10), 'S cone excursion on the LM axis: Displayable ratio']); 
% compute the ratio of displayable values when luminance axis is LM (so when there is S cone excursion)

RatioDispLM = DeltaBMax(3)/DeltaLDMax(3)

disp([char(10), 'L and M cone excursion on the LM axis']); 
% compute the ratio of displayable values when luminance axis is LM in
% order to cancel the L or M cone excursion from L-M axis

RatioDispISOM  = DeltaLMMax(1)/(RatioDispLM*DeltaLDMax(1)) 
% what ratio of L cone excursion on the LM luminance axis? to be cancelled for
% M cone isolating stimulus
RatioDispISOL = DeltaMLMax(1)/(RatioDispLM*DeltaLDMax(2)) 
% what ratio of M cone excursion on the LM luminance axis? to be cancelled for
% L cone isolating stimulus

%% Computes the matrices
% DKL with LMS
ldrgyvMAT=[ 1    1	     dRyv                      
            1    dGrg    dGyv               
   		    1    dBrg    1     ] 
        
% Matrix computed for:
% LMS as first column (R = G = B)
% RG as second column  (dR is max, dG and dB are calculated) 
% YV as third column   (dB is max, dR and dG are calculated)      

% DKL with LM
lmrgyvMAT=[ RatioDispLM - dRyv   1	     dRyv                      
            RatioDispLM - dGyv   dGrg    dGyv               
   		    RatioDispLM - 1      dBrg      1     ]      
    
% Matrix computed for:
% LM as first column (ldrgyv first column ld - ldrgyv third column yv)
% The ratio is applied to the ld column then we do Ratio*ld - yv
% RG as second column   
% YV as third column  

% Cone isolating matrix
lmsMAT=[    -RatioDispISOL*lmrgyvMAT(1,1) + 1         RatioDispISOM*lmrgyvMAT(1,1) - 1 	       dRyv                      
            -RatioDispISOL*lmrgyvMAT(2,1) + dGrg      RatioDispISOM*lmrgyvMAT(2,1) - dGrg      dGyv               
   		    -RatioDispISOL*lmrgyvMAT(3,1) + dBrg      RatioDispISOM*lmrgyvMAT(3,1) - dBrg       1     ]  


          disp('Before normalization');
          disp('M axis');
  ConeContrastMp = lumchrm(Wr+lmsMAT(1,2)*Wr,...
                           Wg+lmsMAT(2,2)*Wg, ...
                           Wb+lmsMAT(3,2)*Wb, L, M, S);
  disp([char(10), 'M axis']); 
  ConeContrastMm = lumchrm(Wr-lmsMAT(1,2)*Wr,...
                           Wg-lmsMAT(2,2)*Wg, ...
                           Wb-lmsMAT(3,2)*Wb, L, M, S);
                       
           disp('L axis');
  ConeContrastLm = lumchrm(Wr-lmsMAT(1,1)*Wr,...
                           Wg-lmsMAT(2,1)*Wg, ...
                           Wb-lmsMAT(3,1)*Wb, L, M, S);
  disp([char(10), 'L axis']); 
  ConeContrastLp = lumchrm(Wr+lmsMAT(1,1)*Wr,...
                           Wg+lmsMAT(2,1)*Wg, ...
                           Wb+lmsMAT(3,1)*Wb, L, M, S); 
                       
                       
% Cone isolating matrix with normalization
lmsMATnorm = [lmsMAT(:,1)/max(abs(lmsMAT(:,1)))  lmsMAT(:,2)/max(abs(lmsMAT(:,2)))   lmsMAT(:,3)]

          disp('After normalization');
          disp('M axis');
  ConeContrastMp = lumchrm(Wr+lmsMATnorm(1,2)*Wr,...% should be 2*Wr
                           Wg+lmsMATnorm(2,2)*Wg, ...
                           Wb+lmsMATnorm(3,2)*Wb, L, M, S);
  disp([char(10), 'M axis']); 
  ConeContrastMm = lumchrm(Wr-lmsMATnorm(1,2)*Wr,... % should be 0
                           Wg-lmsMATnorm(2,2)*Wg, ...
                           Wb-lmsMATnorm(3,2)*Wb, L, M, S);
                       
           disp('L axis');
  ConeContrastLm = lumchrm(Wr-lmsMATnorm(1,1)*Wr,...  % should be 0
                           Wg-lmsMATnorm(2,1)*Wg, ...
                           Wb-lmsMATnorm(3,1)*Wb, L, M, S);
  disp([char(10), 'L axis']); 
  ConeContrastLp = lumchrm(Wr+lmsMATnorm(1,1)*Wr,... % should be 2*Wr
                           Wg+lmsMATnorm(2,1)*Wg, ...
                           Wb+lmsMATnorm(3,1)*Wb, L, M, S);                      
% I end here with same matrix as the new solver 


saveFlag = questdlg('Save DKL Transformation Matrix?', 'save', 'yes', 'no', 'yes');
if saveFlag == 'yes'
	[~, fname, ext] = fileparts(file);
	save([path, fname, '_DKLtransformMat', ext], 'ldrgyvMAT'); 
end 

return
 

%% Internal Functions

function [norm] = cie2lms(a, b)

  x   = a;
  y   = b;
  z   = 1-x-y;
  cie = [x y z]';
  %I have relative z values, it may be better to use X Y and Z
  %values, but since I introduce Luminance (Y) in a later stage, does it really
  %matter? Maybe still better for numerical consideration. Would it be even
  %better with spectral distribution?

  matrix=[ .15514 .54316 -.03286;
		  -.15514 .45684  .03286;
		  0		  0		  .01608];
  
      
  % Is it a S&P matrix or a S&S one?

  
  lms = matrix*cie;
  lms = lms';
  disp(['LMS = ',num2str(lms)]);
  l   = lms(1)/(lms(1)+lms(2)); %   L/(L+M)
  m   = lms(2)/(lms(1)+lms(2)); %   M/(L+M)
  s   = lms(3)/(lms(1)+lms(2)); %   S/(L+M)

  norm=[l m s];
  disp(['lms = ',num2str(norm)]);

  % after normalisation I have l+m = 1

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RB function: gives the two extrema independantly of the luminances. Just
% have to make a subtraction and taking in account luminances in a later
% stage. From Zaidi 1993

%example for the H axis and the relative B gun value variation with regard 
%to R fixed from 0 to 1:
%A = [(Lg+Mg)*Sr - (Lr+Mr)*Sg]/[(Lg+Mg)*Sb - (Lb+Mb)*Sg] and :
%Bh(Rh=0) = 1/2*(1+A) and
%Bh(Rh=1) = 1/2*(1-A)
%Delta Bh is then equal to the difference of both extrema and is actually
%equal to A !!!

function x = newsolvex(a, b, c, d, e, f)
  A  = (a*b - c*d)/(a*e - f*d);
  x1 = 1/2*(1 - A);
  x0 = 1/2*(1 + A);
  x  = [x0, x1];
return

function y = newsolvey(a, b, c, d, e, f)
  B  = (b*f - c*e)/(f*d - a*e);
  y1 = 1/2*(1 - B);
  y0 = 1/2*(1 + B);
  y  = [y0, y1]; 
  
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RB function 2:
% More simple function when I know that L + M = 1 in a cone contrast space
% L + M is constant for all computation here.

% function x = newsolvex2(a, b, c)
%   x1 = (1/2*c - 1/2*a)/(c - b);
%   x0 = (1/2*c + 1/2*a - b)/(c - b);
%   x = [x0, x1];
% return
% 
% function y = newsolvey2(a, b, c)
%   y1 = (1/2*b - 1/2*a)/(b - c);
%   y0 = (1/2*b + 1/2*a - c)/(b - c);
%   y = [y0, y1];  
%   
% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ConeContrastOut = lumchrm(lumr, lumg, lumb, L, M, S)

% This function computes LMS responses and luminance values.
% 1) Shows Y values of R, G, and B guns.
% 2) Compute L, M, and S absolute values of the considered color (Y values of for 
% R, G, B * L, M, or S value)
% 3) Shows L/(L+M) and S/(L+M) cone contrast space values

  TabLetter = ['R'; 'G'; 'B'];  

  lum(1) = lumr;
  lum(2) = lumg;
  lum(3) = lumb;
  bigl   = 0.0;
  bigm   = 0.0;
  bigs   = 0.0;

  for i = 1:3
	 %fprintf('Luminances[%s] = %10.6f \n', TabLetter(i), lum(i));
     disp(['Y ', TabLetter(i,:),' = ', num2str(lum(i)), ' cd/m^2']);
	 lu   = lum(i);
	 bigl = bigl + lu*L(i); % if i = 1 then L, M, S values of R.
	 bigm = bigm + lu*M(i);
	 bigs = bigs + lu*S(i);
  end

  %denom = 1;
  denom = bigl+bigm;
  lout = bigl/denom;
  mout = bigm/denom;  
  sout = bigs/denom;
  
  fprintf('L =%10.6f, M =%10.6f, S =%10.6f\n', bigl, bigm, bigs);
  fprintf('L/(L+M) =%10.6f, M/(L+M) =%10.6f\n , S/(L+M) =%10.6f\n', lout, mout, sout);
  disp(['Ytot = ', num2str(sum(lum))]);

  ConeContrastOut = [lout, mout, sout]; 
  
return
