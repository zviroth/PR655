function [xyYJudd] = s2jJ(SPEC)

 
%spectrum = load('test_spec.txt','-ascii'); %% This is the way for importing the data in the
                                           %% test_spec file preserving the original numbers.
%spectrum=textread('test_spec.txt', '', 'delimiter', ','); %read in the text file with the spectrum in it
% plot(spectrum(:,1), spectrum(:,2), '-ko');
% set(1, 'color', [ 1 1 1]);


SPEC(:,3)=SPEC(:,2); % a clooge because spectra2xyY requires 2 spectra minimum

%JUDDxyY_2=spectra2xyY(SPEC, 'judd');%CHANGE moncie!!! (but where?)
%xyY_2=spectra2xyY(SPEC, 'cie1931');
xyY_1=spectra2xyY(SPEC, 'judd');

%xyYCIE = xyY_2(1,:);
xyYJudd = xyY_1(1,:);
end










