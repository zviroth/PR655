function [] = runPhotometer(numSamples, filename)
%things to check:
% JoshCalibforBL.m --> getPRValues - spec = PR655rawspd(10); so no timeout
% photometer cap off
dirname = '/Users/rothzn/Documents/MATLAB/calibration';

info.computer = 'laptop';
info.display = 'op4';
info.photometer = 'PR655';

info.date = date;
if strcmp(filename(end-3:end))==0
    filename = [filename '.mat'];
end

LumValues = BaseLum255(numSamples);

save([dirname filename],'LumValues','info');%in case createLUT fails

lut = createLUT(LumValues);

save([dirname filename],'LumValues','lut','info');




