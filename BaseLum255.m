function [LumValues] = BaseLum255(numSamples)
global windowPTR gray PRport bitDepth
% Runs through R G B gun values and takes luminance values, returns array with values
% Use in first steps of Calibration 
tic
% whichScreen= 1;%max(Screen('Screens'));
whichScreen= max(Screen('Screens'));
bitDepth = 8;
PRport = 'COM4';
gray = [128 128 129];



% values = (linspace(0,1,15) .^(1/2.4) ) * 255 ;
values = (linspace(0,1,numSamples) .^(1/2.4) ) * 255 ;



Screen('Preference', 'SkipSyncTests', 1);

[windowPTR,screenRect] = Screen('Openwindow',whichScreen,32768,[],32,2);
LumValues = [];
for g = 1:3
reading = 1;
	for i = 1:numel(values)
		if g == 1
		[xyYcie xyYJudd Spectrum] = JoshCalibforBL([values(i) 0 0]);
		LumValues.red(reading,1).gunValue = i;
		LumValues.red(reading,1).xyYcie = xyYcie;
		LumValues.red(reading,1).xyYJudd = xyYJudd;
		LumValues.red(reading,1).Spectrum = Spectrum;
        %{
			if i > 252
			[xyYcie xyYJudd Spectrum] = JoshCalibforBL([values(i) 0 0]);
			LumValues.red(reading+1,1).gunValue = i;
			LumValues.red(reading+1,1).xyYcie = xyYcie;
			LumValues.red(reading+1,1).xyYJudd = xyYJudd;
			LumValues.red(reading+1,1).Spectrum = Spectrum;
			end
        %}
		elseif g == 2
		[xyYcie xyYJudd Spectrum] = JoshCalibforBL([0 values(i) 0]);
		LumValues.green(reading,1).gunValue = i;
		LumValues.green(reading,1).xyYcie = xyYcie;
		LumValues.green(reading,1).xyYJudd = xyYJudd;
		LumValues.green(reading,1).Spectrum = Spectrum;
        %{
			if i > 252
			[xyYcie xyYJudd Spectrum] = JoshCalibforBL([0 255 0]);
			LumValues.green(reading+1,1).gunValue = i;
			LumValues.green(reading+1,1).xyYcie = xyYcie;
			LumValues.green(reading+1,1).xyYJudd = xyYJudd;
			LumValues.green(reading+1,1).Spectrum = Spectrum;
			end
        %}
		elseif g == 3
		[xyYcie xyYJudd Spectrum] = JoshCalibforBL([0 0 values(i)]);
		LumValues.blue(reading,1).gunValue = i;
		LumValues.blue(reading,1).xyYcie = xyYcie;
		LumValues.blue(reading,1).xyYJudd = xyYJudd;
		LumValues.blue(reading,1).Spectrum = Spectrum;
        %{
		if i > 252
			[xyYcie xyYJudd Spectrum] = JoshCalibforBL([0 0 255]);
			LumValues.blue(reading+1,1).gunValue = i;
			LumValues.blue(reading+1,1).xyYcie = xyYcie;
			LumValues.blue(reading+1,1).xyYJudd = xyYJudd;
			LumValues.blue(reading+1,1).Spectrum = Spectrum;
        end
            %}
		end
		reading = reading + 1;
		disp(LumValues)
	end
end
Screen('CloseAll')
toc