function [bsuccess] = CPforBLnoPTB(mx)
% ColourPatch.m
global windowPTR gray


% Open the chosen screen.
% For an explanation of the parameters, type Screen OpenwindowPTR? at the MatLab prompt.

%Screen('LoadNormalizedGammaTable',windowPTR,linspace(0,(255/256),256)'*ones(1,3));
%Screen('LoadNormalizedGammaTable',windowPTR,linspace(0,1,256)'*ones(1,3));

%Clut( 1, : ) = gray(1,:);
%BitsPlusSetClut(windowPTR,Clut);
% Draw the square on both screens.
%[screenWidth, screenHeight]=Screen('windowSize', windowPTR);
pos = get(gcf, 'Position'); %// gives x left, y bottom, width, height
screenWidth = pos(3);
screenHeight = pos(4);
size = min([screenWidth screenHeight]);

squaresize = size / 6;
xpos = screenWidth/2; 
ypos = screenHeight/2;
%Screen('FillRect', windowPTR, gray );
%Screen('FillOval', windowPTR, mx, [ xpos-squaresize/2 ypos-squaresize/2 xpos+squaresize/2 ypos+squaresize/2 ] );

%plot(x,y)
set(gca,'Color',mx/255)
%disp(['Reading #' num2str(i) ': ' num2str(mx(i)*[1 0 0])]);
disp(['Reading #' num2str(1) ': ' num2str(mx(1,:))]);
%Screen( windowPTR, 'Flip', 0, 1 );


bsuccess = 1;

