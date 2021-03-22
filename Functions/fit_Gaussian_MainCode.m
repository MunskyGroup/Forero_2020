
function [fitresult,y,ysem] = fit_Gaussian_MainCode (y_input, ysem_input, rangeToDisplay)

x=-rangeToDisplay:rangeToDisplay;

totalRange = length(y_input);
minVal = round (totalRange/2 - rangeToDisplay);
maxVal = round (totalRange/2 + rangeToDisplay);

%Data using Unph as reference for the minimum
y=y_input(minVal:maxVal);
ysem=ysem_input(minVal:maxVal);

[x, y, ysem] = prepareCurveData( x, y, ysem);
fitresult = fitGaussian (x, y, ysem);

% y_display = y_input (minVal:maxVal);
% ysem_display = ysem_input (minVal:maxVal);

end