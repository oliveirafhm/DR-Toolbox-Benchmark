% Function name....: ROCAdjust
% Date.............: July 07, 2016
% Author...........: Fabio Henrique, (oliveirafhm@gmail.com)
% Description......:
%                    Full ROC vector until reach nSamples using max value
% Parameters.......:
%                    xIn -> False positive rate (array)
%                    yIn -> True positive rate (array)
%                    nSamples -> Out number of samples
%                    toPlot -> 0 (false); 1 (true)
% Return...........:
%                    xOut -> Uniformly spaced x values
%                    yOut -> Interpolated values
% Remarks..........:
%
function [xOut,yOut] = ROCAdjust(xIn,yIn,nSamples,toPlot)

xMax = max(xIn);
yMax = max(yIn);
lIn = length(xIn);
xOut = zeros(nSamples,1);
yOut = zeros(nSamples,1);
if lIn < nSamples
    for i=1:lIn
        xOut(i) = xIn(i);
        yOut(i) = yIn(i);
    end
    for i=lIn+1:nSamples
        xOut(i) = xMax;
        yOut(i) = yMax;
    end
elseif lIn > nSamples
    error('ROC x greater than maximum n samples.');
else
    xOut = xIn;
    yOut = yIn;
end

if toPlot == 1
    figure;
    plot(xIn,yIn,'o',xOut,yOut);
    title('ROC fixed lenght');
    xlabel('False positive rate');
    ylabel('True positive rate');
end
end

