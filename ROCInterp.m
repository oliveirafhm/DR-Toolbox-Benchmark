% Function name....: ROCInterp
% Date.............: July 04, 2016
% Author...........: Fabio Henrique, (oliveirafhm@gmail.com)
% Description......:
%                    Interpolates ROC data using interp1 function
% Parameters.......:
%                    xIn -> False positive rate (array)
%                    yIn -> True positive rate (array)
%                    nSamples -> Out number of samples
%                    xLim -> Array of lower and upper x limits
%                    varargin:
%                      method -> Interpolation method (look intep1 matlab
%                      doc)(string)
%                      toPlot -> 0 (false); 1 (true)
% Return...........:
%                    xOut -> Uniformly spaced x values
%                    yOut -> Interpolated values
%                    uniqueXIn -> Unique values of xIn vector
%                    idxUniqueXIn -> Indexes of uniqueXIn
% Remarks..........:
%                    
function [xOut,yOut,uniqueXIn,idxUniqueXIn] = ROCInterp(xIn,yIn,nSamples,...
    xLim,varargin)

method = varargin{1};
toPlot = varargin{2};
if isempty(method)
    method = 'linear';
end
[uniqueXIn,idxUniqueXIn] = unique(xIn);
xOut = xLim(1):1/(nSamples-1):xLim(2);
yOut = interp1(uniqueXIn,yIn(idxUniqueXIn),xOut,method);
if toPlot == 1
  figure;
  plot(xIn,yIn,'o',xOut,yOut);
  title('ROC Interpolation');
  xlabel('False positive rate');
  ylabel('True positive rate');
end
end

