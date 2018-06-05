%% PolyDetrend - Removal of underlying trends by configurable polynomial
% 
% Function fits an n order polynomial, then subtracts said polynomial from
% input data.
%
% Created by:  D.C. Hartlen, EIT
% Date:        09-May-2018
% Modified by:  
% Date:        

function yyDetrend = PolyDetrend(xx,yy,n)

polyCoefs = polyfit(xx,yy,n);

yyDetrend = yy-polyval(polyCoefs,xx);


end