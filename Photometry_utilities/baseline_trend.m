function [trend,detrended,p]=baseline_trend(x)
x=x(:);
p=fit((501:length(x)-400)',ordfilt2(x(501:end-400),20,true(1000,1),'symmetric'),'exp2');
trend=p(1:length(x));
detrended=x./trend*trend(1);