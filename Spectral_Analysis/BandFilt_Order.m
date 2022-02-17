function [Filtered,FiltHilb,FiltAmp,FiltPh] = BandFilt_Order(eeg,sampFreq,order,low,high)
%Based on code originally written by Ken Harris
%Modified by Joseph M. Stujenske in 2014 to allow for specification of filter order

Nyquist = sampFreq/2;
MyFilt=fir1(order,[low high]/Nyquist);

Filtered = Filter0(MyFilt,eeg);

if (nargout>1)
    FiltHilb = hilbert(Filtered);
    FiltPh = angle(FiltHilb);
    FiltAmp = abs(FiltHilb);
end
