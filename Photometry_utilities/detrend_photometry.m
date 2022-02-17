function [signal_corrected,detrended_signal,detrended_iso,dFF]=detrend_photometry(signal,iso)
%[signal_corrected,detrended_signal,detrended_iso,dFF]=detrend_photometry(signal,iso)
%
%Detrend calcium-dependent trace by an isosbestic trace
%
%INPUT:
%signal - calcium-dependent trace
%iso - isosbestic trace
%
%OUTPUT:
%signal_corrected - signal detrended and with artifacts removed
%detrended_signal - signal just detrended
%detrended_iso - detrended isosbestic trace
%dFF - deltaF/F relative to a running baseline [not used in paper]
%
%Joseph M. Stujenske, Feb 2022
%
limit_length=min(length(signal),length(iso));
[trend,detrended_iso,p]=baseline_trend(iso);
detrended_iso2=detrended_iso(1:limit_length);
trend_sig=p(1:length(signal));
signal_base2=movmax(movmin(medfilt1(signal(11:limit_length-300),100,'truncate'),2000),2000);
poly_coefs_sig=polyfit(trend_sig(11:limit_length-300),signal_base2,1);
trend_sig=polyval(poly_coefs_sig',trend_sig);
poly_coefs=polyfit(trend(1:limit_length),trend_sig(1:limit_length),1);
detrended_signal=signal./trend_sig*trend_sig(1);
detrended_signal2=detrended_signal(1:limit_length);
iso_scale=polyval(poly_coefs',detrended_iso2);
signal_corrected=detrended_signal2-iso_scale+nanmean(iso_scale);
baseline_running=ordfilt2(signal_corrected,5,true(500,1),'symmetric');
dFF=((signal_corrected)-baseline_running)./nanmean(baseline_running);