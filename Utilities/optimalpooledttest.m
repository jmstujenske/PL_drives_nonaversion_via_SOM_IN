function [Tpool,dfpool,pval]=optimalpooledttest(X,Y)
%[Tpool,dfpool,pval]=optimalpooledttest(X,Y)
%Optimal pooled t-test T0pool from Guo B & Yuan Y, Statistical Methods in
%Medical Research 2007, 201 26(3): p 1328.
%
%Inputs:
%X,Y = column or row vectors of the same length with one set of partially
%paired and unpaired observations, where missing data points are marked
%with NaN
%
%Outputs:
%Tpool=the optimal T0pool statistic
%dfpool=the effective degrees of freedom for the t statistic
%p-value for the null hypothesis that X and Y are drawn from the same
%distribution
%
%Joseph M. Stujenske, 2020
%
toremove=isnan(X) & isnan(Y);
X(toremove)=[];
Y(toremove)=[];
pairedvars=(~isnan(X) & ~isnan(Y));
unpairedvars=(isnan(X) | isnan(Y));
n=sum(pairedvars);
if n>0 && sum(unpairedvars)>0
n0=sum(isnan(Y));
n1=sum(isnan(X));
u1=nanmean(X(pairedvars)-Y(pairedvars));
SD=nanvar(X(pairedvars)-Y(pairedvars));
SDx=nanvar(X(isnan(Y)));
SDy=nanvar(Y(isnan(X)));
S1=SD/n;
S2=SDx/n0+SDy/n1;
u2=nanmean(X(unpairedvars))-nanmean(Y(unpairedvars));
wtot=1./S1+1./S2;
w1=1./S1./wtot;
w2=1./S2./wtot;
df1=sum(pairedvars)-1;
df2=((S1).^2)./((SDx/n0).^2./(n0-1)+(SDy/n1).^2./(n1-1));
Su=(1+4*w1.*(1-w1)/df1+4*w2.*(1-w2)/df2)/wtot;
u=w1*u1+w2*u2;
Tpool=u./sqrt(Su);
dfpool=1./(w1.^2/df1+w2.^2/df2);
P=tcdf(Tpool,dfpool);
Pinv=1-P;
pval=min(2*P,2*Pinv);
else
    error('Data is either all paired or all unpaired; Use a regular t-test.')
end