function y = Filter0(b, x)
%Based on code originally written by Ken Harris

if size(x,1) == 1
	x = x(:);
end

if mod(length(b),2)~=1
	error('filter order should be odd');
end

shift = (length(b)-1)/2;

[y0 z] = filter(b,1,x);
y = [y0(shift+1:end,:,:,:) ; z(1:shift,:,:,:)];
