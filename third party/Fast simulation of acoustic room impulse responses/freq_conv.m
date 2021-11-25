function rr = freq_conv(xx,yy)
%freq_conv  Performs convolution in frequency domain
%
% C = freq_conv(X,Y)
%
% This function computes the convolution of X and Y in the frequency domain
% and returns the result in C. Both X and Y must be one-dimensional (line
% or column) vectors. The format of the output C (i.e., line or column
% vector) matches that of the vector X.
%
% This function allows significant savings in execution time compared to
% the time-domain equivalent, i.e., Matlab's conv function.

% Release date: August 2008
% Author: Eric A. Lehmann, Perth, Australia (www.eric-lehmann.com)

xxs = size(xx); yys = size(yy);
if min(xxs)~=1 || min(yys)~=1,
    error('Both input vectors must be one-dimensional.');
end

xx = xx(:); yy = yy(:);
rlen = length(xx)+length(yy)-1;
rlen_p2 = 2^nextpow2(rlen);
XX = fft(xx,rlen_p2);
YY = fft(yy,rlen_p2);
rr = ifft(XX.*YY,'symmetric');
rr = rr(1:rlen);    %column vector

if xxs(1)==1,   % output rr in same format as xx
    rr = rr.';
end
