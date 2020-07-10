function [xout , state] = HaltonSequence( npnt, base , x0 )
% [xout , final_state] = HaltonSequence( npnt, base , x0 )
% OR: [xout , final_state] = HaltonSequence( npnt, initial_state )
% Creates a halton sequence starting from an arbitrary point 0<=x0<1.
% This sequence is a pseudorandom sequence. 
%
% Cite as Wang2000: Wang X, Hickernell FJ: Randomized Halton sequences. Math Comput Model 2000; 32:887�899.
%
% INPUTS:
%   npnt : scalar integer that specifies the number of points to generate.
%   base : A postive integer base for the Halton sequence. To 'uniformly' fill a
%          multi-dimensional space use a different base for the coordinate
%          in each dimension. The bases should be relative prime and small.
%   x0   : starting point (default =0). 
%   initial_state : second output of a previous call to HaltonSequence.
%                   Continue with that sequence.
%
% OUTPUT: 
%   xout : pseudo random sequence; numel(x0)  x npnt
%          with values in the range 0 .. 1
%   final_state : final state which can be used as input to a subsequent
%                 call to HaltonSequence to proceed with the sequence.
%
% Copyright (C)2013  Dirk Poot, TU Delft
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.



if nargin<3
    x0 = 0;
else
    x0 = mod(x0,1);
end;
if isstruct(base)
    state = base;
    base = state.base;
    digits = state.digits;
else
    maxdigits = max( -log(eps)./log(base) );
    mindigits = ceil( log(npnt)./min(log(base)) );
    % decompose x0 into digits
    digits = zeros(numel(x0), mindigits);
    k = 1;
    while any(x0)>0 && k<maxdigits;
        x1 = x0.*base;
        digits(:,k) = floor(x1);
        x0 = x1-digits(:,k);
        k = k+1;
    end;
end;
resid = (base-1-digits)*base.^(0:size(digits,2)-1)';
if resid<npnt
    lastdigit = ceil( log(npnt+base.^size(digits,2)-resid)/log(base) );
    digits(1,end+1:lastdigit)=0;
end;
toValue = base.^(size(digits,2)-1:-1:0)'; % should be all integers.
scale = base.^size(digits,2);
xout = zeros(numel(x0),npnt);
for k=1:npnt
    xout(:,k) = (digits * toValue)/scale; % roundoff only in division by (the integer) scale
    % rightward carry addition, as eq 7 in Wang2000. 
    admsk = ones(size(xout,1),1);
    i = 1;
    while any(admsk) 
        digits(:,i) = digits(:,i)+admsk;
        admsk = digits(:,i)>=base;
        digits(admsk,i)=0;
        i=i+1;
    end;
end;
if nargout>1
    state.base = base;
    state.digits = digits;
end;
    