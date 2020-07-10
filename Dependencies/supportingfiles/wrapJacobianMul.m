function [X] = wrapJacobianMul(J, X, flag, JXfun, JTXfun ,multiplex)
% [JX] = wrapJacobian(J, X, id, JX, JTX, multiplex)
%
% Wrapper function to join separate Jacobian and jacobian transpose
% multiplication functions. 
% Usefull for validateJacobian and lsqnonlin.
%
% INPUTS:
%  J : Jacobian information
%  X : vector (/matrix) to be multiplied by J
%  flag : flag that selects multiplication with J (flag>0), J' (flag<0)  or J'*J (flag==0)
%  JXfun : function that multiplies J with X,  "J*X" = JXfun( J, X )
%  JTXfun : function that multiplies J' with X,  "J'*X" = JTXfun( J, X )
%  multiplex: boolean, default = true; if true JXfun and JTXfun support multiple columns in X, if false
%             JXfun and JTXfun are called for each column of X seperately. 
%
% Example usage:
%  JacMulFul = @(J, X, flag) wrapJacobianMul(J, X, flag, JXfun, JTXfun , multiplex)
%  validateJacobian( fun, x, JacMulFul )
%
% Copyright (C)2016  Dirk Poot, Erasmus Medical Center
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



if nargin<6 || isempty(multiplex)
    multiplex = true;
end;

if flag>=0 
    if multiplex || size(X,2)==1
        X = JXfun(J, X);
    else
        X2 = cell(1,size(X,2));
        for k=1:numel(X2)
            X2{k} = JXfun(J, X(:,k));
        end;
        X = [X2{:}]; % horizontal concatentation
    end;
end;
if flag<=0
    if multiplex || size(X,2)==1
        X = JTXfun(J, X);
    else
        X2 = cell(1,size(X,2));
        for k=1:numel(X2)
            X2{k} = JTXfun(J, X(:,k));
        end;
        X = [X2{:}]; % horizontal concatentation
    end;
end;