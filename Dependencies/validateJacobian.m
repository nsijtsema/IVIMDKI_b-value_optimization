function [out] = validateJacobian(fun, x, jacobmulfun)
% validateJacobian( fun, x [, jacobmulfun])
% Function that validates the jacobian of the function 'fun' in x.
% Also supports model functions for fit_MRI.
% This can be used to check for errors in the jacobian computation.
% NOTE: this routine is slow as it numerically evaluates the jacobian.
%       Only use for debugging, and try to minimize the problem size.
%
% Evaluates [f, J] = fun(x)
% and checks if the numerical jacobian of fun around x is close enough 
% to the numerical jacobian J.
% If jacobmulfun is provided, it is used to compute the jacobian explicitly.
%    J = jacobmulfun(J, eye ,1); 
%    J'= jacobmulfun(J, eye ,-1); 
%   Additionally J and J' are compared.
%   Use wrapJacobianMul if you have separate functions for multiplication 
%   with J and J', call as: 
%      validateJacobian( fun, x, @(J,X,flag) wrapJacobianMul(J, X, flag, Jmul, JTmul ) )
% 
% NOTE:
% This function uses jacobianest from the 'DERIVESTsuite' by John D'Errico:
% http://www.mathworks.com/matlabcentral/fileexchange/13490-automatic-numerical-differentiation
%
% Copyright (C) Dirk Poot, Erasmus MC, 2-2-2012
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


maxJacCols = 1000; % Maximum number of variables for which the full jacobian is numerically evaluated. 
testJacCols = 100; % Number of test directions used if we don't evaluate the full jacobian. (for jacobian multiply function and numerical evaluation of jacobian)
maxJacRows = 1000; % Maximum number of outputs for which the full jacobian is numerically evaluated. 
testJacRows = 101; % Number of test directions used if we don't evaluate the full jacobian. (for jacobian adjoint multiply function )
[f, J] = fun(x);

reduceCols = numel(x) > maxJacCols;
reduceRows = numel(f) > maxJacRows;
if reduceCols
    r = randn(numel(x) , testJacCols );
    [colMatrix, R] = qr( r, 0 );
else
    colMatrix = eye(numel(x));
end;
    
if nargin>=3 && ~isempty(jacobmulfun)
    if reduceRows
        r = randn(numel(f) , testJacRows );
        [rowMatrix, R] = qr( r, 0 );
    else
        rowMatrix = eye(numel(f));
    end;
    % check jacobmulfun:
    J1 = jacobmulfun(J, colMatrix, 1);
    if size(J1,1)~=numel(f)
        error('The size of the output of the jacobian multiplication function (%d rows) does not correspond to the number of elements in the output of fun (%d elements)', size(J1,1), numel(f))
    end;
    JT = jacobmulfun(J, rowMatrix, -1);
    if size(JT,1)~=numel(x)
        error('The size of the output of the adjoint jacobian multiplication function (%d rows) does not correspond to the number of elements in x (%d elements)', size(JT,1), numel(x))
    end;
    if ~isreal(J1)
        JT = complex(JT, jacobmulfun(J, -1i*rowMatrix, -1));
    end;
    if reduceRows
        J1red = rowMatrix'* J1;
    else
        J1red = J1;
    end;
    if reduceCols
        JTred = colMatrix'* JT;
    else
        JTred = JT;
    end;
    Jscale = sqrt(dot(J1red(:),J1red(:))/numel(J1red));
    if max(max(abs(J1red-JTred')./max(abs(J1red+JTred'), Jscale*1e-7)))>1e-7
        % imagebrowse( cat( 3, J1red, JTred', J1red-JTred'));
        error('J and J'' of Jacobmulfun seem to differ.');
    end;
    r = randn(numel(x),1);
    rH = jacobmulfun(J, r, 0);
    if reduceRows || reduceCols
        rHs = jacobmulfun( J , jacobmulfun( J, r ,1 ),-1 );
    else
        rHs = real( JT * (J1 * r) );
    end;
    if max(max(abs(rH-rHs),[],2)./max(abs(rH+rHs),[],2))>1e-6
        error('J''*J of Jacobmulfun seem to differ from multiplying with J'' and J separately.');
    end;
    J = J1;
else
    if size(J,3)==size(x,1) && size(J,2)==size(x(:,:),2) %&& size(J,1)*size(J,2)==numel(x)
        % Assume multi time series fit (as for fit_MRI fitting functions)
        Jf = zeros(size(J,1)*size(J,2), size(J,2)*size(J,3));
        strowid = 1;
        stcolid = 1;
        for k=1:size(x,2)
            edrowid = strowid + size(J,1)-1;
            edcolid = stcolid + size(x,1)-1;
            Jf(strowid:edrowid, stcolid:edcolid) = squeeze(J(:,k,:));
            strowid = edrowid+1;
            stcolid = edcolid+1;
        end;
        J= Jf;
    end;

    if reduceCols
        J = J * colMatrix; % unlikely to be used. If Jacobian is large, a jacobian multiplication function should be implemented instead of explicitly returning the jacobian.
        warning('validateJacobian:Large','The jacobian is large because there are many variables in this function, consider implementing a jacobian multiplication routine instead of explicitly evaluating the jacobian.');
    end;
end;

% Using jacobianest from the 'DERIVESTsuite' by John D'Errico
% http://www.mathworks.com/matlabcentral/fileexchange/13490-automatic-numerical-differentiation
if reduceCols
    fun_red = @(x_red) fun( reshape( colMatrix * x_red, size(x))+x ) ;
    x_red = zeros( size(colMatrix,2) , 1);
    [Jest , errest] = jacobianest(fun_red, x_red);
else
    [Jest , errest] = jacobianest(fun, x);
end;
% if ~isreal(Jest) % original jacobianest gives wrong sign to imaginary part.
%     Jest = conj(Jest);
% end;

% sometimes errest is zero for some elements; this obviously cannot be true
% and causes division by zero in the next comparison. Therefore, adjust
% errest so that each value in each colum is at least 10% of the mean error
% in that column, and at least .1% of the mean over all columns.
errestm = errest +.1*ones(size(errest,1),1)*mean(abs(errest),1)+.001*mean(abs(errest(:)));
% sometimes Jest ==0 and errest ==0; that indicates some kind of failure of
% the jacobian estimation routine to estimate the jacobian properly. I
% don't really want to, but for now it's probably best to just ignore those
% entries. (and assume other checks will actually pick up any error in
% the analytical jacobian implementation)
differel = abs(Jest-J)./(errestm) >10 & abs(J-Jest)./abs(J+Jest)>1e-4 & ~(Jest==0 & errest==0);
if any(any( differel ))
    [rowIdx, colIdx] = find( differel );
    uRowIdx = unique(rowIdx);
    uColIdx = unique(colIdx);
    messagestr = sprintf( 'Numerical and analytic jacobians seem to differ in %d of %d elements.\nSpecifically in %d/%d (reduced) variables: [%s]\nand %d/%d outputs: [%s]', nnz(differel) , numel(J), numel(uColIdx), size(J,2), num2str(uColIdx'), numel(uRowIdx), size(J,1), num2str( uRowIdx' ) );
    if nargout>=1
        warning(messagestr);
    else
        error(messagestr);
    end;
end;
if any(all(Jest==0 & errest==0,1))   
    warning('Numerical jacobian not estimated properly; analytical jacobian not properly verified.');
end;
if nargout>=1
    out.f = f;
    out.J = J;
    out.Jest = Jest;
    if reduceCols
        out.colMatrix = colMatrix;
    else
        out.colMatrix = 'eye';
    end;
    out.errorest = errest;
    out.significantdifferences = differel;
else
    disp('The jacobian appears to be computed correctly.');
end;

function test_validateJacobian_complex
%% Test if validateJacobian treats complex jacobians correctly. Assume x is real and Jx is complex. 
% See complex value as real and imaginary pair that each have entries in the jacobian. 
npar = 3;
Jm = randn(8,npar);
fn = @(x) deal_relaxed( complex( Jm(1:4,:)*x, Jm(5:end,:)*x), Jm );
Jmul = @(J, x) complex( J(1:4,:)*x, J(5:end,:)*x ), 
JTmul = @(J,x) J'*[real(x);imag(x)]
validateJacobian(fn, randn(npar,1), @(J,X,flag) wrapJacobianMul(J, X, flag, Jmul, JTmul) )