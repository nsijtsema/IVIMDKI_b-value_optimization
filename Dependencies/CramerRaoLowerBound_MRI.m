function [CRLB,Il,Jl] = CramerRaoLowerBound_MRI(par, fun, sigma, imagetype, fields, Iprior)
% [CRLB,I,J] = CramerRaoLowerBound_MRI(par, fun, sigma, imagetype, fields, Iprior)
% Computes the Cramer Rao Lower bound of the variance of the parameters par
% of the model specified by fun.
% noiseLvl: default = 1; Specifies the noise standard deviation of the
%    images (sigma in the rice distribution), the resulting CRLB depends on
%    the noise level. 
%
% The resulting CRLB is stored efficiently. For each column of par (i.e. for each voxel), 
% only the upper triangular part of the Cramer Rao matrix is stored.
% restore full CramerRao matrix (covariance matrix) of voxel k by:
% mat = zeros(size(par,1),size(par,1))
% for l = 1 : numel(I); mat( I(l) , J(l)) = CRLB(l,k); mat( J(l) , I(l)) = CRLB(l,k); end;
%
% Also see FisherInformation_MRI, of which the code is(/should be) almost the same.
% The Cramer Rao Lower Bound is the inverse of the Fisher information matrix. 
%
% INPUTS:
% par : N dimensional matrix. Each par(:,k) specifies the parameter vector for one voxel
% fun : function that predicts the MR signal intensity, given a parameter vector par
%       Assumes fun is vectorized (i.e. multiple columns of par can be given at once). 
%       Second output of fun is the derivative of the predicted magnitudes w.r.t. the
%       parameters.
%       So fun is called as follows:
%         [A, dAdpar] = fun( par(:, index_range) )
%       A      = nMRI x numel(index_range) predicted magnitudes
%       dAdpar = nMRI x numel(index_range) x size(par,1) derivative of each A w.r.t. each parameter.
% sigma : default = 1; specifies noise level. Any of:
%                   - scalar   : all voxels in all images have same noise level
%                   - 1 x size(par,2..) : specified per voxel, 1 noiselevel for all images, 
%                   - nMRI x 1 : 1 per image, same for all voxels 
%                   - nMRI x size(par,2..) : noise level specified for each measurement separately (how can you estimate that..). 
%                   - 'lastInTheta' : par(end,:) specifies the noise level for each voxel
% imagetype : {'magnitude'}, 'real'
% fields   : r x size(par,2) matrix with extra parameters. E.g. a background field.
%            if specified and non empty, fun should accept 2 arguments fun( par(:,range) , field(:, range))
% Iprior  : Prior Fisher information on parameters. 
%           Diagonal entries specify inverse a-priori variance.
%           If you have a constrained optimization with around par:  A*delta_par = 0 
%            Iprior = A'*A * c , with c tending to infinity, but small enough to avoid numerical problems.
%
% OUTPUT:
% CRLB : compactified Cramer Rao Lower bound. Each column is the CRLB of one column of par.
% I, J : CRLB(k, l) specifies mat(I(k),J(k)) and mat(J(k),I(k)) of voxel l
%
% Copyright (C) Dirk Poot (d.poot@erasmusmc.nl) Erasmus Medical Center.

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


if nargin<3 || isempty(sigma)
    sigma = 1;
end;
if nargin<4 || isempty(imagetype)
    imagetype = 'magnitude';
end;
sigmaInTheta = isequal( sigma ,'lastInTheta');
if nargin<6 
    Iprior = [];
end;
% Cramer Rao:
% Cov(DT) >= inv(I(DT))
% I(DT)_ij = E[ d{ logLikelihood(x,DT) }/{d(DT_i)} * d{ logLikelihood(x,DT) }/{d(DT_j)} ] = 
%          = E[ (d( logLikelihood(x,A)/d(A)) ^2  * d(A)/d(DT_i) * d(A)/d(DT_j) ]
%          = d(A)/d(DT_i)*d(A)/d(DT_j) * E[ (d( logLikelihood(x,A)/d(A))^2 ]
% E[ (d( logLikelihood(x,A)/d(A))^2 ] = int( ricepdf(x,A) * (d(logricepdf(x,A))/d(A))^2), x= 0.. infinity);

szpar = size(par);
numtr = prod(szpar(2:end));
szSigma = size(sigma);
npar = szpar(1);
mat = triu(ones(npar,npar));
[ Il , Jl ] = find(mat);
if sigmaInTheta
    [ I , J ] = find(mat(1:end-1,1:end-1));
    nparr = npar-1;
else
    I = Il;
    J = Jl;
    nparr = npar;
end;
matset1 = I+(J-1)*npar;
matset2 = J+(I-1)*npar;
matread = Il+(Jl-1)*npar;
CRLB = zeros( [numel(Il), szpar(2:end)] );

progressthreshold = ceil( 10000000/ size(par,1)^2 ) ;
showprogbar = numtr>progressthreshold;
if showprogbar
    progressbar('start',numtr,['Computing CRLB of ' num2str(numtr) ' elements'] );
end;

maxnuminblock = 1000;
selSigma = isequal(szSigma(2:end),szpar(2:end)) && (maxnuminblock<numtr);
stindx = 1;
while stindx<=numtr
    edindx = min(numtr, stindx + maxnuminblock);
    numtrinblock = edindx-stindx+1;
    if sigmaInTheta
        parsel = par(1:end-1,stindx:edindx);
        sigmak = par(end,stindx:edindx);
    else
        parsel = par(:,stindx:edindx);
    end;
    if nargin>=5  && ~isempty(fields)
        [A, dAdpar] = fun( parsel , fields(:,stindx:edindx));
    else
        [A, dAdpar] = fun( parsel );
    end;
    if sigmaInTheta
        if isequal(imagetype ,'magnitude')
            [expectIA, expectd2IAs, expectd2Iss] = riceExpectFuncs_s( A , sigmak) ;
        elseif isequal(imagetype ,'real') || isequal(imagetype, 'complex') % assume normal distribution.
            error('Real and complex images not supported when sigma is last element of theta.');
        else 
            error('wrong imagetype');
        end;
    else
        if ~isequal(sigma,1)
            if selSigma
                sigmak = sigma(:,stindx:edindx);
            else
                sigmak = sigma;
            end;
            A = bsxfun(@times, A, 1./sigmak);
        end;
        if isequal(imagetype ,'magnitude')
            expectIA = rice_information_ExpectVal(A) ;
        elseif isequal(imagetype, 'complex')
            if szSigma(1)==size(A,1)
                sigmak = [sigmak;sigmak];
            end;
            A =[real(A);imag(A)];
            dAdpar = [real(dAdpar);imag(dAdpar)];
            expectIA = ones(size(A));
        elseif isequal(imagetype ,'real') % assume normal distribution.
            expectIA = ones(size(A));
        else 
            error('wrong imagetype');
        end;
        if ~isequal(sigma,1)
            expectIA = bsxfun(@times, expectIA , 1./ sigmak.^2);
        end;
    end;
    dAdpar = permute( reshape(dAdpar,[ size(A,1) , numtrinblock, nparr]),[1 3 2]);
    rep = ones(1,numel(I));
%     doFIMinparts = size(dAdpar,1)*numel(I)>1e6;
    ws = warning('off');
    for k=1:numtrinblock
    %    [A,dAdpar] = DifusionTensor_Apred_m(DT(:,k),gradmdf);
    %    expectIA = rice_information_ExpectVal(A);
    %    FIM = sum(dAdpar(:,I).*dAdpar(:,J).*expectIA(:,ones(1,28)))';
        FIM = sum(dAdpar(:,I,k).*dAdpar(:,J,k).*expectIA(:,k*rep))';
        mat( matset1 ) = FIM(:);
        mat( matset2 ) = FIM(:);
        if sigmaInTheta
            Its = expectd2IAs(:,k)'*dAdpar(:,:,k);
            Iss = sum(expectd2Iss(:,k));
            mat(:,end) = [Its';Iss]; 
            mat(end,1:end-1)=Its;
        end;
        if ~isempty(Iprior)
            if size(Iprior(:,:,:),3)==numtr
                mat = mat + Iprior(:,:,k+stindx-1);  
            else
                mat = mat + Iprior;  
            end;
        end;
        mat = inv(mat);
        CRLB(:,k+stindx-1) = mat( matread );
    end;
    warning(ws);
    stindx = edindx+1;
    if showprogbar
        progressbar(stindx);
    end;
end;
if showprogbar
    progressbar('ready');
end;