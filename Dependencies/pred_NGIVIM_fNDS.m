function [pred, dpred_dtht] = pred_NGIVIM_fNDS( tht, bvalues )
% [pred, dpred] = pred_IVIM( tht, grads )
%
% pred(i) = S0* (f * exp(-b(i) * Dstar) + (1-f)*exp(-b(i) * D + ((1/6)*b^2*D^2*K) )
%
% with:
%  tht = [S0;f;Dstar;D; K]
%
% INPUTS:
%   tht : parameters (to be estimated by fit_MRI) (theta)
%   bvalues : column vector with bvalues.   
%
% OUTPUTS:
%   pred : output array with predicted intensities.
%   dpred_dtht : derivative of pred with respect to tht. 
%
% 
% Copyright (C)2018 Nienke Sijtsema, Erasmus MC
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

%%
computegrad = nargout>1;

% Variable order: NOTE: matches with order in construction of dpred_dtht
S0idx = 1;  
fidx = 2;
Dstaridx = 3;
Didx = 4;
Kidx = 5;

S0 = tht(S0idx,:);
f = tht(fidx,:); 
D = tht(Didx,:);
Dstar = tht(Dstaridx,:);
K = tht(Kidx,:);

%Intialize some easy calculation shortcuts
EDstar = exp(-bvalues .* Dstar);
ED = exp(-bvalues .* D + (1/6) .* bvalues.^2 .* D.^2 .* K);

%Make prediction model
pred = bsxfun(@times, S0 .* f ,  EDstar ) + bsxfun(@times, S0 .* (1-f), ED );

if computegrad
      %Compute derivates with respect to each parameter.  
      dpred_dS0 = bsxfun(@times, f ,  EDstar ) + bsxfun(@times, (1-f), ED );
      
      dpred_df = bsxfun( @times, S0, EDstar) - bsxfun(@times, S0 ,  ED );
      %dpred_df = S0.*EDstar - S0 .* ED;
      
      D_a=  (((-1/3).*bvalues)*(S0.*(f-1))) .* (bvalues*(D.*K)-3) ; %NOTE: ordering to be correct for multi-theta.
      %D_a=(-1/3).*bvalues.*(bvalues.*D.*K-3).*S0.*(f-1); %mistake: was 1-f
      dpred_dD = bsxfun(@times, D_a, ED);
      
      dpred_dDstar = ((-bvalues)*(S0.*f)).* EDstar;
      %dpred_dDstar = bsxfun(@times, -bvalues.*S0.*f, EDstar);
      
      K_a = bvalues.^2 *(S0 .* D.^2 .*((1/6)-(f/6)));
      %K_a = S0 .* bvalues.^2 .* D.^2 .*((1/6)-(f/6));
      dpred_dK = bsxfun(@times, K_a, ED);
      
      
    % construct full derivatives matrix: (assumes a specific variable order)
    dpred_dtht = cat( 3, dpred_dS0, dpred_df, dpred_dDstar, dpred_dD, dpred_dK );
end;

function test
%% Test bi-exponential with CRLB evaluation
bvalues = [10 20 50 100 200 400 600 800 1000 1200]'/1000;
tht_tst = [ones(1,10);      % S0
           .1*ones(1,10);   % f
           6*ones(1,10);    % Dstar
           .2*ones(1,10);   % D
           2*ones(1,10)];   % K
fun3 = @(tht) pred_NGIVIM_fNDS( tht, bvalues);
pred = fun3( tht_tst );
sigma = .1;
[CRLB , I, J] = CramerRaoLowerBound_MRI( tht_tst, fun3, sigma);
plot(bvalues, pred,'*-')
sqrt( CRLB(I==J,:) )
validateJacobian( fun3, tht_tst )