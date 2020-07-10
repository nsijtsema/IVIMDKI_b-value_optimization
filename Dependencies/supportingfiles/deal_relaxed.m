function [varargout] = deal_relaxed( varargin )
% [varargout] = deal_relaxed( varargin )
% This function like deal assigns inputs directy to outputs.
% It only avoids some of the 'limitiations' of deal
% when more outputs than inputs are requested, empty arguments are provided
% for the remaining outputs. Also more inputs than outputs may be given.
% The extra inputs are not used. 
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




if nargin<nargout
    varargout = cell(1,nargout);
    varargout(1:nargin)=varargin(1:nargin);
elseif nargin>nargout
    varargout = varargin(1:nargout);
else
    varargout = varargin;
end;