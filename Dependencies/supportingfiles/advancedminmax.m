function [mn,mx]=advancedminmax(Data)
% [min,max]=advancedminmax(Data);
% Computes both the minimum and maximum of Data. Intended to be faster than
% calling both min and max separately on a large matrix. Also, for complex
% Data 4 minima/maxima are returned: 
% [max(abs(Data)) max(angle(Data)) max(real(Data)) max(imag(Data))]
%
% Copyright (C)2009  Dirk Poot, University of Antwerp.
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


sz = size(Data);
if isreal(Data)
    %mx = max(Data);
    %mn = min(Data);
    if prod(sz)<1e5
        mx = max(Data);
        mn = min(Data);
    else
        if sz(1)<20000
            mx = zeros([1 sz(2:end)]);
            mn = zeros([1 sz(2:end)]);
            for k=1:prod(sz(2:end))
                mx(k) = max(Data(:,k));
                mn(k) = min(Data(:,k));
            end;
        else
            blkst = 1:10000:sz(1)-1000;
            blked = [blkst(2:end)-1 sz(1)];
            mx = zeros([numel(blkst) sz(2:end)]);
            mn = zeros([numel(blkst) sz(2:end)]);
            for k=1:prod(sz(2:end))
                for k2 = 1:numel(blkst)
                    datsel = Data(blkst(k2):blked(k2),k);
                    mx(k2,k) = max( datsel );
                    mn(k2,k) = min( datsel );
                end;
            end;
            mx = max(mx);
            mn = min(mn);
        end;
    end;
else % Complex Data
    %mx = [max(abs(Data)); max(angle(Data)); max(real(Data)); max(imag(Data))]
    %mn = [min(abs(Data)); min(angle(Data)); min(real(Data)); min(imag(Data))]
    if prod(sz)<1e5
        mx = [max(abs(Data)); max(angle(Data)); max(real(Data)); max(imag(Data))];
        mn = [min(abs(Data)); min(angle(Data)); min(real(Data)); min(imag(Data))];
    else
        if sz(1)<20000
            mx = zeros([4 sz(2:end)]);
            mn = zeros([4 sz(2:end)]);
            for k=1:prod(sz(2:end))
                datsel = Data(:,k);
                mx(:,k) = [max(abs(datsel)); max(angle(datsel)); max(real(datsel)); max(imag(datsel))];
                mn(:,k) = [min(abs(datsel)); min(angle(datsel)); min(real(datsel)); min(imag(datsel))];
            end;
        else
            blkst = 1:10000:sz(1)-1000;
            blked = [blkst(2:end)-1 sz(1)];
            mx = -inf([4 sz(2:end)]);
            mn = inf([4 sz(2:end)]);
            for k=1:prod(sz(2:end))
                for k2 = 1:numel(blkst)
                    datsel = Data( blkst(k2):blked(k2) , k);
                    absdt = abs(datsel);
                    angdt = angle(datsel);
                    mx(:,k) = max(mx(:,k),[max(absdt); max(angdt); max(real(datsel)); max(imag(datsel))]);
                    mn(:,k) = min(mn(:,k),[min(absdt); min(angdt); min(real(datsel)); min(imag(datsel))]);
                end;
            end;
        end;
    end;
end;