function out = parse_defaults_optionvaluepairs( defaults , varargin )
% out = parse_defaults_optionvaluepairs( defaults , option1, value1, ...)
% out = parse_defaults_optionvaluepairs( defaults , option_value_struct)
% Helper function to parse option value pairs and update the fields in the defaults structure.
%
% INPUTS:
% defaults : structure with default values for all possible options.
%            There is one special field: 'parse_options' with options
%            provided to parse_defaults_optionvaluepairs. This field is not
%            present in out. 
% option_i : string with name of options, case sensitive matched to the fields in 
%            defaults. Errors if the option is not found in the defaults structure.
%            For scalar structures in defaults, the subfields can be 
%            assigned by specifying 'option.sub_field_name' as option_i.
% value_i  : anything (that can be assigned to a field of a structure).
% option_value_struct : structure with valid fields (each field should also be in defaults)
%            The fields override the default values.
% 
% The following entries in 'parse_options' are supported:
%   case_sensitive : Default true; if false case insensitive matching of field names. 
%   convert_string : Default false; if true try to convert a string into the type of the default value. 
%                    For boolean default: 't', 'f', 'true', 'false', 
%                    For numeric default: str2double
%   obsolete       : structure (like defaults) with fields that are now obsolete. 
%                    The values in the elements of obsolete should be a scalar string or up to 3 element cell array.
%                    The scalar string or first element of the cell array should specify the new name
%                    that should be used instead. Use empty to silently ignore the provided option. 
%                    If present, the second element of the cell array should give a 
%                    description that is displayed when the field is used.
%                    If present the third element should contain a conversion function. 
%                    E.g.: 
%                    default.length_in_meter = 1;
%                    default.parse_options.obsolete.length_in_cm = {'length_in_meter', 'Parameter changed for better SI compatibility.',@(x) x/100 };
%                    out = parse_defaults_optionvaluepairs( default , 'length_in_cm',34 ) 
%                      Warning: Obsolete option field provided. Use length_in_meter instead of length_in_cm.
%                      Parameter changed for better SI compatibility.
%                    out.length_in_meter =0.34
%
% OUTPUTS:
% out      : structure defaults, with all fields specified by the options 
%            replaced by the corresponding values.
%            No checking of the values is performed.
%
%
% Copyrtight Dirk Poot, Erasmus Medical Center
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

% %          There is one special field 'addUnknownFields'.
% %   Reason not to include this option:
% %   (Performance and) We do not want to allow options that we do not know
% %   how to use/interpret as this will cause too many silent bugs (ie.
% %   mistyping an option will cause it to be ignored.)
% addUnknownFields = isfield(defaults,'addUnknownFields') && defaults.addUnknownFields;
% if addUnknownFields
%    defaults = rmfield(defaults,'addUnknownFields');
% end;

parse_options.case_sensitive = true;
parse_options.convert_string = false;
parse_options.obsolete = [];

nondefault_parse_options = isfield( defaults, 'parse_options' );
if nondefault_parse_options
    parse_options = parse_defaults_optionvaluepairs( parse_options, defaults.parse_options );
    out = rmfield( defaults, 'parse_options');
else
    out = defaults;
end;
if mod(nargin,2)~=1
    % if not odd number of input arguments => no well defined option-value pairs; maybe option_value_struct 
    if nargin==2 && isstruct(varargin{1}) && numel(varargin{1})==numel( out )
        % option values pairs are already in a scalar structure.
        fieldnm_opt = fieldnames(varargin{1});
        % Check if all options and match to fieldname in out: 
        if parse_options.case_sensitive
            fieldnm_out = fieldnm_opt;
            isopt = isfield(out, fieldnm_opt);
        else
            fieldnm_out = fieldnames( out );
            [isopt, idx ] = ismember( lower(fieldnm_opt ), lower( fieldnm_out ) );
            fieldnm_out = fieldnm_out( max(1,idx) );
        end;
        if any(~isopt)
            doerror = true;
            % Not all directly match. Check obsolete:
            if nondefault_parse_options && isstruct(parse_options.obsolete) 
                [fieldnm_out, nonmatch, varargin{1}, fieldnm_opt] = parse_obsolete( fieldnm_opt, fieldnm_out, ~isopt, varargin{1},  parse_options );
                doerror = any( nonmatch );
                isopt = ~nonmatch;
            end;
            if doerror
                allfieldnm_out = fieldnames(out);
                errorstring = disperror( fieldnm_opt(~isopt), allfieldnm_out, true);
                error( errorstring );
            end;
        end;
        

        % Do actual parsing:
        for k=1:numel(fieldnm_opt)
            curfieldnm_opt = fieldnm_opt{ k };
            curfieldnm_out = fieldnm_out{ k };
            if isstruct(out.(curfieldnm_out)) && isstruct(varargin{1}.(curfieldnm_opt))
                % if field of default is structure and field in options is structure, merge structures by calling self again:
                if nondefault_parse_options && ~isfield(out.(curfieldnm_out),'parse_options')
                    out.(curfieldnm_out).parse_options = parse_options;
                end;
                out.(curfieldnm_out) = parse_defaults_optionvaluepairs( out.(curfieldnm_out) , varargin{1}.(curfieldnm_opt) );
            else
                if nondefault_parse_options
                    out.( curfieldnm_out ) = assignhelper( out.( curfieldnm_out ) , varargin{1}.(curfieldnm_opt), parse_options);
                else
                    out.( curfieldnm_out ) = varargin{1}.( curfieldnm_opt );
                end;
            end;
        end;
        return;
    else
        error('arguments should come in option+value pairs.');
    end;
end;

% Actual option value pairs:
npairs = (nargin-1)/2;
isopt = [true(1,npairs);false(1,npairs)];
for k=1:npairs
    if ~ischar(varargin{k*2-1}) 
        isopt(1,k) = false;
        break;
    end;
end;
[allfieldnm , nnormal] = expandedfieldnames( out );
isopt_char = isopt;
if parse_options.case_sensitive
    [isopt(isopt), idx ] = ismember( varargin(isopt), allfieldnm);
else
    [isopt(isopt), idx ] = ismember( lower( varargin(isopt) ), lower( allfieldnm ) );
end;
matchedfieldname = allfieldnm( max(1,idx) );
haserror = ~all(isopt(1,:));
if haserror && nondefault_parse_options && isstruct(parse_options.obsolete) 
    % Not all directly match. Check obsolete:
    fieldnm_opt = varargin(1:2:end);
    [matchedfieldname, nonmatch, varargin] = parse_obsolete( fieldnm_opt, matchedfieldname, ~isopt(1,:), varargin,  parse_options );
    haserror = any( nonmatch );
    isopt(1,:) = ~nonmatch;
end;
if haserror
    k2= find( ~isopt(1,:) );
    invalidoptions = varargin( k2*2-1 );
    errorstring = disperror( invalidoptions, allfieldnm, false);
    error( errorstring );
    return; % actually not needed as an error is raised before.
end;

for k=1:numel(idx)
    if idx(k) <= nnormal
        if isstruct( out.( matchedfieldname{ k } ) )  && isstruct( varargin{k*2} )
            if nondefault_parse_options && ~isfield(out.( matchedfieldname{ k } ),'parse_options')
                out.( matchedfieldname{ k } ).parse_options = parse_options;
            end;
            out.( matchedfieldname{ k } ) = parse_defaults_optionvaluepairs( out.( matchedfieldname{ k } ) , varargin{k*2} );
        else
            if nondefault_parse_options
                out.( matchedfieldname{ k } ) = assignhelper( out.( matchedfieldname{ k } ) , varargin{k*2}, parse_options);
            else
                out.( matchedfieldname{ k } ) = varargin{k*2};
            end;
        end;
    else
        if nondefault_parse_options
            eval( ['out.' matchedfieldname{ k } '=assignhelper( out.' matchedfieldname{ k } ',varargin{k*2}, parse_options);']);
        else
            eval( ['out.' matchedfieldname{ k } '=varargin{k*2};']);
        end;
    end;
end;

function out = assignhelper( default, value, parse_options )
if parse_options.convert_string && ischar( value ) && ~ischar( default )
    if islogical( default )
        if isequal(lower(value),'true') || isequal(lower(value),'t') || isequal(str2double(value),1)
            value = true;
        elseif isequal(lower(value),'false') || isequal(lower(value),'f') || isequal(str2double(value),0)
            value = false;
        else
            error(['cannot convert ''' value ''' to a logical value']);
        end;
    elseif isnumeric( default ) 
        if numel(default)==1
            value = str2double( value) ;
        else
            error('currenly can only convert scalars.')
        end;
    else
        error('default value is of a type to which I cannot convert a string')
    end;
end;
out = value;

function [list, norig] = expandedfieldnames( s )
list = fieldnames( s );
norig = numel(list);
adlist = {};
for k=1:norig
    if isstruct(s.(list{k})) && numel(s.(list{k}))==1
        adlist{end+1} = expandedfieldnames( s.(list{k}) );
        for k2 = 1:numel(adlist{end})
            adlist{end}{k2} = [list{k} '.' adlist{end}{k2}];
        end;
    end;
end;
if ~isempty(adlist)
    list = vertcat(list,adlist{:});
end;

function errorstring = disperror( invalidoptions, alloptions, isstruct)
allfieldsstr = sprintf('"%s", ', alloptions{:});
alloptsstr = cell(1,numel(invalidoptions));
reploptstr = cell(1,numel(invalidoptions));
for k2= 1:numel(invalidoptions)
    if ischar(invalidoptions{k2})
        alloptsstr{k2} = ['"' invalidoptions{k2} '", '];
        if exist('levenshtein_distance.m','file')==2
            [dist] =  levenshtein_distance( invalidoptions{k2},alloptions,'caseCost',.2);
            replopt = ( dist <= 1.1*min(dist) );
            reploptstr{k2} = sprintf('"%s", ',alloptions{replopt});
        end;
    else
        alloptsstr{k2} = sprintf('class of option %d: %s, ',k2, class(invalidoptions{k2}));
    end;
end;
alloptsstr= [alloptsstr{:}];
reploptstr = [reploptstr{:}];
%     error('Each option should be a string identifying a valid option.\nAll valid options: %s\nProvided options : %s', allfieldsstr(1:end-2), alloptsstr(1:end-2));
if isstruct
    startstr = 'Each field of your option structure should identify a valid option.' ;
else
    startstr = 'Each option should be a string identifying a valid option.' ;
end;
if isempty(reploptstr)
    errorstring = sprintf( '%s\n\nProvided invalid options : %s\nAll valid options        : %s', startstr, alloptsstr(1:end-2), allfieldsstr(1:end-2) );
else
    errorstring = sprintf( '%s\n\nProvided invalid options : %s\nLikely valid options     : %s\n\nAll valid options        : %s', startstr, alloptsstr(1:end-2),reploptstr(1:end-2), allfieldsstr(1:end-2) );
end;

function [fieldnm_out, nonmatch, optval,  fieldnm_opt] = parse_obsolete( fieldnm_opt, fieldnm_out, nonmatch, optval,  parse_options )
fieldnm_obs = fieldnames( parse_options.obsolete ) ;
fieldnm_nonmatch = fieldnm_opt( nonmatch );
if ~parse_options.case_sensitive
    fieldnm_obs = lower(fieldnm_obs);
    fieldnm_nonmatch = lower(fieldnm_nonmatch);
end;
[isobs_opt, idx_obs ] = ismember( fieldnm_nonmatch, fieldnm_obs );
if any(isobs_opt)
    % at least one of the non direct matches is an obsolete field. Parse this and give 'Obsolete' warning. 
    nonmatch_idx = find(nonmatch);
    delete = false(size(nonmatch_idx));
    for id = find( isobs_opt )
        obsval = parse_options.obsolete.( fieldnm_obs{ idx_obs(id) } );
        if ischar(obsval)
            obsval = {obsval};
        end;
        % obsval = { new_field_name, message, converter function };
        adstr = '';
        if numel(obsval)>=2
            adstr = obsval{2};
        end;
        obsoletename = fieldnm_opt{ nonmatch_idx(id) } ;
        if isempty( obsval{1} )
            warning('ParseOptions:ObsoleteValue','Obsolete option field provided. Field ''%s'' is no longer available and is ignored. %s',obsoletename, adstr);
        else
            warning('ParseOptions:ObsoleteValue','Obsolete option field provided. Use ''%s'' instead of ''%s''. %s',obsval{1},obsoletename, adstr);
        end;
        if numel(obsval)>=3 %&& isa(obsval{3},'function_handle')
            % parse value of input:
            if isstruct(optval)
                optval.( obsoletename ) = obsval{3}( optval.( obsoletename ) );
            else
                optval{ nonmatch_idx(id) *2 } = obsval{3}( optval{ nonmatch_idx(id) *2 } );
            end;
        end;
        % set to write the (converted) value to the new field.
        % when case sensitive fieldnm_out matches to fieldnm_opt
        if isempty(obsval{1})
            delete( id ) = true;
        else
            fieldnm_out( nonmatch_idx( id  ) ) = obsval(1);
        end;
    end;
    nonmatch( nonmatch_idx( isobs_opt) ) = false; % we found obsolete fields that match. 
    if any(delete)
        deleteidx = nonmatch_idx( delete ) ;
        fieldnm_out(deleteidx)=[];
        nonmatch(deleteidx)=[];
        fieldnm_opt(deleteidx)=[];
    end;
end;
