function out=maxloc(x,dim)
%function out=maxloc(x,dim)
%
%  Replicates the functionality of maxloc in fortran.
%
%  INPUTS:       x -> data variable    
%              dim -> (optional) dimension upon which to operate
%                     (same in Fortran and Matlab)
%
% OUTPUTS:     out -> set to second output of Matlab's max() if dim is specified
%                     otherwise, is the maximum location for the entire array
%                     returned as a cell array resulting from ind2sub
%                     this could then be used as a subscript list
%
% EXAMPLES:   a=magic(5); b=maxloc(a); a(b{:}),b
%                         b=maxloc(a,1)
%
% author: Ben Barrowes 3/2008, barrowes@alum.mit.edu


if nargin<2, dim=0; end

xdim=length(size(x));
if dim==0
 [dumvar,tempOut]=max(x(:));
 out=cell(1,xdim);
 [out{:}]=ind2sub(size(x),tempOut);
 %out=[out{:}];
else
 if ~isnumeric(dim), error('dim must be anumber in maxloc'); end
 if dim~=round(dim), error('dim must be an integer in maxloc'); end
 if dim>xdim || dim<1, error('dim must be 0<dim<=length(size(x))'); end
 [dumvar,out]=max(x,[],dim);
end

