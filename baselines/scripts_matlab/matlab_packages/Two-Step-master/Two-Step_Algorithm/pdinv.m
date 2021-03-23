function Ainv = pdinv(A);

% PDINV Computes the inverse of a positive definite matrix
% Copyright (c) 2010-2011  ...
% All rights reserved.  See the file COPYING for license terms.
numData = size(A, 1);
try
  U = chol(A);
  invU = eye(numData)/U;
  Ainv = invU*invU'; 
catch
  [void, errid] = lasterr;
  if strcmp(errid, 'MATLAB:posdef')
    warning(['Matrix is not positive definite in pdinv, inverting' ...
	     ' using svd'])
    [U, S, V] = svd(A);
    Ainv = V*diag(1./diag(S))*U';
    return
  else
    error(lasterr)
  end
end

