function state = rng(s)

if     nargin == 0 % save
    
    assert(nargout == 1,'bad syntax');
    state.rand  = rand('twister');
    state.randn = randn('state'); 
    
elseif nargin == 1 % seed/restore
    
    assert(nargout == 0,'bad syntax');
    if isstruct(s)
        seed = s;
    else
        seed.rand  = s;
        seed.randn = s;
    end
    rand('twister',seed.rand);
    randn('state', seed.randn);

else
    error('bad number of arguments');
end
