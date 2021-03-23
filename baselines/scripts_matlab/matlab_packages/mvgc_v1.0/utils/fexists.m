function res = fexists(fhandle)

try
    feval(fhandle);
    %fprintf('no exception\n');
    res = true;
catch except
    %fprintf('exception: ''%s'', ''%s''\n', except.identifier,except.message);
    res = ~strcmpi(except.identifier,'MATLAB:UndefinedFunction');
end
