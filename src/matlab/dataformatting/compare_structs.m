function d = compare_structs(s1, s2, tol)

if ~exist('tol', 'var'), tol = 1e-10; end

fnames = union(fieldnames(s1), fieldnames(s2));
d = struct;

% loop over the fields
for i = 1:length(fnames)
  fname = fnames{i};
  
  % check if data is present
  if isfield(s1, fname)
    val1 = s1.(fname);
  else
    val1 = 'Data not present';
  end
  if isfield(s2, fname)
    val2 = s2.(fname);
  else
    val2 = 'Data not present';
  end

  % convert funcs to strings
  if isa(val1, 'function_handle')
    val1 = func2str(val1);
  end
  if isa(val2, 'function_handle')
    val2 = func2str(val2);
  end

  % recursion for structs and cells
  if isstruct(val1) && isstruct(val2)
    nested_d = compare_structs(val1, val2, tol);
    if ~isempty(fieldnames(nested_d))
      d.(fname) = nested_d;
    end
  
  % check equality
  else   
      if isequalntol(val1, val2, tol)
        continue
      else
        d.(fname) = {val1, val2};
      end
  end
end
end


%%
% ---- infile function -----

% Returns True if A and B are within tolerance, ignoring nans and infs
function ok = isequalntol(A,B,tol)

% data types are not equal
if ~isequal(class(A), class(B))
  ok = 0;
  return
end

% data sizes are not equal
if size(A) ~= size(B)
  ok = 0;
  return
end

% For non-numeric data types, check equality using the built-in isequaln
if ~isnumeric(A)

  % recursively enter cell arrays 
  % (TODO: recursively enter structs as well) 
  if iscell(A)
    if isempty(A) && isempty(B)
      ok = 1;
      return
    else
      for i = 1:numel(A)
        ok = isequalntol(A{i}, B{i}, tol);
      end
    end
  else
    ok = isequaln(A,B);
  end

% For numeric data types, apply special processing to handle nans, infs,
% and a numeric tolerance
else
  
  % ignore nans
  idx_nan1 = isnan(A);
  idx_nan2 = isnan(B);
  if ~isequal(idx_nan1, idx_nan2)
    ok = 0;
    return
  else
    A(idx_nan1) = 0;
    B(idx_nan2) = 0;
  end  
  
  % ignore infs (but preserve sign of inf)
  idx_inf1 = isinf(A);
  idx_inf2 = isinf(B);
  if ~isequal(idx_inf1, idx_inf2)
    ok = 0; 
    return
  else
    % preserve sign of inf, force larger than tol so that (-inf != inf)
    A(idx_inf1) = sign(A(idx_inf1)) * tol * 10;  
    B(idx_inf1) = sign(B(idx_inf1)) * tol * 10;
  end
  
  % compute relative difference
  rel_diff = (A(:) - B(:));
  
  % normalize by magnitude
  rel_diff = rel_diff / max([norm(A(:)), norm(B(:)), 1]);
  
  % check if close
  ok = all(abs(rel_diff) < tol);
end
end
