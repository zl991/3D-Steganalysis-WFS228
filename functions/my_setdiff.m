function a = my_setdiff(a,b)
% MYSETDIFF Set difference of two sets of positive integers (much faster than built-in setdiff)
% C = my_setdiff(A,B)
% C = A \ B = { things in A that are not in B }

% 
    bits = true(1,max(a));

    bits(b) = false;
    a = a(bits(a));

% if isempty(A)
%     C = [];
%     return;
% elseif isempty(B)
%     C = A;
%     return; 
% else % both non-empty
%     bits = zeros(1, max(max(A), max(B)));
%     bits(A) = 1;
%     bits(B) = 0;
%     C = A(logical(bits(A)));
end