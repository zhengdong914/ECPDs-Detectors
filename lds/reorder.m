% function Y=reorder(X,m);
% reorder the columns of X.
%
% If X is a n by k matrix, Y will be an n by k matrix with the columns
% reordered modulo m.
%
% for example
% reorder([1 2 3 4 5 6 7 8 9 10 11 12],4)
% ans =
%
%  1   5   9   2   6  10  3   7  11   4   8  12
%

function Y=reorder(X,m);

[n,k]=size(X);

if rem(k,m)~=0
  disp('Error in reorder: X must have a multiple of m rows');
  return;
end;
l=k/m;

if (m==1 | m==k)
  Y=X;
else
  i=1:k;
  j=rem(i-1,m)*l+fix((i-1)/m)+1;
  Y(:,j)=X(:,i);
end;