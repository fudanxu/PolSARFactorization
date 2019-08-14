function X = randsphere_complex(m,n)
%
%   generate a complex 3D vector
%
%   n-D vector
%
X = randn(m,n)+1i*randn(m,n);
X = X./repmat(X(:,1),1,n);  %subtract absolute phase
X = X./repmat(sqrt(sum(abs(X).^2,2)),1,n);