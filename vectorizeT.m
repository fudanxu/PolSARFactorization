function X = vectorizeT(T)
%
%   vectorize 3*3*n*m T to 9*(n*m) X matrix
%   of inverse transform X back to T
%

[q,p,n,m] = size(T);

if q==3 && p==3
    X = zeros(9,n*m);
    X(1,:) = T(1,1,:);
    X(2,:) = T(2,2,:);
    X(3,:) = T(3,3,:);
    X(4,:) = sqrt(2)*real(T(1,2,:));
    X(5,:) = sqrt(2)*imag(T(1,2,:));
    X(6,:) = sqrt(2)*real(T(1,3,:));
    X(7,:) = sqrt(2)*imag(T(1,3,:));
    X(8,:) = sqrt(2)*real(T(2,3,:));
    X(9,:) = sqrt(2)*imag(T(2,3,:));    
elseif q==9 && n==1 && m==1
    X2 = T; %input is X matrix
    T2 = zeros(3,3,p);
    T2(1,1,:) = max(X2(1,:),0);
    T2(2,2,:) = max(X2(2,:),0);
    T2(3,3,:) = max(X2(3,:),0);
    T2(1,2,:) = (X2(4,:)+1i*X2(5,:))/sqrt(2);
    T2(2,1,:) = (X2(4,:)-1i*X2(5,:))/sqrt(2);
    T2(1,3,:) = (X2(6,:)+1i*X2(7,:))/sqrt(2);
    T2(3,1,:) = (X2(6,:)-1i*X2(7,:))/sqrt(2);
    T2(2,3,:) = (X2(8,:)+1i*X2(9,:))/sqrt(2);
    T2(3,2,:) = (X2(8,:)-1i*X2(9,:))/sqrt(2);
    X = T2;
else
    error('incorrect dimension');
end