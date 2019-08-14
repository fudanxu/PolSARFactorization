function T = C2T(C)
%
%   convert C to T
%
%   by Feng Xu
%   Fudan University, EMW Lab
%   fengxu@fudan.edu.cn
%   GNU licence

[p,q] = size(C);
if p==3 && q ==3
    U3 = [1,0,1;
          1,0,-1;
          0,sqrt(2),0];
    T = U3*C*U3'/2;
elseif p==4 && q==4
    U4 = [1,0,0,1;
          1,0,0,-1;
          0,1,1,0;
          0,1j,-1j,0];
    T = U4*C*U4'/2;
end