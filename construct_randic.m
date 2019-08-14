function [kp, Tdic, pars]= construct_randic(ndic)
%
%   construct a random dictionary
%
%   by Feng Xu
%   Fudan University, EMW Lab
%   fengxu@fudan.edu.cn
%   GNU licence

kp = randsphere_complex(ndic,3)';

Tdic(1,:) = kp(1,:).*conj(kp(1,:));
Tdic(2,:) = kp(2,:).*conj(kp(2,:));
Tdic(3,:) = kp(3,:).*conj(kp(3,:));
Tdic(4,:) = sqrt(2)*real(kp(1,:).*conj(kp(2,:)));
Tdic(5,:) = sqrt(2)*imag(kp(1,:).*conj(kp(2,:)));
Tdic(6,:) = sqrt(2)*real(kp(1,:).*conj(kp(3,:)));
Tdic(7,:) = sqrt(2)*imag(kp(1,:).*conj(kp(3,:)));
Tdic(8,:) = sqrt(2)*real(kp(2,:).*conj(kp(3,:)));
Tdic(9,:) = sqrt(2)*imag(kp(2,:).*conj(kp(3,:)));

pars = zeros(6,size(kp,2));
for i=1:size(kp,2)
    pars(:,i) = kdeo(kp(:,i));
end