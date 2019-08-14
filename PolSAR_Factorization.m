%
%   Over-complete dictionary assembly and 
%   simutanously pixel clustering via NMF
%
%   by Feng Xu
%   Fudan University, EMW Lab
%   fengxu@fudan.edu.cn
%   GNU General Public License v3.0
%
%   Notation:
%   C: Input Covariance Matrix
%   T: Coherency Matrix
%   H: Height of the image
%   W: Width of the image
%   K: Size of the dictionary of the atom scatterers
%   U: over-complete dictionary
%   V: distribution maps
%
%   
%   Reference:
%   Feng Xu, Qian Song, and Ya-Qiu Jin,
%   "Polarimetric SAR Image Factorization", 
%   IEEE Transactions on Geoscience and Remote Sensing, 
%   Vol. 55, No. 9, pp. 5026-5041, 2017.  

close all
clear 

%% Read in data and data preprocessing
%  If the image size is too large, 
%  downsampling operation is conducted in advance.
disp('Read in data');
load C_TestData2.mat;
[~,~,H,W] = size(C);
for i=1:H
    for j=1:W
        T0(:,:,i,j) = C2T(C(:,:,i,j));
    end
end
clear C;

HDSRate = 1; WDSRate = 1;
if (H > 300)
    HDSRate = floor(H/300);
end
if (W > 300)
    WDSRate = floor(W/300);
end

T = T0(:,:,1:HDSRate:H,1:WDSRate:W);
[~,~,Hsub,Wsub] = size(T);

Tdata = vectorizeT(T);
[~,~,n,m] = size(T);
npix = n*m;

%% Construct dictionary for decomposition
disp('Construct dictionary for decomposition');
rng(568)
[kp, Tdic, pars] = construct_randic(100);
ndic = size(kp,2);
Tspan = sum(Tdata(1:3,:),1);
Tdata = Tdata ./ (ones(9,1) * Tspan);
gsum = Tspan;

%% Redundant coding matrix
disp('Generate redundant coding matrix');
RCM = zeros(ndic,npix);
for i=1:npix
    T = vectorizeT(Tdata(:,i));
    kpTinv = kp'/T;
    pdf = exp(-real(sum(kpTinv.*kp.',2)));
    RCM(:,i) = pdf/sum(pdf);
    disp(i);
end

%% NMF learning
disp('NFM learning');
options = [];
% options.maxIter = 1000;
options.alpha = 0;
options.nRepeat = 5;
%when alpha = 0, GNMF boils down to the ordinary NMF.
K = 8;
rng(8900);
tic
[U,V, nIter_final, objhistory_final] = GNMF_KL(RCM,K,[],options); %'
% options.WeightMode = 'Binary';
% W = constructW(RCM',options);
% options.maxIter = 100;
% options.alpha = 100;
% [U,V] = GNMF_KL(RCM,K,W,options); %'
toc

%normalize U
Unorm = sum(U,1);
U = U./(ones(ndic,1)*Unorm);
V = V.*(ones(npix,1)*Unorm);

disp(['Obj Function = ',num2str(objhistory_final)]);
snr = -10*log10(mean(mean((RCM-U*V').^2))/mean(mean(RCM.^2)));
disp(['Residue SNR = ',num2str(snr)]);

%assemble T matrix for each cluster
Tclur = Tdic*U;
Trecover = Tdic*U*V';

save U.mat U;

%% Generate the distribution map of the original data (without downsampling)
if (HDSRate>1 || WDSRate>1)    
    disp('Generate distribution map of the original image');
    RCM = reshape(RCM, ndic, Hsub, Wsub);
    for i=1:ndic
        RCMi = reshape(RCM(i,:,:), Hsub, Wsub);
        RCM0(i,:,:) = kron(RCMi, ones(HDSRate, WDSRate));
    end
    RCM0 = reshape(RCM0, ndic, H*W);
    
    Ttmp = T0;
    T0 = vectorizeT(Ttmp);
    T0 = reshape(T0, 9, H*W);
    V0 = zeros(H, W, K);
    V = reshape(V, Hsub, Wsub, 8);
    for i=1:K
        Vi = V(:,:,i);
        V0(:,:,i) = kron(Vi, ones(HDSRate, WDSRate));
    end
    V0 = reshape(V0, H*W, K);
    save V.mat V0;
    
    options = [];
    options.alpha = 0;
    options.nRepeat = 1;
    options.origImg = 1;
    K = 8;
    rng(8900);
    tic
    [U1, V1, nIter_final, objhistory_final] = GNMF_KL(RCM0, K, [], options);
    toc
    
    Unorm = sum(U1,1);
    U1 = U1./(ones(ndic,1)*Unorm);
    V1 = V1.*(ones(H*W,1)*Unorm);
end

%% Result illustration
figure;
for i=1:K
    subplot(2,4,i);
    visualizeT([],Tclur(:,i));
end

V = reshape(V, Hsub*Wsub, 8);
figure;
for i=1:K
    subplot(2,4,i);
    mesh((flipud(reshape(V(:,i), Hsub, Wsub)))); axis([1 Wsub 1 Hsub]); axis equal tight
%     imagesc(10*log10(reshape(V(:,i)'.*gsum,[n,m])));caxis([-30,20]);axis equal tight
end
% Wormap(flipud(Wormap(bone)))

figure;
for i=1:9
    subplot(3,3,i);
    imagesc(reshape(Tdata(i,:).*gsum,[n,m]));axis equal tight
    if i<4
        caxis([0,2]);
    else
        caxis([-1,1]);
    end
end
% Wormap(flipud(Wormap(bone)))

figure;
[P,Tr,Tg,Tb] = imrgb(reshape(Tdata(2,:).*gsum,[n,m]),reshape(Tdata(3,:).*gsum,[n,m]),reshape(Tdata(1,:).*gsum,[n,m]));
figure;
[P,Tr,Tg,Tb] = imrgb(reshape(Trecover(2,:).*gsum,[n,m]),reshape(Trecover(3,:).*gsum,[n,m]),reshape(Trecover(1,:).*gsum,[n,m]));

% The original image
if (HDSRate>1 || WDSRate>1)
    figure;
    for i=1:K
        subplot(2,4,i);
        mesh((flipud(reshape(V1(:,i), H, W)))); axis([1 W 1 H]); axis equal tight
    end
end
