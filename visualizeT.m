function visualizeT(kp,T)
%
%   plot T's PDF on the Extended alpha-beta hemisphere
%
if isempty(kp)
    %generate the standard mesh
    if size(T,1) == 9
        T = vectorizeT(T);
    elseif ~(size(T,1) == 3 && size(T,2) == 3)
        error('unrecognized input T');
    end
    
    n = 1000;
    [x,y,z] = sphere(n);
    x = x(n/2+1:n+1,:);
    y = y(n/2+1:n+1,:);
    z = z(n/2+1:n+1,:);
    kp = [z(:)';x(:)';y(:)'];

    %estimate pdf
    kpTinv = kp'/(T+diag([1,1,1])*1e-4);
    pdf = real(sum(kpTinv.*kp.',2)).^(-4);%exp(-real(sum(kpTinv.*kp.',2)));
    pdf = pdf/sum(pdf);
    
    cpdf = fix(127*(pdf-min(pdf))/(max(pdf)-min(pdf)))+1;
    cpdf = reshape(cpdf,size(x));
    surf(x,y,z,cpdf);shading interp
    colormap((colormap(jet(128))));
    view(0,90);xlim([-1 1]);ylim([-1 1]);zlim([0 1])
    axis equal tight; axis off;
    hold on;
    plot3(x(251,:),y(251,:),z(251,:),'color',[1,1,1],'linewidth',2)

%     n = 50;
%     [x,y,z] = sphere(n);
%     x = x(n/2+1:n+1,:);
%     y = y(n/2+1:n+1,:);
%     z = z(n/2+1:n+1,:);
%     mesh(x,y,z,ones(size(z))*64)
    
else
    
    [~,ndic] = size(kp);
    if numel(T) == ndic
        pdf = T;
    else
        if size(T,1) == 9
            T = vectorizeT(T);
        elseif ~(size(T,1) == 3 && size(T,2) == 3)
            error('unrecognized input T');
        end
        %estiamte pdf
        kpTinv = kp'/(T+diag([1,1,1])*1e-4);
        pdf = real(sum(kpTinv.*kp.',2)).^(-4);%exp(-real(sum(kpTinv.*kp.',2)));
        pdf = pdf/sum(pdf);
    end

    n = 24;
    [x,y,z] = sphere(n);
    x = x(n/2+1:n+1,:);
    y = y(n/2+1:n+1,:);
    z = z(n/2+1:n+1,:);

    % %intepolate
    % vq = griddata(abs(kp(2,:)),abs(kp(3,:)),pdf,x,y);

    clrs = flipud(colormap(hot(128)));
    cpdf = fix(127*(pdf-min(pdf))/(max(pdf)-min(pdf)))+1;

    mesh(x,y,z,ones(size(x))*10,'linewidth',0.001)  % sphere centered at origin
    hold on;
    for i=1:ndic
        plot3(real(kp(2,i)),real(kp(3,i)),real(kp(1,i)),'.',...
            'MarkerFaceColor',clrs(cpdf(i),:),'MarkerEdgeColor',clrs(cpdf(i),:));
    end
    view(0,90);xlim([0 1]);ylim([0 1]);zlim([0 1])
    axis equal tight;
    
end

