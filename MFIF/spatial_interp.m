function out = spatial_interp(in, warp, str, transform, nx, ny)

if (strcmp(transform,'affine')||strcmp(transform,'euclidean'))
   if size(warp,1)==2
       warp=[warp;zeros(1,3)];
   end
end

if strcmp(transform,'translation')
    warp = [eye(2) warp];
    warp = [warp; zeros(1,3)];
end

[xx yy] = meshgrid(nx, ny);
xy=[xx(:)';yy(:)';ones(1,length(yy(:)))];


A = warp;
A(3,3) = 1;

% new coordinates
xy_prime = A * xy;

if strcmp(transform,'homography')

    % division due to homogeneous coordinates
    xy_prime(1,:) = xy_prime(1,:)./xy_prime(3,:);
    xy_prime(2,:) = xy_prime(2,:)./xy_prime(3,:);
end

% Ignore third row
xy_prime = xy_prime(1:2,:);

% Subpixel interpolation
out = interp2(in, xy_prime(1,:), xy_prime(2,:), str);
out(isnan(out))=0;%replace Nan
out=reshape(out,length(ny),length(nx));
