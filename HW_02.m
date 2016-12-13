clc
clear all
[X,x] = matches();
% b = [b11; b12; b22; b13; b23; b33;];
equations = [];
for i = 1 : length(x)
    [H hnorm invhnorm] = compute_homography(x{i}, X);
    h1 = H(:,1);
    h2 = H(:,2);
%     h3 = H(:,3);
    v12 = [h1(1)*h2(1), h1(1)*h2(2)+h1(2)*h2(1), h1(2)*h2(2), h1(3)*h2(1)+h1(1)*h2(3), h1(3)*h2(2)+h1(2)*h2(3), h1(3)*h2(3)];
    v11 = [h1(1)*h1(1), h1(1)*h1(2)+h1(2)*h1(1), h1(2)*h1(2), h1(3)*h1(1)+h1(1)*h1(3), h1(3)*h1(2)+h1(2)*h1(3), h1(3)*h1(3)];
    v22 = [h2(1)*h2(1), h2(1)*h2(2)+h2(2)*h2(1), h2(2)*h2(2), h2(3)*h2(1)+h2(1)*h2(3), h2(3)*h2(2)+h2(2)*h2(3), h2(3)*h2(3)];
    
    tmp = [v12; (v11-v22);];
    equations = [equations;tmp;];
end
B = svd(equations,0);
B11 = B(1);
B12 = B(2);
B22 = B(3);
B13 = B(4);
B23 = B(5);
B33 = B(6);
%% Intrinsic Parameters
cy = (B12*B13-B11*B23)/(B11*B22-B12*B12)
lambda = B33 - (B13*B13 + cy*(B12*B13 - B11*B23))/B11
fx = sqrt(lambda/B11)
fy = sqrt((lambda*B11)/(B11*B22-B12*B12))
cx = -B13*fx*fx/lambda
K = [fx 0 cx;0 fy cy;0 0 1]
B = [1/fx.^2 0 -cx/fx.^2;0 1/fy.^2 -cy/fy.^2;-cx/fx.^2 -cy/fy.^2 cx/fx.^2+cy/fy.^2+1;];
%% Extrinsic Parameters
files = dir('images/*.jpg');
for i = 1 : length(x)
    [H hnorm invhnorm] = compute_homography(x{i}, X);
    h1 = H(:,1);
    h2 = H(:,2);
    h3 = H(:,3);
    r1 = lambda*inv(K)*h1;
    r2 = lambda*inv(K)*h2;
    t = lambda*inv(K)*h3;
%     r3 = r1*r2;
    R(i).data = [r1 r2];
    disp(['R(',num2str(i),'):'])
    disp(R(i).data);
    T(i).data = t;
    disp(['T(',num2str(i),'):'])
    disp(T(i).data);
    points = x{i};
    err = K*[r1 r2 t]*[X; zeros(1,size(X,2))];
    nx = [points; zeros(1,size(points,2))] - err;
    img = imread(strcat('images/',files(i).name));
    figure
    imshow(img);
    hold on
    plot(nx(1,:),nx(2,:),'ro','MarkerSize',5);
    hold off
end