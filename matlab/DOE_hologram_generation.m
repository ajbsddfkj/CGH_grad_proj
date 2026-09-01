clear all

% input parameters
% wavelength
lambda = 0.520;

% reconstruction distance
z = 1e6; z0 = z;
% 4F system parameters (near-eye display)
F1 = 100e3; F2 = 35e3;  M4F = F2/F1;

% SLM parameters
dxs = 3.74; Bfs = 1/dxs;
Nx1 = 2048*1; dx1 = 3.74*M4F;  Ny1 = Nx1/1;
Nx1 = 4160; Ny1 = 2464; Nxslm = 4160; Nyslm = 2464;
Nxo = 2048+0*1024; Nyo = Nxo/2;

% create grids in space and frequency domains
x1 = (-Nx1/2:Nx1/2-1)*dx1;    dfx1 = 1/Nx1/dx1; dfy1 = 1/Ny1/dx1; y1 = (-Ny1/2:Ny1/2-1)*dx1;
fx1 = (-Nx1/2:Nx1/2-1)*dfx1;    fy1 = (-Ny1/2:Ny1/2-1)*dfy1;
[FX1, FY1] = meshgrid(fx1,fy1);
[X1, Y1] = meshgrid(x1,y1);
[Xs, Ys] = meshgrid((-Nx1/2:Nx1/2-1)*dxs,(-Ny1/2:Ny1/2-1)*dxs);

% set FOV
FoV = asin(lambda/2/dx1);
z = 1e6;
Bx = 2*z*tan(FoV)*0.99;
By = 2*z*tan(FoV)*0.99;
ymax = By/2;

% occulsion parametes
wielk_siatki=3;
kernel=1;
start = 5.03;

% read cloud
A=(pcread('Boat_4M.ply'));
% A=(pcread('Boat_20M.ply'));

% show cloud
pcshow(A.Location,A.Color);

% prepare position of the cloud
x_width= max(A.Location(:,1))-min(A.Location(:,1));
x_Off=0.8*x_width;
Loc_shifted=[A.Location(:,1)-x_Off,A.Location(:,2),A.Location(:,3)];
Loc_rot=rotate(Loc_shifted,13,start,180);
Loc_rot(:,1) = Loc_rot(:,1) - mean(Loc_rot(:,1));
pcshow(Loc_rot,A.Color);

% take cloud parameters
xo = single(Loc_rot(:,1)); yo = single(Loc_rot(:,2)); zo = single(Loc_rot(:,3));
c = single([A.Color(:,1),A.Color(:,2),A.Color(:,3)]);

xomax = (max(max(abs(xo))));
yomax = (max(max(abs(yo))));

% cloud scaling
offset_xy=0;
lambda_def = 0.52;
scale_ratio = asind(lambda/2/dx1)/asind(lambda_def/2/dx1);
scale = 8.089992819529476e+03 * scale_ratio/4;
xo = (xo-0.14*xomax).*scale + offset_xy; yo = (yo).*scale+ymax*0.05+offset_xy; zo = zo.*scale+z;

pointsx = sortrows([xo,yo,zo,c],3);
xoo = pointsx(:,1);
yoo = pointsx(:,2);
zoo = pointsx(:,3);
aoo = pointsx(:,4);

% calculate occlusions - delete hidden points
ii = myFPAS_sph_MT2_OCC_return_points(xoo,yoo,zoo,aoo,Nxo,Nyo,dx1,lambda,128*2*kernel,wielk_siatki,z);
xoo = (pointsx(ii,1)); yoo = (pointsx(ii,2)); zoo = (pointsx(ii,3));
aoo = ((pointsx(ii,4)));

ci=1;
% calculate hologram
h = wa_cgh_SpectrumRecordingPlaneMemory1_v2(Nx1,Ny1,dx1,lambda,xoo,yoo,zoo,aoo);
figure
imagesc(abs(fftshift(fft2(h))))
max_z = max(max(zoo))
min_z = min(min(zoo))


% preparing for SLM
Nxslm = 4160; Nyslm = 2464;
uslm = zeros(Nyslm,Nxslm);
uslm(:,(Nxslm-Nx1)/2+1:(Nxslm+Nx1)/2) = h((Ny1-Nyslm)/2+1:(Ny1+Nyslm)/2,:);
        
uslm = angle(uslm);
uslm = mat2gray(imcomplement(uslm));

% save hologram
folder = 'Boat_DOE_rotate_off_axis';
if not(isfolder(folder))
    mkdir(folder)
end
imwrite(uslm, strcat(folder,'\Boat_4M_DOE_',num2str(rot),'.bmp'));

% reconstruct hologram
q = 1.125;  % pading of input

fooz = 1;
[uo,dx2,dy2] = WA_BL_accurate_prop2d_v2(h,z,dx1,lambda,q);
%[uo,dx2,dy2] = WA_BL_accurate_prop2d_PartialCenterCenter_v2(h,z,dx1,lambda,q);
figure(2)
uo = rot90(uo,2);

% show reconstruction
mshow(abs(uo).^1,2);