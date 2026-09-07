function [mx,my,fit_hyper_coef_4_fxy] =  get_mx_my_2d(lambda,fx1,fy1,Nx1,Ny1)
% input parameters

% % finding model with fitting
% z = .5e6;
% lambda = 0.633;
% dx1 = 0.5;Nx1 = 2048*4/1; Ny1 = 2048*2/1;

%fx_max  = 1/2/dx1; fy_max  = 1/2/dx1;
fx_max  = max(fx1); fy_max  = max(fy1);
max_lfxy_4_aprox =  lambda^2*fx_max.^2/(1-lambda^2*fy_max.^2);
lam_fxy = linspace(-max_lfxy_4_aprox,max_lfxy_4_aprox,51);
hyper_coef = hypergeom([1/2, 3/4], 3/2, lam_fxy);
fit_hyper_coef_4_fxy = polyfit(lam_fxy,hyper_coef,9); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation of x y magnifications
% dfx1 = 1/Nx1/dx1; dfy1 = 1/Ny1/dx1; 
% fx1 = (-Nx1/2:Nx1/2-1)*dfx1+dfx1/2; fy1 = (-Ny1/2:Ny1/2-1)*dfy1+dfx1/2; 
[FX1,FY1] = meshgrid(fx1,fy1);

FX1_ = reshape(FX1,[Ny1*Nx1,1]) ; FY1_ = reshape(FY1,[Ny1*Nx1,1]) ;
argument_4_fx = lambda^2*FX1_.^2./(1-lambda^2*FY1_.^2);
fit_hyper_4_fx_ = polyval(fit_hyper_coef_4_fxy,argument_4_fx);

argument_4_fy = lambda^2*FY1_.^2./(1-lambda^2*FX1_.^2);
fit_hyper_4_fy_ = polyval(fit_hyper_coef_4_fxy,argument_4_fy);

fit_hyper_4_fx = reshape(fit_hyper_4_fx_,[Ny1,Nx1]);
fit_hyper_4_fy = reshape(fit_hyper_4_fy_,[Ny1,Nx1]);
mx = (1-lambda^2*FY1.^2).^(-1/4).*fit_hyper_4_fx;
my = (1-lambda^2*FX1.^2).^(-1/4).*fit_hyper_4_fy;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % testing magnification
%  FX1_m = FX1.*mx;
% [d_fx1_m_0,~]=gradient(FX1_m,dfx1);
% r_l_fx = (1-lambda^2*FY1.^2)*z./sqrt(1 - lambda^2*FX1.^2-lambda^2*FY1.^2).^3;
% rl_m_x = (r_l_fx).*d_fx1_m_0.^-2;
% mshow(rl_m_x)
% 
% FY1_m = FY1.*my;
% [~,d_fy1_m_0]=gradient(FY1_m,dfx1,dfy1);
% r_l_fy = (1-lambda^2*FX1.^2)*z./sqrt(1 - lambda^2*FX1.^2-lambda^2*FY1.^2).^3;
% rl_m_y = (r_l_fy).*(d_fy1_m_0.^-2);
% mshow(rl_m_y)
% sdfs

