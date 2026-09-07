function [Fabs_u2o,dx2o,dy2o] = WA_BL_accurate_prop2d_v2(ui,z,dx1,lambda,q)

% 7 - ccd of capture, 2 - SLM magnified, 3 - image reconstructed from slm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geometry
% SLM planex
[Ny1,Nx1] = size(ui);
%    x1 = (-Nx1/2:Nx1/2-1)*dx1;    y1 = (-Ny1/2:Ny1/2-1)*dx1;    
dfx1 = 1/Nx1/dx1; dfy1 = 1/Ny1/dx1; 
% fx1 = (-Nx1/2:Nx1/2-1)*dfx1; fy1 = (-Ny1/2:Ny1/2-1)*dfy1;
% [X1,Y1] = meshgrid(x1,y1); 

mode_accurrate_output_fit = 0;
mode_aprox_output_fit = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generration of the object wave
Bx2 = 2*z*tan(asin(lambda/2/dx1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dermining image domain coordinates and corresponing magnification
Nx1q = Nx1*q; Ny1q = Ny1*q;   dfx1q = dfx1/q;    dfy1q = dfy1/q;
fx1e = (-Nx1q/2:Nx1q/2-1)*dfx1q; fy1e = (-Ny1q/2:Ny1q/2-1)*dfy1q; 

[FX1em,FY1em] = meshgrid(fx1e,fy1e); 
% m = sqrt(1-fx1e.^2*lambda^2).^(1/p_m);
% mmx = sqrt(1-FX1em.^2*lambda^2-FY1em.^2*lambda^2).^(1/p_m);
% mmy = sqrt(1-FX1em.^2*lambda^2-FY1em.^2*lambda^2).^(1/p_m);

[mmx,mmy,fit_hyper_coef_4_fxy] =  get_mx_my_2d(lambda,fx1e,fy1e,Nx1q,Ny1q);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % cooordinate at output
% y0 = 0;
% f0x = x0/lambda/(sqrt(z^2+x0^2+y0^2)); 
% f0y = y0/lambda/(sqrt(z^2+x0^2+y0^2)) ;
% xi = x0*z/(sqrt(z^2+x0^2+y0^2))*(1-lambda^2*f0y.^2)^(-1/4).*hypergeom([1/2, 3/4], 3/2, lambda^2*f0x.^2/(1-lambda^2*f0y.^2)) 
% yi = y0*z/(sqrt(z^2+x0^2+y0^2))*(1-lambda^2*f0x.^2)^(-1/4).*hypergeom([1/2, 3/4], 3/2, lambda^2*f0y.^2/(1-lambda^2*f0x.^2)) 
% 
% y0 = x0;
% f0x = x0/lambda/(sqrt(z^2+x0^2+y0^2)); 
% f0y = y0/lambda/(sqrt(z^2+x0^2+y0^2)) ;
% xi2 = x0*z/(sqrt(z^2+x0^2+y0^2))*(1-lambda^2*f0y.^2)^(-1/4).*hypergeom([1/2, 3/4], 3/2, lambda^2*f0x.^2/(1-lambda^2*f0y.^2)) 
% yi2 = y0*z/(sqrt(z^2+x0^2+y0^2))*(1-lambda^2*f0x.^2)^(-1/4).*hypergeom([1/2, 3/4], 3/2, lambda^2*f0y.^2/(1-lambda^2*f0x.^2)) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applying magnifications with padding
ui = padarray(ui,[(Ny1q-Ny1)/2,(Nx1q-Nx1)/2],'both');
% w = tukeywin(Ny1q,0.25)*tukeywin(Nx1*q,0.25).';
%ftuo =  fftshift(fft2(fftshift(uo.*w)));
ftuo =  fftshift(fft2(fftshift(ui)));
FX1m = FX1em.*mmx;    FY1m = FY1em.*mmy;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % to change for redusing the pixel pitch at the output
% P = 0;
% Nx1qq = Nx1+256*P;    qqx= Nx1qq/Nx1;
% Ny1qq = Ny1+128*P;    qqy= Ny1qq/Ny1;
% fx1o = (-Nx1qq/2:Nx1qq/2-1)*dfx1; dx1o = dx1/qqx;  x1o = (-Nx1qq/2:Nx1qq/2-1)*dx1o;
% fy1o = (-Ny1qq/2:Ny1qq/2-1)*dfy1; dy1o = dx1/qqy; y1o = (-Ny1qq/2:Ny1qq/2-1)*dy1o;
% [FX1o,FY1o] = meshgrid(fx1o,fy1o); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NO inetrpolation and enlarge output for signal output to fit all image
Nx1qq = Nx1;     Ny1qq = Ny1;    
dx1o = lambda*z/Bx2;  dy1o = lambda*z/Bx2;
dfx1o  = 1/Nx1qq/dx1o; dfy1o  = 1/Ny1qq/dy1o;
fx1o = (-Nx1qq/2:Nx1qq/2-1)*dfx1o; x1o = (-Nx1qq/2:Nx1qq/2-1)*dx1o;
fy1o = (-Ny1qq/2:Ny1qq/2-1)*dfy1o; y1o = (-Ny1qq/2:Ny1qq/2-1)*dy1o;
[FX1o,FY1o] = meshgrid(fx1o,fy1o); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolation part
FX1m_ = reshape(FX1m,[Ny1q*Nx1q,1]) ; FY1m_ = reshape(FY1m,[Ny1q*Nx1q,1]) ;
re_ftuo = reshape(real(ftuo),[Ny1q*Nx1q,1]) ;     imz_ftuo = reshape(imag(ftuo),[Ny1q*Nx1q,1]) ;

Fre = scatteredInterpolant(FX1m_,FY1m_,re_ftuo,'linear','none');
disp('step 1')
reftuio_ = Fre(FX1o,FY1o) ;
clear Fre ftuo
disp('step 2')
reftuio_(isnan(reftuio_)) = 0; 
Fimz = scatteredInterpolant(FX1m_,FY1m_,imz_ftuo,'linear','none');
disp('step 3')
imzftuio_ = Fimz(FX1o,FY1o) ;
disp('step 4')
clear Fimz
imzftuio_(isnan(imzftuio_)) = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fresnel processing
uo =  ifftshift(ifft2(ifftshift(reftuio_ + imzftuio_*1i)));


[X1o,Y1o] = meshgrid(x1o,y1o); 
u2o = fftshift(fft2(ifftshift(uo.*exp(1j*pi/lambda/z*((X1o.^2)+(Y1o.^2))) )));



dx2o = abs(lambda*z/dx1o/Nx1qq); x2 = (-Nx1qq/2:Nx1qq/2-1)*dx2o;
dy2o = abs(lambda*z/dy1o/Ny1qq); y2 = (-Ny1qq/2:Ny1qq/2-1)*dy2o;
[X2,Y2] = meshgrid(x2,y2);

%imagesc(x2,y2,(abs(u2o)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% accurate output fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mode_accurrate_output_fit == 1
    tic
    [sInt_X2p,sInt_Y2p] =  estimateOutInterpolants(max(1.025*Bx2/2),max(1.025*Bx2/2),z,lambda,fit_hyper_coef_4_fxy);
    X2R = sInt_X2p(X2,Y2);
    Y2R = sInt_Y2p(X2,Y2);
    plot(X2R(Ny1qq/2+1,:),abs(u2o(Ny1qq/2+1,:)))
    
    X2R_ = reshape(X2R,[Ny1qq*Nx1qq,1]) ; Y2R_ = reshape(Y2R,[Ny1qq*Nx1qq,1]) ;
    abs_u2o = reshape((abs(u2o)),[Ny1qq*Nx1qq,1]) ;
    
    Fabs_u2o = scatteredInterpolant(X2R_,Y2R_,abs_u2o,'linear','none');
    disp('step 5')
    Fabs_u2o = Fabs_u2o(X2,Y2) ;
    toc
    imagesc(x2,y2,Fabs_u2o)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% approximate output fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mode_aprox_output_fit == 1
    tic
%    X2R = X2;   Y2R = Y2;
    X2_ = reshape(X2,[Nx1qq*Ny1qq,1]) ; Y2_ = reshape(Y2,[Nx1qq*Ny1qq,1]) ;
    X2R_ = reshape(X2,[Nx1qq*Ny1qq,1]) ; Y2R_ = reshape(Y2,[Nx1qq*Ny1qq,1]) ;
%    x0 = xi2;  y0 = yi2;
    for foo = 1:4
    
    
        R2_ = sqrt(X2R_.^2+Y2R_.^2+z^2);%add
        FX2o_ = X2R_/lambda./R2_; FY2o_ = Y2R_/lambda./R2_; %add
        
        %FX2o_ = reshape(FX2o,[Nx1qq*Ny1qq,1]) ; FY2o_ = reshape(FY2o,[Nx1qq*Ny1qq,1]) ;
        argument_4_fx = lambda^2*FX2o_.^2./(1-lambda^2*FY2o_.^2);
        fit_hyper_4_fx_ = polyval(fit_hyper_coef_4_fxy,argument_4_fx);
        
        argument_4_fy = lambda^2*FY2o_.^2./(1-lambda^2*FX2o_.^2);
        fit_hyper_4_fy_ = polyval(fit_hyper_coef_4_fxy,argument_4_fy);
        
        %fit_hyper_4_fx = reshape(fit_hyper_4_fx_,[Noy,Nox]);
        %fit_hyper_4_fy = reshape(fit_hyper_4_fy_,[Noy,Nox]);
        mxo = (1-lambda^2*FY2o_.^2).^(-1/4).*fit_hyper_4_fx_;
        myo = (1-lambda^2*FX2o_.^2).^(-1/4).*fit_hyper_4_fy_;
        
%         X2R_p  = X2R_;
%         Y2R_p  = Y2R_;
        X2R_  = X2_./z.*R2_./mxo;
        Y2R_  = Y2_./z.*R2_./myo;
%         plot((X2R_p-X2R_)/dx2o)
%         mshow(reshape((X2R_p-X2R_)/dx2o,[Ny1qq,Nx1qq]))
%         contourf(reshape((X2R_p-X2R_)/dx2o,[Ny1qq,Nx1qq]),10)
    end
    abs_u2o = reshape((abs(u2o)),[Ny1qq*Nx1qq,1]) ;
    
    Fabs_u2o = scatteredInterpolant(X2R_,Y2R_,abs_u2o,'linear','none');

    disp('step 5')
    Fabs_u2o = Fabs_u2o(X2,Y2) ;
     toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %   imagesc(x2,y2,Fabs_u2o)
    
%     iix0 = round(x0/dx2o);
%     Nicut = 40;
%     img_cut = Fabs_u2o(Ny1qq/2+1-Nicut/2:Ny1qq/2+Nicut/2,Nx1qq/2+1+iix0-Nicut/2:Nx1qq/2+iix0+Nicut/2);
%     img_cut_r =  resample2(img_cut,16);
end

end

function [sInt_X2p,sInt_Y2p] =  estimateOutInterpolants(x2o_max,y2o_max,z,lambda,fit_hyper_coef_4_fxy)
Nox = 512; Noy = 256;


x2o = linspace(-x2o_max,x2o_max,Nox); y2o = linspace(-y2o_max,y2o_max,Noy);
[X2o,Y2o] = meshgrid(x2o,y2o);
R0_ = sqrt(X2o.^2+Y2o.^2+z^2);%add
FX2o = X2o/lambda./R0_; FY2o = Y2o/lambda./R0_; %add

FX2o_ = reshape(FX2o,[Nox*Noy,1]) ; FY2o_ = reshape(FY2o,[Nox*Noy,1]) ;
argument_4_fx = lambda^2*FX2o_.^2./(1-lambda^2*FY2o_.^2);
fit_hyper_4_fx_ = polyval(fit_hyper_coef_4_fxy,argument_4_fx);

argument_4_fy = lambda^2*FY2o_.^2./(1-lambda^2*FX2o_.^2);
fit_hyper_4_fy_ = polyval(fit_hyper_coef_4_fxy,argument_4_fy);

%fit_hyper_4_fx = reshape(fit_hyper_4_fx_,[Noy,Nox]);
%fit_hyper_4_fy = reshape(fit_hyper_4_fy_,[Noy,Nox]);
mxo = (1-lambda^2*FY2o_.^2).^(-1/4).*fit_hyper_4_fx_;
myo = (1-lambda^2*FX2o_.^2).^(-1/4).*fit_hyper_4_fy_;
%X2R = X0_*z./R0_.*mxo;
X2_nonuniform  = z*lambda*FX2o_.*mxo;
Y2_nonuniform  = z*lambda*FY2o_.*myo;

sInt_X2p = scatteredInterpolant(X2_nonuniform ,Y2_nonuniform ,reshape(X2o,Noy*Nox,1) ,'linear','nearest');
sInt_Y2p = scatteredInterpolant(X2_nonuniform ,Y2_nonuniform ,reshape(Y2o,Noy*Nox,1),'linear','nearest');


end