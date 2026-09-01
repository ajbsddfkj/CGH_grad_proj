function [h] = wa_cgh_SpectrumRecordingPlaneMemory1_v2(Nx1,Ny1,dx1,lambda,xo,yo,zo,ao)

    k = 2*pi/lambda;    %lambda_m2 = lambda^-2; lambda_2 = lambda^2;
    dfx1 = 1/Nx1/dx1;    dfy1 = 1/Ny1/dx1; 
    fx1 = (-Nx1/2:Nx1/2-1)*dfx1; fy1 = (-Ny1/2:Ny1/2-1)*dfy1;
    ft_hre = zeros(Ny1,Nx1);
    ft_himz = zeros(Ny1,Nx1);
    [FX1,FY1] = meshgrid(fx1,fy1);
    fi_FZ_noZ = k*sqrt(1-lambda^2*FX1.^2-lambda^2*FY1.^2);
    FX1_2pi = FX1*2*pi;
    FY1_2pi = FY1*2*pi;
%    IIY = 1:Ny1;  IIX = 1:Nx1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uncomment for smaller memory usage
%    clear 'FX1' 'FY1'
%    Nx1s = 48;  Ny1s = 24; 
%    fx1s = (-Nx1s/2:Nx1s/2-1)*dfx1; fy1s = (-Ny1s/2:Ny1s/2-1)*dfy1;
%    [FX1s,FY1s] = meshgrid(fx1s,fy1s);
%     FX1_2piS = FX1s*2*pi;
%     FY1_2piS = FY1s*2*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     fi_h = zeros(Ny1s,Nx1s);
    
    Bx1_2_o_lam = Nx1/2*dx1/lambda;  By1_2_o_lam = Ny1/2*dx1/lambda;
    Nxo2p1 = Nx1/2+1;   Nyo2p1 = Ny1/2+1;
    
    Npo = length(xo);
    
%     IIY = 1:Nx1*Ny1;
%     IIY = reshape(IIY,[Ny1,Nx1]);
abs_ao = abs(ao);
ph_ao = angle(ao);
tic
    for foo = 1 : Npo
    
%         yo(foo) = Bx2 * 0.0001;
%         xo(foo) = Bx2 * 0.35;
%         zo(foo)  = z;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calc. f0xy
        zo2 = zo(foo)^2;    xo2 = xo(foo)^2;    yo2 = yo(foo)^2;
        %R02 = (xo2+yo2+zo2);
        R0 = sqrt(xo2+yo2+zo2);
%         rl0x = -R0*(R02-yo2)/(zo2);
%         rl0y = -R0*(R02-xo2)/(zo2);
        rl0x = -R0*(zo2+xo2)/(zo2);
        rl0y = -R0*(zo2+yo2)/(zo2);
%         rl0x = -R0*(1+xo2/(zo2));
%         rl0y = -R0*(1+yo2/(zo2));
        f0x = (xo(foo))/R0/lambda;
        f0y = (yo(foo))/R0/lambda;

%         R02 = (xo2+yo2+zo2);
%         
%         R0 = sqrt(R02);
%         rl0x2 =  (1-yo2/R02)*zo(foo)/sqrt(1-yo2/R02-xo2/R02)^3;
%         rl0x2 = R0.*(R02-yo2)*zo(foo)/sqrt(R02-yo2-xo2)^3
%         rl0x2 = R0.*(R02-yo2)*zo(foo)/sqrt(zo2)^3
        
% 
%         R0_lam_2 = (xo2+yo2+zo2)*lambda_2;
%         f0x2 = xo2/R0_lam_2;
%         f0y2 = yo2/R0_lam_2;
%         f0x = sign(xo(foo))*sqrt(f0x2);
%         f0y = sign(yo(foo))*sqrt(f0y2);
%         K = sqrt(1-lambda^2*f0x2-lambda^2*f0y2).^3;
%         rl0x =  -(1-lambda^2*f0y2)*zo(foo)/K;
%         rl0y =  -(1-lambda^2*f0x2)*zo(foo)/K;

    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fq. limts for x
        %x0_plus = x0+Bx1_2  x0_minus = x0-Bx1_2
%         xo_plus =   xo(foo)+Bx1_2;  xo_minus = xo(foo)-Bx1_2;
%         xo_plus2  = (xo_plus)^2;    xo_minus2 = (xo_minus)^2;
%     
%         %fxlim_minus = sign(x0_minus)*sqrt((x0_minus^2*lambda^-2-x0_minus^2*f0y^2)/(x0_minus^2+z^2));
%         fxlim_minus = sign(xo_minus)*sqrt((xo_minus2*lambda_m2-xo_minus2*f0y2)/(xo_minus2+zo2));
%         %fxlim_plus = sign(x0_plus)*sqrt((x0_plus^2*lambda^-2-x0_plus^2*f0y^2)/(x0_plus^2+z^2));
%         fxlim_plus = sign(xo_plus)*sqrt((xo_plus2*lambda_m2-xo_plus2*f0y2)/(xo_plus2+zo2));
        fxlim_minus = +Bx1_2_o_lam/rl0x+f0x;
        fxlim_plus = -Bx1_2_o_lam/rl0x+f0x;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fq. limts for y
%         % y0_plus  = y0 + By1_2;
%         % y0_minus = y0 - By1_2;
%         yo_plus =   yo(foo)+By1_2;  yo_minus = yo(foo)-By1_2;
%         yo_plus2  = (yo_plus)^2;    yo_minus2 = (yo_minus)^2;
%     
%         %fylim_minus = sign(y0_minus)*sqrt((y0_minus^2*lambda^-2-y0_minus^2*f0x^2)/(y0_minus^2+z^2));
%         fylim_minus = sign(yo_minus)*sqrt((yo_minus2*lambda_m2-yo_minus2*f0x2)/(yo_minus2+zo2));
%         %fylim_plus = sign(y0_plus)*sqrt((y0_plus^2*lambda^-2-y0_plus^2*f0x^2)/(y0_plus^2+z^2));
%         fylim_plus  = sign(yo_plus)*sqrt((yo_plus2*lambda_m2-yo_plus2*f0x2)/(yo_plus2+zo2));
        fylim_minus = +By1_2_o_lam/rl0y+f0y;
        fylim_plus = -By1_2_o_lam/rl0y+f0y;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rounding fq. limts for xy
%          ii_fxlim_minus  = max(Nxo2p1+round(fxlim_minus/dfx1),1);
%          ii_fxlim_plus   = min(Nxo2p1+round(fxlim_plus/dfx1),Nx1);
%          ii_fxlim_minus  = Nxo2p1+ceil(fxlim_minus/dfx1);
%          ii_fxlim_plus   = Nxo2p1+floor(fxlim_plus/dfx1);
         ii_fxlim_minus  = Nxo2p1+round(fxlim_minus/dfx1);
         ii_fxlim_plus   = Nxo2p1+round(fxlim_plus/dfx1);
    
%        ii_fylim_minus = max(Nyo2p1+round(fylim_minus/dfy1),1);
%        ii_fylim_plus =  min(Nyo2p1+round(fylim_plus/dfy1),Ny1);
%        ii_fylim_minus = Nyo2p1+ceil(fylim_minus/dfy1);
%        ii_fylim_plus = Nyo2p1+floor(fylim_plus/dfy1);
         ii_fylim_minus = Nyo2p1+round(fylim_minus/dfy1);
         ii_fylim_plus = Nyo2p1+round(fylim_plus/dfy1);
    
         iiy = ii_fylim_minus:ii_fylim_plus;
         iix = ii_fxlim_minus:ii_fxlim_plus;
         
         if ii_fylim_minus < 1 ; continue; end
         if ii_fxlim_minus < 1 ; continue; end
         if ii_fylim_plus > Ny1 ; continue; end
         if ii_fxlim_plus > Nx1 ; continue; end 
         
%          iiyS = 1:length(iiy);
%          iixS = 1:length(iix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uncomment for smaller memory usage
%         iiyS = 1:ii_fylim_plus-ii_fylim_minus+1;
%         iixS = 1:ii_fxlim_plus-ii_fxlim_minus+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        %   fi_h = (fi_FZ(iiy,iix)+2*pi*FX1(iiy,iix)*x0+2*pi*FY1(iiy,iix)*y0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uncomment for smaller memory usage
%         fi_h = (fi_FZ_noZ(iiy,iix)*zo(foo)+FX1_2piS(iiyS,iixS)*xo(foo)+FY1_2piS(iiyS,iixS)*yo(foo)+ph_ao(foo));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comment for smaller memory usage
         fi_h = (fi_FZ_noZ(iiy,iix)*zo(foo)+FX1_2pi(iiy,iix)*xo(foo)+FY1_2pi(iiy,iix)*yo(foo)+ph_ao(foo));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ft_hre(iiy,iix)  = ft_hre(iiy,iix) + cos(-fi_h)*abs_ao(foo);
         ft_himz(iiy,iix) = ft_himz(iiy,iix) + sin(-fi_h)*abs_ao(foo);
    end

     h = fftshift(ifft2(fftshift(ft_hre+1i*ft_himz)));
toc
end