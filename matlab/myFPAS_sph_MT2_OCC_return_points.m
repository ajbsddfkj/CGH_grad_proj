function [ outidx ] = myFPAS_sph_MT2_OCC_return_points( xop,yop,zop,uop,Nx,Ny,dx,lambda,Nxs,q,zr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FPAS_CGH_3DPS_2d_dev
% method for development of PAS computer generated hologram
% Implementation of paper: Hoonjong Kang, Takeshi Yamaguchi, and Hiroshi
% Yoshikawa, "Accurate phase-added stereogram to improve the coherent
% stereogram," Appl. Opt. 47, D44-D54 (2008)
% % Kang H, Fujii T, Yamaguchi T, Yoshikawa H; Compensated phase-added
% with elemnt from chapter on holographic printing
% Advance over CPAS_CGH_3DPS_1d_dev: more accurate sampling of plane waves
% from 3D point source data
% 2D method
% !!!!   NOTE  !!!!:
% 1. It works for squered segments of size Nxs x Nxs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Zaux,Iz] = sort(zop);
[Xaux]=xop(Iz);
[Yaux]=yop(Iz);
Caux = uop(Iz);

xo=Xaux;
yo=Yaux;
zo=Zaux;
uo=Caux/(max(Caux));

Np = length(xo); k0 = 2*pi/lambda;
if(Np ~= mean(length(yo),length(zo)))
    error('ERR: wrong input data')
end
x = (-Nx/2 :Nx/2-1)*dx;  y = (-Ny/2:Ny/2-1)*dx;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%definition of double segment
Nxs2 = q*Nxs;
Nosx = Nx/Nxs;  Nosy = Ny/Nxs; % Noumber of segments in X and Y directions


Ts = Nxs*dx; Ts2 = Nxs2*dx; dfxs2=1/Ts2;
fxs2 =  (-Nxs2/2:Nxs2/2-1)*dfxs2; % frequency spectrum for double segment

x0seg = x(1)+(0:Nosx-1)*Ts+Ts/2;% location of segments center for X
y0seg = y(1)+(0:Nosy-1)*Ts+Ts/2;% location of segments center for Y

nseg_bx = 1+(0:Nosx-1)*Nxs;% location of segments startst for X
nseg_by = 1+(0:Nosy-1)*Nxs;% location of segments startst for X



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hologram computation
ppool = gcp;
worknum = ppool.NumWorkers;
% seg_per_work = ceil(length(SegIdx)/worknum);
divid = zeros(worknum,1);




pcell = cell(Nosy, Nosx);

zo2 = zo.^2;
% SegIdx = [SegIdx; zeros(worknum*seg_per_work-length(SegIdx),2)];
jj=0; % 0 for the clasical Fresnel generation of the complex field OR* (juan addition);
% 1 for the Fourier generation of the complex field OR* (juan addition);
parfor foosy = 1 : Nosy %%%%%Oryginalnie najpierf for foosx potem for foosy
    for foosx = 1: Nosx
        fths = zeros(q*Nxs,q*Nxs);
        rrs = zeros(Nxs2,Nxs2);
        pidx = zeros(Nxs2,Nxs2);
        yr = y0seg(foosy);
        xr = x0seg(foosx);
        rr = sqrt(zr.^2 + xr.^2 + yr.^2);
        
        fxr = 1*xr/rr/lambda;
        fyr = 1*yr/rr/lambda;
        cc = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Points action
        for foo = 1 : Np
            
            yp = yo(foo) - yr;
            xp = xo(foo) - xr;
            rp = sqrt(zo2(foo) + xp.^2 + yp.^2);
            
            fxp = xp/rp/lambda;
            fyp = yp/rp/lambda;
            
            fxp = fxp + jj*fxr;
            fyp = fyp + jj*fyr;
            
            iifx = round(fxp/dfxs2)+Nxs2/2+1;
            iify = round(fyp/dfxs2)+Nxs2/2+1;
            
            iixlogic = (iifx<=0 | iifx>Nxs2);
            iiylogic = (iify<=0 | iify>Nxs2);
            
            if iixlogic || iiylogic
                iifx = Nxs2/2+1;
                iify = Nxs2/2+1;
            end
            
            C_eta=(fxp-fxs2(iifx));
            C_ecsi=(fyp-fxs2(iify));
            rr0 =(rp-jj*rr);
            [rrs,pidx]=selector_s_occ(foo,rrs,rr0,iifx,iify,cc,pidx);
            
        end
        %figure, nimage(pidx)
        %         temp = pidx(:);
        %         temp = temp(temp~=0);
        %         id_num = [id_num; temp];
        %         id_num = unique(id_num);
        pidx = pidx(:);
        
        pcell{foosy,foosx} = pidx(pidx~=0);
        disp(foosx + (foosy-1)*Nosx)
    end
end

% h = zeros(Ny,Nx);
% for fooh = 1 : worknum
%     for foopw = 1: seg_per_work
%         foos = (fooh-1)*seg_per_work+foopw;
%         if foos > length(SegIdx)
%             break
%         end
%         h(nseg_byy(foos):nseg_byy(foos)+Nxs-1, nseg_bxx(foos) : nseg_bxx(foos)+Nxs-1 ) = hcell{1,fooh}((foopw-1)*Nxs+1:foopw*Nxs,:);
%     end
% end
outidx = [];
for foosy = 1: Nosy
    for foosx = 1:Nosx
        outidx = [outidx; pcell{foosy,foosx}];
        outidx = unique(outidx);
    end
end
outidx = sort(unique(outidx));

