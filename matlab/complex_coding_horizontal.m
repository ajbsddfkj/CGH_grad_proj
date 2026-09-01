function [uout] = complex_coding(uc, mm_th, mm_max, Nx, dx,filtering)
% filtering == 1  for images when we have to use filter and apply modulation frequency
% filtering == 0  for points when we have band limited signal to right or left frequencies

% 
%------------ COMPLEX CODING -----------------------------
nonlinear_corection = 1;


%reuh = zeros(Nxout,Nxout);
dfx = 1/Nx/dx;
fmod = 1/4/dx;
dfx = 1/Nx/dx;
x = (-Nx/2:Nx/2-1)*dx;
[X,Y] = meshgrid(x);

%No Filter
fGmin = -1/2/dx;    fGmax =  1/2/dx-dfx;
fBmin = -1/2/dx;    fBmax =  1/2/dx-dfx;
fRmin = -1/2/dx;    fRmax =  1/2/dx-dfx;

flim(1,:)=[fBmin,fBmax]; flim(2,:)=[fGmin,fGmax];  flim(3,:)=[fRmin,fRmax];

% create frequency filter
filterX = [zeros(1,Nx/4),tukeywin(Nx/2,0.1)',zeros(1,Nx/4)];
nf_min = round(flim(1,1)/dfx)+Nx/2+1;
nf_max = round(flim(1,2)/dfx)+Nx/2+1;
filterY = [zeros(1,nf_min),tukeywin(nf_max-nf_min,0.01)',zeros(1,Nx-nf_max)];
filter  =  filterX.'*filterY;

% frequency filtering
if filtering == 1
    uho = ifft2(ifftshift(fftshift(fft2(uc)).*filter));
elseif filtering == 0
    uho = uc;
end

% amplitude correction
Betamax = 1.84; betam = -Betamax:0.001:Betamax; betamr = real(besselj(1,betam));
max_a_uhbl = max(max(abs(uho)));
if max_a_uhbl > 0
    a_uhbl = abs(uho)./max_a_uhbl*(mm_max+mm_th);
    a_uhbl(a_uhbl>mm_max) = mm_max; % application of tresholding
    a_uhbl1 = a_uhbl;
    if nonlinear_corection == 1
        a_uhbl = interp1(real(betamr),betam,a_uhbl);
    end
    
    uhblc =  a_uhbl.*exp(i*angle(uho));
    fmod_c = fmod;
    if filtering == 1
        %create symetrical signal
         uhoo = 0.5*(uhblc.*exp(1i*2*pi*Y*fmod_c) ...
                     + conj(uhblc).*exp(-1i*2*pi*Y*fmod_c));
    elseif filtering == 0
        uhoo = uhblc + conj(uhblc);
    end
%     uhoo = 0.5*(uhblc.*exp(1i*2*pi*X*fmod_c));% ...
%               %  + conj(uhblc).*exp(-1i*2*pi*X*fmod_c));
    
end
% uhoo = uho + conj(uho);

uout = exp(1j*uhoo);

end
