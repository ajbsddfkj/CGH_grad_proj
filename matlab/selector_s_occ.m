function [rrs, pidx]=selector_s_occ(foo,rrs,rr0,iifx,iify,cc,pidx)

%             if ~any(rrs(iify,iifx))
            if rrs(iify,iifx)==0

                    rrs(iify,iifx) =  rr0;
                   % fths(iify,iifx) =  uo(foo).*exp(-1i*(k0*rr0 + 2*pi*((fxs2(iifx)+fxs2(iify)))*(Ts/2)+2*pi*((C_eta+C_ecsi)))) ;
%                     cloud(cc,1)=xo(foo);
%                     cloud(cc,2)=yo(foo);
%                     cloud(cc,3)=zo(foo);
%                     cloud(cc,4)=uo(foo);
                    pidx(iify,iifx) = foo;% - only for surviving points
                   % extraction
                    cc=cc+1;
                    
                
            elseif rrs(iify,iifx) > rr0                 
                   rrs(iify,iifx) =  rr0;
 %                   fths(iify,iifx) = uo(foo).*exp(-1i*(k0*rrs(iify,iifx) + 2*pi*((fxs2(iifx)+fxs2(iify)))*(Ts/2)+2*pi*((C_eta+C_ecsi)))) ;
%                     cloud(cc,1)=xo(foo);
%                     cloud(cc,2)=yo(foo);
%                     cloud(cc,3)=zo(foo);
%                     cloud(cc,4)=uo(foo);
                     pidx(iify,iifx) = foo;
                    cc=cc+1;
%                 end    
            end  