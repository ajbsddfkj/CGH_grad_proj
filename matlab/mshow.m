function mshow(p1,p2,p3,p4)
 
if nargin==1
    minv = min(min(p1)); maxv = max(max(p1));
    if minv ~= maxv
        imagesc(p1,[minv,maxv])
        disp(['minv =  ',num2str(minv),' maxv = ',num2str(maxv)]);
    else
        disp(['minv = maxv = ',num2str(minv)]);
        disp('image disply error');
    end
elseif nargin==3  
    %m = p3;
    imagesc(p1,p2,p3)
    minv = min(min(p3)); maxv = max(max(p3));
    disp(['minv =  ',num2str(minv),' maxv = ',num2str(maxv)]);

% colormap(gray)
% colormap default
elseif nargin==4    
    %m = p3;
    imagesc(p1,p2,p3,[min(min(p3)),max(max(p3)/p4)])
    minv = min(min(p3)); maxv = max(max(p3));
    disp(['minv =  ',num2str(minv),' maxv = ',num2str(maxv)]);
%    imagesc(p2,p3,m )
% colormap(gray)
%colormap default
elseif nargin==2
    
    p1(p1 > max(max(p1))/p2) = max(max(p1))/p2;
    p1 = imresize(p1,[1024*2.0,1024*3], 'bilinear');
    imagesc(p1)
    minv = min(min(p1)); maxv = max(max(p1));
    disp(['minv =  ',num2str(minv),' maxv = ',num2str(maxv)]);
else 
    error('ktzerror: wrong params!')
end



impixelinfo

%colormap gray