function M=rotate(N,tetax,tetay,tetaz)

rotx=[1 0 0; 0 cosd(tetax) -sind(tetax);0 sind(tetax) cosd(tetax)];
roty=[cosd(tetay) 0 sind(tetay);0 1 0; -sind(tetay) 0 cosd(tetay)];
rotz=[cosd(tetaz) -sind(tetaz) 0; sind(tetaz) cosd(tetaz) 0; 0 0 1];


if size(N,2)==9||size(N,2)==10
M1=N(:,1:3);
M2=N(:,4:6);
M3=N(:,6:9);
M1obr=M1*roty*rotx*rotz;
M2obr=M2*roty*rotx*rotz;
M=cat(2,M1obr,M2obr,M3);

elseif size(N,2)==6
    
M1=N(:,1:3);
M2=N(:,4:6);
M1obr=M1*roty*rotx*rotz;

M=cat(2,M1obr,M2);

elseif size(N,2)==3
    
    M1=N(:,1:3);
    M=M1*roty*rotx*rotz;
    
end
