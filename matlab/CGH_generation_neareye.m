clear

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
Bx = 2*z*tan(FoV)*0.8;
By = 2*z*tan(FoV)*0.8;
ymax = By/2;

% occulsion parametes
wielk_siatki=3;
kernel=1;
start = 5.03;

nazwa_folderu = {'zmienNazwe11','zmienNazwe12','zmienNazwe13','zmienNazwe21','zmienNazwe22','zmienNazwe23','zmienNazwe31','zmienNazwe32','zmienNazwe33'};

nazwy = {'multiplex1.bmp', 'multiplex2.bmp','multiplex3.bmp','multiplex4.bmp','multiplex5.bmp','multiplex6.bmp','multiplex7.bmp','multiplex8.bmp','multiplex9.bmp'};


rot = [0  0  0;  % Vertex 1
       4  0  0;  % Vertex 2
       0  4  0;  % Vertex 3
       4  4  0;  % Vertex 4
       -4  0  0;  % Vertex 5
       0  -4  0;  % Vertex 6
       -4  -4  0;  % Vertex 7
       -4  4  0;  % Vertex 8
       4  -4  0];  % Vertex 9

count = 0;

%% zacznij obliczenia

for i = 1:9
    fprintf('Current iteration: %d\n', i);
    % Sprawdź czy folder istnieje. Jeśli nie (warunek tyldy '~'), utwórz go.
    if ~exist(nazwa_folderu{i}, 'dir')
        mkdir(nazwa_folderu{i});
    end
    
    % read and preprocess cloud - ZMIEN NAZWE
    A=(pcread('zmienNazwe_gpu_z_amplituda.ply'));
    
    % show cloud
    %pcshow(A.Location,A.Color);
    
    % show cloud no color
    %pcshow(A.Location);
    
    % 3. Przygotowanie wizualizacji
    % figure('Name', 'Weryfikacja Obrotu 3D', 'Color', 'w');
    % hold on; 
    % grid on;
    
    % X_model = A.Location(:,1);
    % Y_model = A.Location(:,2);
    % Z_model = A.Location(:,3);
    % scatter3(X_model, Y_model, Z_model, 2, 'r', 'filled');
    
    
    
    % prepare position of the cloud
    x_width= max(A.Location(:,1))-min(A.Location(:,1));
    x_Off=0.8*x_width;
    %Loc_shifted=[A.Location(:,1)-x_Off,A.Location(:,2),A.Location(:,3)];
    
    % location without shift
    Loc_shifted=[A.Location(:,1),A.Location(:,2),A.Location(:,3)];
    
    xyz = rot(i,:)

    x_rot = xyz(1);
    y_rot = xyz(2);
    z_rot = xyz(3);
    
    Loc_rot=rotate(Loc_shifted,x_rot,y_rot,z_rot);
    
    %Loc_rot=Loc_shifted;
    Loc_rot(:,1) = Loc_rot(:,1) - mean(Loc_rot(:,1));
    Loc_rot(:,2) = Loc_rot(:,2) - mean(Loc_rot(:,2));
    Loc_rot(:,3) = Loc_rot(:,3) - mean(Loc_rot(:,3));
    
    
    % oś z - to obrót wokół azymutu - czyli jak gdyby obrót patrząc do dołu -
    % trzeba obracać układ w osiach x i y
    
    % -- Rysowanie ORYGINALNEGO modelu (niebieski) --
    % X_orig = Loc_rot(:, 1);
    % Y_orig = Loc_rot(:, 2);
    % Z_orig = Loc_rot(:, 3);
    % % Używamy scatter3: (X, Y, Z, rozmiar_punktu, kolor, 'filled' do wypełnienia)
    % scatter3(X_orig, Y_orig, Z_orig, 2, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
    
    % % 4. Ustawienia wykresu (kluczowe dla oceny obrotu)
    % axis equal; % BARDZO WAŻNE: Wymusza równe proporcje osi, aby obrót nie wyglądał jak spłaszczenie
    % xlabel('Oś X');
    % ylabel('Oś Y');
    % zlabel('Oś Z');
    % legend('Oryginał', 'Po obrocie', 'Location', 'best');
    % 
    % % Narysowanie osi układu współrzędnych w punkcie (0,0,0) dla lepszej orientacji
    % plot3([0 15], [0 0], [0 0], 'k-', 'LineWidth', 1.5); % Oś X
    % plot3([0 0], [0 15], [0 0], 'k-', 'LineWidth', 1.5); % Oś Y
    % plot3([0 0], [0 0], [0 15], 'k-', 'LineWidth', 1.5); % Oś Z
    % 
    % view(3); % Ustawia domyślny, czytelny widok 3D
    % hold off;
    
    %pcshow(Loc_rot,A.Color);
    
    % show cloud no color
    pcshow(Loc_rot);
    
    % Ustawienie kamery: view(azymut, elewacja) w stopniach
    view(45, -90);
    % Zapisz plik bezpośrednio do tego folderu
    exportgraphics(gcf, fullfile(nazwa_folderu{i}, 'chmura.png'), 'Resolution', 300, 'BackgroundColor', 'current');
    
    % take cloud parameters
    xo = single(Loc_rot(:,1)); yo = single(Loc_rot(:,2)); zo = single(Loc_rot(:,3));
    c = single([A.Color(:,1),A.Color(:,2),A.Color(:,3)]);
    
    % no color
    %c = ones(length(xo), 1, 'single');
    
    xomax = (max(max(abs(xo))));
    yomax = (max(max(abs(yo))));
    
    % cloud scaling
    offset_xy=0;
    lambda_def = 0.52;
    scale_ratio = asind(lambda/2/dx1)/asind(lambda_def/2/dx1);
    scale = 0.489992819529476e+03 * scale_ratio/4;
    %xo = (xo-0.14*xomax).*scale + offset_xy; yo = (yo).*scale+ymax*0.05+offset_xy; zo = (zo.*scale) +z;
    
    xo = (xo-0.14*xomax).*scale + offset_xy; yo = (yo).*scale+ymax*0.05+offset_xy; zo = zo - mean(zo); zo = zo .* scale; zo = zo - min(zo); zo = zo + 1e6;
    
    % new cloud scaling
    %xo = xo .* scale + offset_xy; yo = yo .* scale + offset_xy; zo = zo .* scale + z; % Teraz środek ciężkości osi Z wyląduje idealnie w odległości z (1 metr)
    
    pointsx = sortrows([xo,yo,zo,c],3);
    xoo = pointsx(:,1);
    yoo = pointsx(:,2);
    zoo = pointsx(:,3);
    aoo = pointsx(:,4);
    
    % calculate occlusions - delete hidden points
    %ii = myFPAS_sph_MT2_OCC_return_points(xoo,yoo,zoo,aoo,Nxo,Nyo,dx1,lambda,128*2*kernel,wielk_siatki,z);
    %xoo = (pointsx(ii,1)); yoo = (pointsx(ii,2)); zoo = (pointsx(ii,3));
    %aoo = ((pointsx(ii,4)));
    
    ci=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preview processed point cloud
    % Wyświetlenie chmury po samym przeskalowaniu
    figure(10); % Otwieramy nowe okno, żeby nie nadpisać poprzednich wykresów
    pcshow([xo, yo, zo], c); 
    title('Chmura punktów po skalowaniu (po okluzji)');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    
    % Zapisz plik bezpośrednio do tego folderu
    exportgraphics(gcf, fullfile(nazwa_folderu{i}, 'chmura-skal.png'), 'Resolution', 300, 'BackgroundColor', 'current');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % calculate hologram
    h = wa_cgh_SpectrumRecordingPlaneMemory1_v2(Nx1,Ny1,dx1,lambda,xoo,yoo,zoo,aoo);
    figure
    imagesc(abs(fftshift(fft2(h))))
    % Zapisz plik bezpośrednio do tego folderu
    %exportgraphics(gcf, fullfile(nazwa_folderu, 'rekonstrukcja.png'));
    
    max_z = max(max(zoo))
    min_z = min(min(zoo))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % hologram saving
    % preparing for SLM
    Nxslm = 4160; Nyslm = 2464;
    uslm = zeros(Nyslm,Nxslm);
    uslm(:,(Nxslm-Nx1)/2+1:(Nxslm+Nx1)/2) = h((Ny1-Nyslm)/2+1:(Ny1+Nyslm)/2,:);
            
    uslm = angle(uslm);
    uslm = mat2gray(imcomplement(uslm));
    
    % save hologram
    imwrite(uslm, fullfile(nazwa_folderu{i},'hologram.bmp'));
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %% import hologram NIECZYNNE
    % 
    % % read the image and convert to 0.0 to 1.0 range
    % img = im2double(imread(fullfile(nazwa_folderu,'hologram.bmp')));
    % 
    % % scale and inversion back to radians (-pi to pi)
    % phase_slm = pi - (img *2 * pi);
    % 
    % % removing the padding
    % phase_cropped = phase_slm(:, (Nxslm-Nx1)/2+1 : (Nxslm+Nx1)/2);
    % 
    % % recreating the complex wavefront (remove later)
    % h = exp(1i * phase_cropped);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reconstruct hologram
    
    q = 1.125;  % pading of input
    fooz = 1;
    [uo,dx2,dy2] = WA_BL_accurate_prop2d_v2(h,z,dx1,lambda,q);
    %[uo,dx2,dy2] = WA_BL_accurate_prop2d_PartialCenterCenter_v2(h,z,dx1,lambda,q);
    figure(2)
    %uo = rot90(uo,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % show reconstruction and save
    mshow(abs(uo),2);
    
    % przeskaluj macierz na obraz w skali szarości
    absuo_norm = mat2gray(abs(uo));

    % Zapis macierzy typu uint8 (0-255) lub uint16 (0-65535) bezpośrednio do BMP:
    imwrite(absuo_norm , fullfile(nazwa_folderu{i}, 'wynik_czysty.bmp'));
    
    %exportgraphics(gcf, fullfile(nazwa_folderu,'rekonstrukcja.png'))

end


%% multipleksing

close all; 
clear;



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
Bx = 2*z*tan(FoV)*0.8;
By = 2*z*tan(FoV)*0.8;
ymax = By/2;

% occulsion parametes
wielk_siatki=3;
kernel=1;
start = 5.03;

Nxslm = 4160; Nyslm = 2464;

% odczytaj folder
mltpx_img = zeros(Nyslm,Nxslm);


%%

% odczytaj folder

for i = 1:9



    nazwa_folderu = {'zmienNazwe11','zmienNazwe12','zmienNazwe13','zmienNazwe21','zmienNazwe22','zmienNazwe23','zmienNazwe31','zmienNazwe32','zmienNazwe33'}; 
    %nazwa_folderu = {'interceptor11','interceptor12','interceptor13','interceptor21','interceptor22','interceptor23','interceptor31','interceptor32','interceptor33'}; 
    
    nazwy = {'multiplex1.bmp', 'multiplex2.bmp','multiplex3.bmp','multiplex4.bmp','multiplex5.bmp','multiplex6.bmp','multiplex7.bmp','multiplex8.bmp','multiplex9.bmp'};
    
    % Sprawdź czy folder istnieje. Jeśli nie (warunek tyldy '~'), utwórz go.
    if ~exist(nazwa_folderu{i}, 'dir')
        mkdir(nazwa_folderu{i});
    end
    
    %wczytaj obraz
    img = im2double(imread(fullfile(nazwa_folderu{i},'wynik_czysty.bmp')));
    
    
    %nie skaluj obrazy tylko odrazu sumuj amplitudy
    
    mltpx_img = mltpx_img + img;
    
    
    averaged_img = mltpx_img/i;
    
    
    imshow(averaged_img);
    
    save=mat2gray(averaged_img);
    
    wynik_folder = 'multiplexing';
    
    if ~exist(wynik_folder, 'dir')
        mkdir(wynik_folder);
    end
    
    imwrite(save, fullfile(wynik_folder, nazwy{i}));
end
