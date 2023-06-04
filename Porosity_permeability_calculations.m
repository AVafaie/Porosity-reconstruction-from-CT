clc;  clearvars;  close all


tic


currentfolder = pwd;
addpath(currentfolder)

image_folder = [currentfolder,'\Stack\Binary\'];
save_image_folder = [currentfolder,'\Stack\'];


% Increasing CT porosity to the effective porosity value
% by randomly distributing the sub-resolution pore on the solid part
% Y: yes 
% N: no
Porosity_Adjus='Y';


%% Initial rock properties
%initial average porosity (from Helium porosimeter)
Por_Eff=0.28;

% input for CRUNCHFLOW CODE
%initial average permeability of the rock
Perm_Eff=3.8e-14;

%pump zone permeability
Perm_Pump=1.0e-9;

%minimum porosity
Por_Min=0.001;

%maximum porosity
Por_Max=1;

%permeability of inert zone
Perm_Inert=1.0e-22;

%permeability of the ghost zones
Perm_Ghost=1.0e-22;

% power of porosity in permeability-porosity relationship to calculate
% permeability of each cell
permeability_power=3;


%% Griding information
%number of grids in the X direction
Cell_X_No=20;

% number of grids in the Y direction
Cell_Y_No=20;

% number of grids in the Z direction
Cell_Z_No=15;

% core length and diameter (width of CT cross-section in mm)
Core_Diameter=25.3;
Core_length=53;

% width of each cell in millimeter (x and y direction)
Cell_width=Core_Diameter/Cell_X_No;

% length of each cell in millimeter (z direction)
Cell_length=Core_length/Cell_Z_No;


%% Loading images and calculate initial CT porosity

for z=1:Cell_Z_No
    % loading image names for each grid set in the Z direction
    image_information=dir([image_folder,'Z',num2str(z),'\','*.tif']);
    
    numslice=numel(image_information);
    
    for i=1:numslice
        % reading binary images associated to each Z
        t=imread([image_folder,'Z',num2str(z),'\', image_information(i).name]);
        

        % convert the image to 0,1
        t=(t~=0);

        % 3D stack of brithness matrices 
        % can be viewed in 3D using volumeViewer(IMG)
        % 1: solid, 0: pore
        IMG(:,:,i)=t;
        
	end;
     
    % equal in both directions here
    Pixel_No = length(IMG);

    % pixel width in millimeter
    Pixel_Width=Core_Diameter/Pixel_No;

    
    % Defining the core cross section on the loaded images
    % we put the origin at the left-uppermost point
    
    Section_Domain=zeros(Pixel_No,Pixel_No);
    for m=1:Pixel_No
        for n=1:Pixel_No
            X_Pixel=(n-1/2)*Pixel_Width;
            Y_Pixel=(m-1/2)*Pixel_Width;
            
            if (X_Pixel-Core_Diameter/2)^2+(Y_Pixel-Core_Diameter/2)^2<=(Core_Diameter/2)^2
                % 1 for the part of images containing the rock
                Section_Domain(m,n)=1;
            end
        end;
    end;

    % initial porosity of all slices along each Z grid
    % a matrix of (number_of_zgrid * number_of_slice) 
    for s=1:numslice
        Porosity_slice_total(s,z)=(Pixel_No*Pixel_No-sum(sum(IMG(:,:,s))))/sum(sum(Section_Domain(:,:)));
    end;

    % save intensity matrices for the whole sample (separated by zgrid)
    if exist([save_image_folder,'Intensity_matrices\'])==7
        save( [save_image_folder,'Intensity_matrices\','Intensity_','Z',num2str(z)],'IMG');
    else
        mkdir(save_image_folder,'Intensity_matrices');
        save( [save_image_folder,'Intensity_matrices\','Intensity_','Z',num2str(z)],'IMG');
    end;
end;

% average CT porosity of each Z grid
Porosity_zgrid_total=mean(Porosity_slice_total);

% total sample porosity from CT
Porosity_total=mean(mean(Porosity_zgrid_total));
    
% save porosity for all slices, Sectiodomain,
% Porosity_zgrid_total and total porosity
if exist([currentfolder,'\Results'])==7
    save( [currentfolder,'\Results','\Initial_CT_porosity'],'Porosity_slice_total','Section_Domain','Porosity_zgrid_total','Porosity_total');
else
    mkdir(currentfolder,'Results');
    save( [currentfolder,'\Results','\Initial_CT_porosity'],'Porosity_slice_total','Section_Domain','Porosity_zgrid_total','Porosity_total');
end;


toc


%% Porosity adjustment from CT to effective
if Porosity_Adjus=='Y';
tic 
    if exist([currentfolder,'\Results','\Initial_CT_porosity.mat'])==2
        load([currentfolder,'\Results','\Initial_CT_porosity.mat'])
    else
        disp('Porosity calculations have to be performed first')
    end;
        
    for z=1:Cell_Z_No
        if exist([save_image_folder,'Intensity_matrices\','Intensity_','Z',num2str(z),'.mat'])==2
            load([save_image_folder,'Intensity_matrices\','Intensity_','Z',num2str(z),'.mat'])
        else
            disp('Porosity calculations have to be performed first')
        end;
        
        [Pixel_No Pixel_No numslice]=size(IMG);
        
        % number of pixels that has to be added to each CT cross section to reproduce
        % effective porosity of the rock
        Delta_Pixel=fix((Por_Eff-Porosity_total)*sum(sum(Section_Domain(:,:))));
    
        % We add "Delta_pixel" pixels to all CT slices of all Zgrids
        for s=1:numslice
            % Matrix containing pixel positions for the solid part of the rock
            % 1 for the matrix part
            Solid_Pixels=IMG(:,:,s).*Section_Domain;

            % Returns rows and columns that have value 1 (Solid) in Solid_pixels  
            [R_mat,C_mat] = find(Solid_Pixels==1);

            % select randomly "Delta_pixel" pixels
            idx = randsample(length(C_mat),Delta_Pixel);

            % In the brightness matrix IMG, replace matrix points with 0 (add porosity)
            for I=1:length(idx)
                IMG(R_mat(idx(I)),C_mat(idx(I)),s)=0;
            end;
            % final porosity of each section after adding random pores
            Porosity_slice_total_updated(s,z)=(Pixel_No*Pixel_No-sum(sum(IMG(:,:,s))))/sum(sum(Section_Domain(:,:))); 
        end;
        
        % save updated intensity matrices for the whole sample 
        if exist([save_image_folder,'Intensity_matrices_updated'])==7
            save( [save_image_folder,'Intensity_matrices_updated\','Intensity_','Z',num2str(z),'_updated'],'IMG');
        else
            mkdir(save_image_folder,'Intensity_matrices_updated');
            save( [save_image_folder,'Intensity_matrices_updated\','Intensity_','Z',num2str(z),'_updated'],'IMG');
        end;
    end;  
    
    % average CT porosity of each Z grid (updated)
    Porosity_zgrid_total_updated=mean(Porosity_slice_total_updated);

    % total sample porosity from CT (updated)
    Porosity_total_updated=mean(mean(Porosity_zgrid_total_updated));

    % save updated porosity for all slices and total porosity
    save( [currentfolder,'\Results','\Updated_CT_porosity'],'Porosity_slice_total_updated','Porosity_zgrid_total_updated','Porosity_total_updated');

    toc  
end;

  

%% buliding permeability files in each direction (CrunchFlow code)
tic
%cell_X_no+2 number of cells plus 2 ghost cells in X direction
% 4: we need a matrix of (cell_X_no+2)*(cell_Y_no)*(Cell_Z_No) rows and 4
% columns
Single_Column_Perm_X=zeros((Cell_X_No+2)*(Cell_Y_No)*(Cell_Z_No+1),4)+Perm_Ghost;

for Z=1:Cell_Z_No+1
    for Y=1:Cell_Y_No
        for X=1:Cell_X_No+2
            Single_Column_Perm_X(((Cell_X_No+2)*(Cell_Y_No))*(Z-1)+((Y-1)*(Cell_X_No+2))+X,1)=X-1;
            Single_Column_Perm_X(((Cell_X_No+2)*(Cell_Y_No))*(Z-1)+((Y-1)*(Cell_X_No+2))+X,2)=Y;
            Single_Column_Perm_X(((Cell_X_No+2)*(Cell_Y_No))*(Z-1)+((Y-1)*(Cell_X_No+2))+X,3)=Z;
        end 
    end
end

%cell_Y_no+2 number of cells plus 2 ghost cells in Y direction
% 4: we need a matrix of (cell_X_no)*(cell_Y_no+2)*(numimage) rows and 4
% columns
Single_Column_Perm_Y=zeros((Cell_X_No)*(Cell_Y_No+2)*(Cell_Z_No+1),4)+Perm_Ghost;

for Z=1:Cell_Z_No+1
    for Y=1:Cell_Y_No+2
        for X=1:Cell_X_No
            Single_Column_Perm_Y(((Cell_X_No)*(Cell_Y_No+2))*(Z-1)+((Y-1)*(Cell_X_No))+X,1)=X;
            Single_Column_Perm_Y(((Cell_X_No)*(Cell_Y_No+2))*(Z-1)+((Y-1)*(Cell_X_No))+X,2)=Y-1;
            Single_Column_Perm_Y(((Cell_X_No)*(Cell_Y_No+2))*(Z-1)+((Y-1)*(Cell_X_No))+X,3)=Z;
        end 
    end
end

%numimage+2 number of cells plus 2 ghost cells in Z dircetion
% 4: we need a matrix of (cell_X_no+2)*(cell_Y_no)*(numimage+2) rows and 4
% columns
Single_Column_Perm_Z=zeros((Cell_X_No)*(Cell_Y_No)*(Cell_Z_No+3),4)+Perm_Ghost;

for Z=1:Cell_Z_No+3
    for Y=1:Cell_Y_No
        for X=1:Cell_X_No
            Single_Column_Perm_Z(((Cell_X_No)*(Cell_Y_No))*(Z-1)+((Y-1)*(Cell_X_No))+X,1)=X;
            Single_Column_Perm_Z(((Cell_X_No)*(Cell_Y_No))*(Z-1)+((Y-1)*(Cell_X_No))+X,2)=Y;
            Single_Column_Perm_Z(((Cell_X_No)*(Cell_Y_No))*(Z-1)+((Y-1)*(Cell_X_No))+X,3)=Z-1;
        end 
    end
end

toc

%% Calculating porosity values in all grid cells
tic
% 2D matrix file of porosity for permeability calculations
Porosity_vector=zeros(Cell_Y_No*Cell_X_No,Cell_Z_No);

% 3D array for plotting digitized porosity cross sections
Porosity_2D=zeros(Cell_Y_No,Cell_X_No,Cell_Z_No);
    
% 3D array for plotting digitized porosity in 3D
Porosity_3D=zeros(Cell_Y_No,Cell_X_No,Cell_Z_No);

% 2D matrix of grid coordinates (X, Y)
Grid_coordinates=zeros(Cell_Y_No,Cell_X_No);
    

for z=1:Cell_Z_No
    % read binary images from image folders
    if Porosity_Adjus=='Y';
        if exist([save_image_folder,'Intensity_matrices_updated'])==7
            load([save_image_folder,'Intensity_matrices_updated\','Intensity_','Z',num2str(z),'_updated','.mat'])
        else
            disp('Porosity updating has to be performed first')
        end;
    else
        if exist([save_image_folder,'Intensity_matrices'])==7
            load([save_image_folder,'Intensity_matrices\','Intensity_','Z',num2str(z),'.mat'])
        else
            disp('Porosity calculations have to be performed first')
        end;
    end;
    
    [Pixel_No Pixel_No numslice]=size(IMG);

    % I increasing from left to right and from bottom to top
    for I=1:Cell_X_No*Cell_Y_No
        Cell_M=ceil(I/Cell_X_No);  %row number of cell I 

        Cell_N=rem(I,Cell_X_No);   %column number of cell I
        
        if Cell_N==0
            Cell_N=Cell_X_No;
        end;

        % number of the pixel located on the top row of the cells
        Pixel_Row_Up=Pixel_No-ceil(Cell_M*Cell_width/Pixel_Width);
        
        if Pixel_Row_Up<=0
            Pixel_Row_Up=1;
        end;
        
        % number of the pixel located on the bottom row of the cells
        Pixel_Row_Down=Pixel_No-ceil((Cell_M-1)*Cell_width/Pixel_Width);
        
        % number of the pixel located on the left column of the cells
        if Cell_N==1
            Pixel_Column_Left=1;
        else
            Pixel_Column_Left=ceil((Cell_N-1)*Cell_width/Pixel_Width);
        end;
        
        % number of the pixel located on the right column of the cells
        Pixel_Column_Right=ceil(Cell_N*Cell_width/Pixel_Width);
        
        if Pixel_Column_Right>=Pixel_No
            Pixel_Column_Right=Pixel_No;
        end;

        % exporting grid coordinates for plotting (X,Y)
        % we assume that the coordinates are the same along all z grids
        % X=0 on the left
        % Y=0 at the bottom
        if z==1
            Grid_coordinates(Cell_M,Cell_N,1)=(Cell_N-1)*Cell_width+Cell_width/2;
            Grid_coordinates(Cell_M,Cell_N,2)=(Cell_M-1)*Cell_width+Cell_width/2;
            
            % see if the grid is on rock or the surrounding media
            % 0 means surrounding media
            Section_domain_status(Cell_M,Cell_N)=sum(sum(Section_Domain(Pixel_Row_Up:Pixel_Row_Down,Pixel_Column_Left: Pixel_Column_Right)));
        end;

        %porosity of cell I
        for s=1:numslice
            Cell_Phi_slice(s)=(((Pixel_Row_Down-Pixel_Row_Up+1)*(Pixel_Column_Right-Pixel_Column_Left+1))-sum(sum(IMG(Pixel_Row_Up:Pixel_Row_Down,Pixel_Column_Left: Pixel_Column_Right,s))))/((Pixel_Row_Down-Pixel_Row_Up+1)*(Pixel_Column_Right-Pixel_Column_Left+1));
        end;
        
        Cell_Phi(I)=mean(Cell_Phi_slice);
        clear('Cell_Phi_slice');
        
        % Porosity file for permeability measurements
        % columns: Z grid
        % rows: porosity of grids (Cell_X_No*Cell_Y_No)
        Porosity_vector(I,z)=Cell_Phi(I);
        
        % Convert porosity to brightness for plotting purposes
        % Note that the order of cell rows are inverted for 2D plotting
        Porosity_2D(Cell_Y_No-Cell_M+1,Cell_N,z)=Cell_Phi(I);
        
        % Cell_Y_No*Cell_X_No matrix of porosity (for 3D porosity plot)
        Porosity_3D(Cell_M,Cell_N,z)=Cell_Phi(I);
        
        % the area or volume surrounding the rock is eliminated from plots
        if Section_domain_status(Cell_M,Cell_N)==0
            Porosity_2D(Cell_Y_No-Cell_M+1,Cell_N,z)=NaN;
            Porosity_3D(Cell_M,Cell_N,z)=NaN;
        end 
        
    end
 end;
 
% save 2D and 3D grid porosity
if Porosity_Adjus=='Y'
    Porosity_2D_adjusted=Porosity_2D;
    Porosity_3D_adjusted=Porosity_3D;
    
    save( [currentfolder,'\Results','\2D_grid_porosity_adjusted'],'Porosity_2D_adjusted','Cell_X_No','Cell_Y_No','Cell_Z_No');
    save( [currentfolder,'\Results','\3D_grid_porosity_adjusted'],'Porosity_3D_adjusted','Cell_X_No','Cell_Y_No','Cell_Z_No','Cell_width','Cell_length');

else
    Porosity_2D_raw=Porosity_2D;
    Porosity_3D_raw=Porosity_3D;
    
    save( [currentfolder,'\Results','\2D_grid_porosity_raw'],'Porosity_2D_raw','Cell_X_No','Cell_Y_No','Cell_Z_No');
    save( [currentfolder,'\Results','\3D_grid_porosity_raw'],'Porosity_3D_raw','Cell_X_No','Cell_Y_No','Cell_Z_No','Cell_width','Cell_length');
end;

toc

%% Permeability calculations from the porosity map
tic
% output according to the CRUNCHFLOW code format
for z=1:Cell_Z_No
    for I=1:Cell_X_No*Cell_Y_No

        % Assigning a lower porosity threshold for numerical simulations
        % of flow
        if Porosity_vector(I,z)<Por_Min
            Porosity_vector(I,z)=Por_Min;
        end;
        
        % Determines in which column we are
        R=rem(I,Cell_X_No);

        % how many rows are filled
        Q=(I-R)/Cell_X_No;
        
        % permeability in each direction by considering ghost cells 
        if R==0
           Q=Q-1;
           R=Cell_X_No;
        end;
        
        % calculating permeability using a powerlaw equation (pump zone)
        if z==1
            if Porosity_vector(I,z)==Por_Min
               Single_Column_Perm_X((z-1)*((Cell_X_No+2)*Cell_Y_No)+Q*(Cell_X_No+2)+R+1,4)=Perm_Inert;
            else 
               Single_Column_Perm_X((z-1)*((Cell_X_No+2)*Cell_Y_No)+Q*(Cell_X_No+2)+R+1,4)=Perm_Pump; 
            end;

            if Porosity_vector(I,z)==Por_Min
                Single_Column_Perm_Y((z-1)*(Cell_X_No*(Cell_Y_No+2))+(Q+1)*(Cell_X_No)+R,4)=Perm_Inert;
              else 
                Single_Column_Perm_Y((z-1)*(Cell_X_No*(Cell_Y_No+2))+(Q+1)*(Cell_X_No)+R,4)=Perm_Pump; 
            end;

            if Porosity_vector(I,z)==Por_Min
               Single_Column_Perm_Z((z)*(Cell_X_No*(Cell_Y_No))+I,4)=Perm_Inert;
            else 
               Single_Column_Perm_Z((z)*(Cell_X_No*(Cell_Y_No))+I,4)=Perm_Pump; 
            end;
        end;
        
        % calculating permeability using a powerlaw equation (other cells)
        if Porosity_vector(I,z)==Por_Min
           Single_Column_Perm_X((z)*((Cell_X_No+2)*Cell_Y_No)+Q*(Cell_X_No+2)+R+1,4)=Perm_Inert;
        else 
           Single_Column_Perm_X((z)*((Cell_X_No+2)*Cell_Y_No)+Q*(Cell_X_No+2)+R+1,4)=Perm_Eff*((Porosity_vector(I,z)/Por_Eff)^permeability_power)*((1-Por_Eff)^2/(1-Porosity_vector(I,z))^2); 
        end;

        if Porosity_vector(I,z)==Por_Min
            Single_Column_Perm_Y((z)*(Cell_X_No*(Cell_Y_No+2))+(Q+1)*(Cell_X_No)+R,4)=Perm_Inert;
          else 
            Single_Column_Perm_Y((z)*(Cell_X_No*(Cell_Y_No+2))+(Q+1)*(Cell_X_No)+R,4)=Perm_Eff*((Porosity_vector(I,z)/Por_Eff)^permeability_power)*((1-Por_Eff)^2/(1-Porosity_vector(I,z))^2); 
        end;

        if Porosity_vector(I,z)==Por_Min
           Single_Column_Perm_Z((z+1)*(Cell_X_No*(Cell_Y_No))+I,4)=Perm_Inert;
        else 
           Single_Column_Perm_Z((z+1)*(Cell_X_No*(Cell_Y_No))+I,4)=Perm_Eff*((Porosity_vector(I,z)/Por_Eff)^permeability_power)*((1-Por_Eff)^2/(1-Porosity_vector(I,z))^2); 
        end;
     
        % Single_Column_Perm_Z(1:cell_X_no*cell_Y_no,4)=Single_Column_Perm_Z((cell_X_no*cell_Y_no+1):2*cell_X_no*cell_Y_no,4);
        Single_Column_Perm_Z(((Cell_Z_No+2)*Cell_X_No*Cell_Y_No+1):(Cell_Z_No+3)*Cell_X_No*Cell_Y_No,4)=Single_Column_Perm_Z(((Cell_Z_No+1)*Cell_X_No*Cell_Y_No+1):(Cell_Z_No+2)*Cell_X_No*Cell_Y_No,4);

    end;
end

% save 3D grid permeability (only for plotting)
Permeability_3D=Perm_Eff*((Porosity_3D/Por_Eff).^permeability_power).*((1-Por_Eff)^2./(1-Porosity_3D).^2);


if Porosity_Adjus=='Y'
    Permeability_3D_adjusted=Permeability_3D;
    save( [currentfolder,'\Results','\3D_grid_permeability_adjusted'],'Permeability_3D_adjusted');
else
    Permeability_3D_raw=Permeability_3D;
    save( [currentfolder,'\Results','\3D_grid_permeability_raw'],'Permeability_3D_raw');
end;


% Single column porosity ofr plotting in surfer
for i=1:Cell_Z_No
    Porosity_3D_surfer(:,:,i)=transpose(Porosity_3D(:,:,i));
end;
% only the rock section
Single_Column_Porosity=reshape(Porosity_3D_surfer,[Cell_X_No*Cell_Y_No*Cell_Z_No,1,1]);

Single_Column_Porosity(isnan(Single_Column_Porosity))=Por_Min;

% save single column porosity
save([currentfolder,'\Results','\single_column_porosity'],'Single_Column_Porosity');


%% Exporting permeability data into Excel files

% permeability files in each direction
if Porosity_Adjus=='Y'
    xlswrite( [currentfolder,'\Results','\PermX_adjusted'],Single_Column_Perm_X);
    xlswrite( [currentfolder,'\Results','\PermY_adjusted'],Single_Column_Perm_Y);
    xlswrite( [currentfolder,'\Results','\PermZ_adjusted'],Single_Column_Perm_Z);
else
    xlswrite( [currentfolder,'\Results','\PermX_raw'],Single_Column_Perm_X);
    xlswrite( [currentfolder,'\Results','\PermY_raw'],Single_Column_Perm_Y);
    xlswrite( [currentfolder,'\Results','\PermZ_raw'],Single_Column_Perm_Z);
end;

toc


memory


Mem_data=whos;

for z=1:length(Mem_data)
    Mem_size(z,1)=Mem_data(z).bytes;
end;

sum(Mem_size)



