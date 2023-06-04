clc
clear
close all

currentfolder = pwd;
addpath(currentfolder)

image_folder = [currentfolder,'\Stack\'];
Grid_data_folder = [currentfolder,'\Results\'];

% loading porosity data at the grid level, initial and updated intensity
% matrices
load([Grid_data_folder,'2D_grid_porosity_adjusted.mat']);
load([Grid_data_folder,'Initial_CT_porosity.mat']);
load([Grid_data_folder,'Updated_CT_porosity.mat']);

% Plot properties
image_position=[12 5 16.5 10];
background_color=[1 1 1];
image_title_fontsize=11;
subplot_title_fontsize=10;
colorbar_fontsize=10;
colorbar_title_fontsize=11;


% design the plot structure
if Cell_Z_No<=10
    row=2;
    column=ceil(Cell_Z_No/row);
elseif Cell_Z_No<20
    row=3;
    column=ceil(Cell_Z_No/row);
elseif Cell_Z_No<30
    row=4;
    column=ceil(Cell_Z_No/row);
end

% load colormap
load('mycmap');



%% 2D plots of raw grayscale images
% plot one slide (the first one in the stack) as an example
Raw_grayscale_images=figure('Color',[1 1 1]);
Raw_grayscale_images.Units = 'centimeters';
Raw_grayscale_images.Position = image_position;

clf;  

for z=1:Cell_Z_No
    image_information=dir([image_folder,'Grayscale\','Z',num2str(z),'\','*.bmp']);
    
    Original_image=imread([image_folder,'Grayscale\','Z',num2str(z),'\', image_information(1).name]);
    
    subplot(row,column,z);
    imagesc(Original_image); axis equal tight; 
    title(['Grid ',num2str(z)],'FontName', 'Helvetica', ...
        'FontSize',subplot_title_fontsize, 'interpreter', 'latex');
    colormap(gray)
    drawnow;
    
    set(gca,'DataAspectRatio',[1 1 1],...
    'XTick',zeros(1,0),'YTick',zeros(1,0));
    
end;

%annotation(Raw_grayscale_images,'textbox',...
annotation('textbox',...
    [0.36 0.9 0.1 0.1],...
    'String',{'Original CT images'},...
    'FontName', 'Helvetica', ...
    'interpreter', 'latex',...
    'fontsize',image_title_fontsize,...
    'EdgeColor','none');


%% 2D plots of raw binary CT images
% plot one slide (the first one in the stack) as an example
Binary_images_raw=figure('Color',[1 1 1]);;
Binary_images_raw.Units = 'centimeters';
Binary_images_raw.Position = image_position;

clf;  

for z=1:Cell_Z_No
    load([image_folder,'Intensity_matrices\','Intensity_Z',num2str(z),'.mat']);
    subplot(row,column,z);
    imagesc(IMG(:,:,1)); axis equal tight; 
    title(['Grid ',num2str(z),' ($\phi$ = ',num2str(Porosity_zgrid_total(z),'%5.3f'),')'],...
        'FontName', 'Helvetica', ...
        'FontSize',subplot_title_fontsize, 'interpreter', 'latex');
    colormap(gray)
    drawnow;
    
    set(gca,'DataAspectRatio',[1 1 1],...
        'XTick',zeros(1,0),'YTick',zeros(1,0));
    
end;

% annotation(Binary_images_raw,'textbox',...
annotation('textbox',...
    [0.36 0.9 0.1 0.1],...
    'String',{'Raw binary images'},...
    'FontName', 'Helvetica', ...
    'interpreter', 'latex',...
    'fontsize',image_title_fontsize,...
    'EdgeColor','none');

   
%% 2D plots of updated binary CT images
% plot one slide (the first one in the stack) as an example
Binary_images_updated=figure('Color',[1 1 1]);;
Binary_images_updated.Units = 'centimeters';
Binary_images_updated.Position = image_position;

clf;  

for z=1:Cell_Z_No
    load([image_folder,'Intensity_matrices_updated\','Intensity_Z',num2str(z),'_updated.mat']);
    subplot(row,column,z);
    imagesc(IMG(:,:,1)); axis equal tight; 
    title(['Grid ',num2str(z),' ($\phi$ = ',num2str(Porosity_zgrid_total_updated(z),'%5.3f'),')'],...
        'FontName', 'Helvetica', ...
        'FontSize',subplot_title_fontsize, 'interpreter', 'latex');
    colormap(gray)
    drawnow;
    
    set(gca,'DataAspectRatio',[1 1 1],...
        'XTick',zeros(1,0),'YTick',zeros(1,0));
    
end;

%annotation(Binary_images_updated,'textbox',...
annotation('textbox',...
    [0.34 0.9 0.1 0.1],...
    'String',{'Updated binary images'},...
    'FontName', 'Helvetica', ...
    'interpreter', 'latex',...
    'fontsize',image_title_fontsize,...
    'EdgeColor','none');


%% Plot 2D digitized porosity maps
    Digitized_porosity=figure('Color',[1 1 1]);;
    Digitized_porosity.Units = 'centimeters';
    Digitized_porosity.Position = image_position;

    clf;

for z=1:Cell_Z_No
    imAlpha=ones(Cell_X_No);
    imAlpha(isnan(Porosity_2D_adjusted(:,:,z)))=0;
    
    subplot(row,column,z);

    imagesc(Porosity_2D_adjusted(:,:,z),'AlphaData',imAlpha);  axis equal tight; 
    
    title(['Grid ',num2str(z)],'FontName', 'Helvetica', ...
        'FontSize',subplot_title_fontsize, 'interpreter', 'latex');
    
    map=mycmap;     
    colormap(map)       
    drawnow;
    
    set(gca,'DataAspectRatio',[1 1 1],...
        'XTick',zeros(1,0),'YTick',zeros(1,0));
    
end;

% annotation(Digitized_porosity,'textbox',...
annotation('textbox',...
    [0.34 0.9 0.1 0.1],...
    'String',{'Digitized CT map'},...
    'FontName', 'Helvetica', ...
    'interpreter', 'latex',...
    'fontsize',image_title_fontsize,...
    'EdgeColor','none');

cbar = colorbar('Location','Eastoutside',...
    'FontName', 'Helvetica', ...
    'FontSize',colorbar_fontsize);
set(get(cbar,'ylabel'),'string','Porosity, \phi [-]',...
    'FontName', 'Helvetica', ...
    'fontsize',colorbar_title_fontsize);


