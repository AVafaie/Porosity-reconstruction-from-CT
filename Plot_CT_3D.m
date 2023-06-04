clc
clear
close all

% Select a plot type for illustrating porosity distribution
% Y: Yes
% N: No
phi_full='N';
phi_3quarter='Y';
phi_cross='N';
phi_slice='N';
porosity_histogram='Y';

% Select a plot type for illustrating permeability distribution
% Y: Yes
% N: No
permeability_full='N';
perm_3quarter='N';
permeability_slice='N';
permeability_histogram='Y';

% parameters for setting logarithmic colorbar of the permeability map
% note that in Matlab version 2010 in which the code is written, assigning
% a logarithmic colorbar is not straight forward and requires some
% calculations. To do so, the minimum and maximum values of permeability on
% the colorbar range have to be specified

% min permeability value assigned to the light color
color_min=1e-17;

% max permeability value assigned to the dark color
color_max=1e-11;



currentfolder = pwd;
addpath(currentfolder)

image_folder = [currentfolder,'\Stack\'];
Grid_data_folder = [currentfolder,'\Results\'];

load([Grid_data_folder,'3D_grid_porosity_adjusted.mat']);
load([Grid_data_folder,'3D_grid_permeability_adjusted.mat']);

load('mycmap');


%% full 3D plot of the porosity map
if phi_full=='Y'
    Fig_3D_phi_full=figure('Color',[1 1 1]);
    Fig_3D_phi_full.Units = 'centimeters';
    Fig_3D_phi_full.Position = [12 5 20 10];


    for k = 1:Cell_Z_No
        for j = 1:Cell_X_No
            for i = 1:Cell_Y_No
                Point_x=-[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0] * Cell_width + j*Cell_width;
                Point_y=-[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1] * Cell_width + i*Cell_width;
                Point_z=-[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1] * Cell_length + k*Cell_length;

                color = ones(4,6)*Porosity_3D_adjusted(i,j,k);

                if (isnan(Porosity_3D_adjusted(i,j,k))==0)
                    for l=1:6
                        patch(Point_x(:,l),Point_z(:,l),Point_y(:,l),color(:,l));
                    end
                end
            end
        end
    end

    view(45, 20);
    hold off
    axis equal
    axis off

    
    cbar = colorbar('Location','Eastoutside',...
        'FontName', 'Helvetica', ...
        'FontSize',14);
    set(get(cbar,'ylabel'),'string','Porosity, $\phi$ [-]',...
        'Rotation',90,...
        'FontName', 'Helvetica', ...
        'interpreter', 'latex',...
        'fontsize',16);
    
    caxis([0 0.7])
    
    map=mycmap;     
    colormap(map)      % parula  imcomplement(gray)  hot 
end;


%% 3Quarter 3D plot of the porosity map
if phi_3quarter=='Y'
    Fig_3D_phi_3quarter=figure('Color',[1 1 1]);
    Fig_3D_phi_3quarter.Units = 'centimeters';
    Fig_3D_phi_3quarter.Position = [12 5 20 10];

    for k = 1:Cell_Z_No
        for j = 1:Cell_X_No
            for i = 1:ceil(Cell_Y_No/2)
                Point_x=-[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0] * Cell_width + j*Cell_width;
                Point_y=-[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1] * Cell_width + i*Cell_width;
                Point_z=-[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1] * Cell_length + k*Cell_length;

                color = ones(4,6)*Porosity_3D_adjusted(i,j,k);

                if (isnan(Porosity_3D_adjusted(i,j,k))==0)
                    for l=1:6
                        patch(Point_x(:,l),Point_z(:,l),Point_y(:,l),color(:,l));
                    end
                end
            end
        end
    end

   for k = 1:Cell_Z_No
        for j = 1:ceil(Cell_X_No/2)
            for i = (ceil(Cell_Y_No/2)+1):Cell_Y_No
                Point_x=-[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0] * Cell_width + j*Cell_width;
                Point_y=-[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1] * Cell_width + i*Cell_width;
                Point_z=-[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1] * Cell_length + k*Cell_length;

                color = ones(4,6)*Porosity_3D_adjusted(i,j,k);

                if (isnan(Porosity_3D_adjusted(i,j,k))==0)
                    for l=1:6
                        patch(Point_x(:,l),Point_z(:,l),Point_y(:,l),color(:,l));
                    end
                end
            end
        end
    end

    view(50, 20);
    hold off
    axis equal
    axis off

    cbar = colorbar('Location','Eastoutside',...
        'FontName', 'Arial', ...
        'FontSize',10);
    set(get(cbar,'ylabel'),'string','Porosity, \phi [-]',...
        'Rotation',90,...
        'FontName', 'Arial', ...
        'interpreter', 'tex',...
        'fontsize',10);
    caxis([0 0.7])
    
    map=mycmap;     
    colormap(map)      % parula  imcomplement(gray)  hot
end;

%% 3D cross plot of the porosity map
if phi_cross=='Y'
    Fig_3D_phi_cross=figure('Color',[1 1 1]);
    Fig_3D_phi_cross.Units = 'centimeters';
    Fig_3D_phi_cross.Position = [12 5 20 10];


    for k = 1:Cell_Z_No
        for j = ceil(Cell_X_No/2)
            for i = 1:Cell_Y_No
                Point_x=-[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0] * Cell_width + j*Cell_width;
                Point_y=-[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1] * Cell_width + i*Cell_width;
                Point_z=-[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1] * Cell_length + k*Cell_length;

                color = ones(4,6)*Porosity_3D_adjusted(i,j,k);

                if (isnan(Porosity_3D_adjusted(i,j,k))==0)
                    for l=1:6
                        patch(Point_x(:,l),Point_z(:,l),Point_y(:,l),color(:,l));
                    end
                end
            end
        end
    end

    for k = 1:Cell_Z_No
        for j = 1:Cell_X_No
            for i = ceil(Cell_Y_No/2)
                Point_x=-[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0] * Cell_width + j*Cell_width;
                Point_y=-[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1] * Cell_width + i*Cell_width;
                Point_z=-[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1] * Cell_length + k*Cell_length;

                color = ones(4,6)*Porosity_3D_adjusted(i,j,k);

                if (isnan(Porosity_3D_adjusted(i,j,k))==0)
                    for l=1:6
                        patch(Point_x(:,l),Point_z(:,l),Point_y(:,l),color(:,l));
                    end
                end
            end
        end
    end

    view(45, 15);
    hold off
    axis equal
    axis off


    cbar = colorbar('Location','Eastoutside',...
        'FontName', 'Helvetica', ...
        'FontSize',14);
    set(get(cbar,'ylabel'),'string','Porosity, $\phi$ [-]',...
        'Rotation',90,...
        'FontName', 'Helvetica', ...
        'interpreter', 'latex',...
        'fontsize',16);
    
    caxis([0 0.7])
    
    map=mycmap;     
    colormap(map)      % parula  imcomplement(gray)  hot
end;



%% slice plot of the porosity map
if phi_slice=='Y'
    Fig_3D_phi_cross=figure('Color',[1 1 1]);
    Fig_3D_phi_cross.Units = 'centimeters';
    Fig_3D_phi_cross.Position = [12 5 20 10];


    for k = 1:Cell_Z_No
        for j = ceil(Cell_X_No/2-1)
            for i = 1:Cell_Y_No
                Point_x=-[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0] * Cell_width + j*Cell_width;
                Point_y=-[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1] * Cell_width + i*Cell_width;
                Point_z=-[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1] * Cell_length + k*Cell_length;

                color = ones(4,6)*Porosity_3D_adjusted(i,j,k);

                if (isnan(Porosity_3D_adjusted(i,j,k))==0)
                    for l=1:6
                        
                        % to draw a vertical cross section
                        patch(Point_x(:,l),Point_y(:,l),-Point_z(:,l),color(:,l),...
                            'LineWidth',0.1);                        
                    end
                end
            end
        end
    end



    view(45, 15);
    hold off
    axis equal
    axis off

   
     % colorbar for the methodological overview figure 
    cbar = colorbar('Position',...
    [[0.643200320225391 0.543516666666667 0.0244085427689588 0.383116666666667]],...
        'FontName', 'Arial', ...
        'FontSize',10);

    caxis([0 1])
    
    map=mycmap;     
    colormap(map)      % parula  imcomplement(gray)  hot
end;



%% Plot porosity histogram
if porosity_histogram=='Y'
    porosity_vector=reshape(Porosity_3D_adjusted,Cell_X_No*Cell_Y_No*Cell_Z_No,1,1);
        
    porosity_hist=figure('Color',[1 1 1]);
    porosity_hist.Units = 'centimeters';
    porosity_hist.Position = [12 2 6.7 6.3];

    bins = (0:0.05:1);
    hist(porosity_vector, bins);
    xlim([-0.025 1.025]);
    
    set(gca,'FontName','Arial','FontSize',10);
    
    ylabel('Frequency','FontSize',10,'FontName','Arial'); 
    xlabel('Porosity, \phi (-)','FontSize',10,'FontName','Arial');
    
    % Calculate standard deviation
    count=0;
    for i=1:length(porosity_vector)
        if (isnan(porosity_vector(i))==0)
            count=count+1;
            porosity_vector_value(count)=porosity_vector(i);
        end;
    end;
    
    disp(['Standard deviation of porosity= ',num2str(std(porosity_vector_value))])
    disp(['Average porosity= ',num2str(mean(porosity_vector_value))])
end


%% full 3D plot of the permeability map
if permeability_full=='Y'
    Fig_3D_perm_full=figure('Color',[1 1 1]);
    Fig_3D_perm_full.Units = 'centimeters';
    Fig_3D_perm_full.Position = [12 5 20 10];
    
    % as mentioned above, a permeability transformation is required to
    % assign the logarithmic colorbar
    transformed_color=(10.^(((log10(Permeability_3D_adjusted))-(log10(color_min)))/(log10(color_max)-log10(color_min))))*(color_max/10);

    for k = 1:Cell_Z_No
        for j = 1:Cell_X_No
            for i = 1:Cell_Y_No
                Point_x=-[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0] * Cell_width + j*Cell_width;
                Point_y=-[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1] * Cell_width + i*Cell_width;
                Point_z=-[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1] * Cell_length + k*Cell_length;

                color = ones(4,6)*transformed_color(i,j,k);

                if (isnan(Permeability_3D_adjusted(i,j,k))==0)
                    for l=1:6
                        patch(Point_x(:,l),Point_z(:,l),Point_y(:,l),color(:,l));
                    end
                end
            end
        end
    end

    view(45, 20);
    hold off
    axis equal
    axis off

    cbar = colorbar('Location','Eastoutside',...
        'FontName', 'Helvetica', ...
        'YScale','log',...
        'FontSize',14);
    set(get(cbar,'ylabel'),'string','Permeability, $k$ [$m^2$]',...
        'Rotation',90,...
        'FontName', 'Helvetica', ...
        'interpreter', 'latex',...
        'fontsize',16);
        
    caxis([color_min color_max])
    
    map=mycmap;     
    colormap(map)      % parula  imcomplement(gray)  hot
end



%% 3Quarter 3D plot of the permeability map
if perm_3quarter=='Y'
    Fig_3D_perm_3quarter=figure('Color',[1 1 1]);
    Fig_3D_perm_3quarter.Units = 'centimeters';
    Fig_3D_perm_3quarter.Position = [12 5 20 10];

    % as mentioned above, a permeability transformation is required to
    % assign the logarithmic colorbar
    transformed_color=(10.^(((log10(Permeability_3D_adjusted))-(log10(color_min)))/(log10(color_max)-log10(color_min))))*(color_max/10);

    for k = 1:Cell_Z_No
        for j = 1:Cell_X_No
            for i = 1:1:ceil(Cell_Y_No/2)
                Point_x=-[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0] * Cell_width + j*Cell_width;
                Point_y=-[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1] * Cell_width + i*Cell_width;
                Point_z=-[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1] * Cell_length + k*Cell_length;

                color = ones(4,6)*transformed_color(i,j,k);

                if (isnan(Permeability_3D_adjusted(i,j,k))==0)
                    for l=1:6
                        patch(Point_x(:,l),Point_z(:,l),Point_y(:,l),color(:,l),'LineWidth',0.3);
                    end
                end
            end
        end
    end
    
    for k = 1:Cell_Z_No
        for j = 1:ceil(Cell_X_No/2)
            for i = (ceil(Cell_Y_No/2)+1):Cell_Y_No
                Point_x=-[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0] * Cell_width + j*Cell_width;
                Point_y=-[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1] * Cell_width + i*Cell_width;
                Point_z=-[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1] * Cell_length + k*Cell_length;

                color = ones(4,6)*transformed_color(i,j,k);

                if (isnan(Permeability_3D_adjusted(i,j,k))==0)
                    for l=1:6
                        patch(Point_x(:,l),Point_z(:,l),Point_y(:,l),color(:,l),'LineWidth',0.3);
                    end
                end
            end
        end
    end

    view(45, 20);
    hold off
    axis equal
    axis off

    cbar = colorbar('Location','Eastoutside',...
        'FontName', 'Helvetica', ...
        'YScale','log',...
        'FontSize',14);
    
    set(get(cbar,'ylabel'),'string','Permeability, $k$ [$m^2$]',...
        'Rotation',90,...
        'FontName', 'Helvetica', ...
        'interpreter', 'latex',...
        'fontsize',16);
          
    caxis([color_min color_max])
    
    map=mycmap;     
    colormap(map)      % parula  imcomplement(gray)  hot
end;


%% slice plot of the permeability map
if permeability_slice=='Y'
    Fig_3D_perm_full=figure('Color',[1 1 1]);
    Fig_3D_perm_full.Units = 'centimeters';
    Fig_3D_perm_full.Position = [12 5 20 10];
    
    % as mentioned above, a permeability transformation is required to
    % assign the logarithmic colorbar
    transformed_color=(10.^(((log10(Permeability_3D_adjusted))-(log10(color_min)))/(log10(color_max)-log10(color_min))))*(color_max/10);


    for k = 1:Cell_Z_No
        for j = ceil(Cell_X_No/2)
            for i = 1:(Cell_Y_No)
                Point_x=-[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0] * Cell_width + j*Cell_width;
                Point_y=-[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1] * Cell_width + i*Cell_width;
                Point_z=-[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1] * Cell_length + k*Cell_length;

                color = ones(4,6)*transformed_color(i,j,k);

                if (isnan(Permeability_3D_adjusted(i,j,k))==0)
                    for l=1:6
%                         patch(Point_x(:,l),Point_z(:,l),Point_y(:,l),color(:,l));
                        
                        % to draw a vertical cross section
                        patch(Point_x(:,l),Point_y(:,l),-Point_z(:,l),color(:,l),...
                            'LineWidth',0.1);
                    end
                end
            end
        end
    end

    view(45, 20);
    hold off
    axis equal
    axis off
    

    cbar = colorbar('Location','Eastoutside',...
        'FontName', 'Helvetica', ...
        'YScale','log',...
        'FontSize',14);
    
    set(get(cbar,'ylabel'),'string','Permeability, $k$ [$m^2$]',...
        'Rotation',90,...
        'FontName', 'Helvetica', ...
        'interpreter', 'latex',...
        'fontsize',16);
          
    caxis([color_min color_max])
    
    map=mycmap;     
    colormap(map)      % parula  imcomplement(gray)  hot
end;


%% Plot permeability histogram
if permeability_histogram=='Y'
    permeability_vector=reshape(Permeability_3D_adjusted,Cell_X_No*Cell_Y_No*Cell_Z_No,1,1);
        
    permeability_hist=figure('Color',[1 1 1]);
    permeability_hist.Units = 'centimeters';
    permeability_hist.Position = [12 2 6.7 6.3];

    bins = 10.^(-22:0.3:-10);
    hist(permeability_vector, bins);
    xlim([1e-20 1e-10]);
   
    set(gca,'FontName','Arial','FontSize',10,'XMinorTick','on','XScale','log');
    
    ylabel('Frequency','FontSize',10,'FontName','Arial'); 
    xlabel('Permeability (m^2)','FontName','Arial');

        % Calculate standard deviation
    count=0;
    for i=1:length(permeability_vector)
        if (isnan(permeability_vector(i))==0)
            count=count+1;
            permeability_vector_value(count)=permeability_vector(i);
        end;
    end;
    
    disp(['Standard deviation of permeability= ',num2str(std(permeability_vector_value)),'(m2)'])
    disp(['Average permeability= ',num2str(mean(permeability_vector_value)),'(m2)'])
end



     
