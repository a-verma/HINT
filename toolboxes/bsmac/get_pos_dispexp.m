function ddat = get_pos_dispexp(command)

% command is the current axis being clicked on

% GET_POS: get current cursor positon, and display data information in the
% activity map display GUI.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Lijun Zhang, Ph.D                        %
% Center for Biomedical Imaging Statistics (CBIS)  %
% Dept. of Biostatistics and Bioinformatics        %
% Emory University                                 %
% Atlanta, GA, 30322                               %
% Email: l.zhang@emory.edu                         %
% http://www.sph.emory.edu/bios/CBIS/about.html    %
% Jan, 2009, Last revision: Dec. 1, 2011           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified by Joshua Lukemire to extend to multiple populations

global ddat;
cp = get(get(gcf,'CurrentAxes'),'CurrentPoint');
mri = ddat.mri_struct;
[sagDim, corDim, axiDim] = size(mri.img);
coordindex = 1;

pp = ddat.roi_voxel;

% Loop over each of the population images. For most cases this will just be
% one population. The exception is the subpopulation display window.
for i = 1:ddat.nCompare

    switch command
        % User clicked on the sagittal image
        case 'sag'
            % Get current coronal and axial positions.
            cor = round(cp(1,1));
            axi = round(cp(1,2));
            % Make sure nothing out of bounds
            if axi>axiDim
                axi = axiDim;
            elseif axi <1
                axi = 1;
            end
            if cor > corDim
                cor = corDim;
            elseif cor <1
                cor = 1;
            end
            % Update the data handle
            ddat.cor = cor;
            ddat.axi = axi;

            set(ddat.sagittal_xline{i},'Ydata',[axi axi]);
            set(ddat.sagittal_yline{i},'Xdata',[cor cor]);
            set(ddat.coronal_xline{i},'Ydata',[axi axi]);
            set(ddat.axial_xline{i},'Ydata',[cor cor]);

            for cl = 1:3
                Saxial(:, :, cl) = squeeze(ddat.combinedImg{i}(cl).combound(:, :, axi))';
                Scor(:, :, cl) = squeeze(ddat.combinedImg{i}(cl).combound(:, cor,:))';
            end
            set(ddat.axial_image{i},'CData',Saxial);
            set(ddat.coronal_image{i},'CData',Scor);
            set(findobj('Tag', 'AxiSlider'),'Value',axi);
            set(findobj('Tag', 'CorSlider'),'Value',cor);
            if coordindex == 1,
                set(findobj('Tag', 'crosshairPos'), 'String', sprintf('%7.0d %7.0d %7.0d',ddat.sag, cor, axi));
            elseif coordindex == 2,
                cor_mm = (cor-ddat.origin(2))*ddat.pixdim(2);
                axi_mm = (axi-ddat.origin(3))*ddat.pixdim(3);    
                set(ddat.crosshairPos, 'String', sprintf('%7.1f %7.1f %7.1f',ddat.sag, cor_mm, axi_mm));
            end 

        case 'cor'
            sag = round(cp(2,1));
            axi = round(cp(2,2));
            if axi>axiDim
                axi = axiDim;
            elseif axi <1
                axi = 1;
            end
            if sag > sagDim
                sag = sagDim;
            elseif sag <1
                sag = 1;
            end
            ddat.sag = sag;
            ddat.axi = axi;

            set(ddat.coronal_xline{i},'Ydata',[axi axi]);
            set(ddat.coronal_yline{i},'Xdata',[sag sag]);
            set(ddat.axial_yline{i},'Xdata',[sag sag]);
            set(ddat.sagittal_xline{i},'Ydata',[axi axi]); 
            for cl = 1:3
                Saxial(:, :, cl) = squeeze(ddat.combinedImg{i}(cl).combound(:, :, axi))';
                Ssag(:, :, cl) = squeeze(ddat.combinedImg{i}(cl).combound(sag, :, :))';
            end
            set(ddat.axial_image{i},'CData',Saxial);
            set(ddat.sagittal_image{i},'CData',Ssag);

            set(findobj('Tag', 'AxiSlider'),'Value',axi);
            set(findobj('Tag', 'SagSlider'),'Value',sag); 

            if coordindex == 1,
                set(findobj('Tag', 'crosshairPos'), 'String', sprintf('%7.0d %7.0d %7.0d',sag, ddat.cor, axi));
            elseif coordindex == 2,
                sag_mm = (sag-ddat.origin(1))*ddat.pixdim(1);
                axi_mm = (axi-ddat.origin(3))*ddat.pixdim(3);    
                set(findobj('Tag', 'crosshairPos'), 'String', sprintf('%7.1f %7.1f %7.1f',sag_mm, ddat.cor, axi_mm));
            end 

        case 'axi'
            sag = round(cp(2,1));
            cor = round(cp(2,2));
            if cor>corDim
                cor = corDim;
            elseif cor <1
                cor = 1;    
            end
            if sag > sagDim
                sag = sagDim;
            elseif sag <1
                sag = 1;    
            end
            ddat.sag = sag;
            ddat.cor = cor;

            set(ddat.axial_yline{i},'Xdata',[sag sag]);
            set(ddat.axial_xline{i},'Ydata',[cor cor]);
            set(ddat.coronal_yline{i},'Xdata',[sag sag]);
            set(ddat.sagittal_yline{i},'Xdata',[cor cor]); 

            for cl = 1:3
                Scor(:, :, cl) = squeeze(ddat.combinedImg{i}(cl).combound(:, cor,:))';
                Ssag(:, :, cl) = squeeze(ddat.combinedImg{i}(cl).combound(sag, :, :))';
            end
            
            set(ddat.coronal_image{i},'CData',Scor);
            set(ddat.sagittal_image{i},'CData',Ssag);
            set(findobj('Tag', 'SagSlider'),'Value',sag);
            set(findobj('Tag', 'CorSlider'),'Value',cor);
            if coordindex == 1,
                set(findobj('Tag', 'crosshairPos'), 'String', sprintf('%7.0d %7.0d %7.0d',sag, cor, ddat.axi));
            elseif coordindex == 2,
                sag_mm = (sag-ddat.origin(1))*ddat.pixdim(1);
                cor_mm = (cor-ddat.origin(2))*ddat.pixdim(2);    
                set(findobj('Tag', 'crosshairPos'), 'String', sprintf('%7.1f %7.1f %7.1f',sag_mm, cor_mm, ddat.axi));
            end 
end

if pp == 0
    set(findobj('Tag', 'crosshairVal'),'String',...
        sprintf('PP: %4.2f',...%mri.img(ddat.sag, ddat.cor, ddat.axi),...
                ddat.region_ppmap(ddat.sag, ddat.cor, ddat.axi)));
elseif pp == 1;
    set(findobj('Tag', 'crosshairVal'),'String',...
        sprintf('PP: %4.2f',...%mri.img(ddat.sag, ddat.cor, ddat.axi),...
                ddat.voxel_ppmap(ddat.sag, ddat.cor, ddat.axi)));
end

valId = cell2mat(ddat.total_region_name(:, 1));
curIdVal = ddat.region_struct.img(ddat.sag, ddat.cor, ddat.axi);

curIdPos = find(ismember(valId, curIdVal)); 

% Check to be sure that position is updated properly
    if get(findobj('Tag', 'viewZScores'), 'Value') == 0
    set(findobj('Tag', 'crosshairVal'),'String',...
        sprintf('Value at Voxel: %4.2f',...%ddat.mri_struct.img(ddat.sag, ddat.cor, ddat.axi),
                    ddat.img{i}(ddat.sag, ddat.cor, ddat.axi)));
    elseif get(findobj('Tag', 'viewZScores'), 'Value') == 1
        set(findobj('Tag', 'crosshairVal'),'String',...
        sprintf('Z = %4.2f',...%ddat.mri_struct.img(ddat.sag, ddat.cor, ddat.axi),
                    ddat.img{i}(ddat.sag, ddat.cor, ddat.axi)));
    end

% josh blocked this off xxx, re add if re add
if curIdPos
    set(findobj('Tag', 'curInfo'), 'ForegroundColor','g',...
        'FontSize', 10, 'HorizontalAlignment', 'left', 'String', ...
         sprintf('Current crosshair is located in the region: %s', ...
                  ddat.total_region_name{curIdPos, 2}));
else
    set(findobj('Tag', 'curInfo'), 'String', '');
end

talsag = (ddat.sag-ddat.origin(1))*ddat.pixdim(1);
talcor = (ddat.cor-ddat.origin(2))*ddat.pixdim(2);
talaxi = (ddat.axi-ddat.origin(3))*ddat.pixdim(3); 
outpoints = round(mni2tal([talsag, talcor, talaxi]));

%set(ddat.talairachText,'ForegroundColor','r',...
 %   'FontSize', 10, 'HorizontalAlignment', 'left', 'String',sprintf('Current Crosshair Talairach Coordinates: %s, %s, %s',num2str(outpoints(1)), num2str(outpoints(2)), num2str(outpoints(3))));

 
end


        if ddat.nCompare == 1
            if get(findobj('Tag', 'viewZScores'), 'Value') == 0
            set(findobj('Tag', 'crosshairVal1'),'String',...
                sprintf('Value at Voxel: %4.2f', ddat.img{1}(ddat.sag, ddat.cor, ddat.axi)));
            elseif get(findobj('Tag', 'viewZScores'), 'Value') == 1
                set(findobj('Tag', 'crosshairVal1'),'String',...
                sprintf('Z = %4.2f', ddat.img{1}(ddat.sag, ddat.cor, ddat.axi)));
            end
        else
            if get(findobj('Tag', 'viewZScores'), 'Value') == 0
                for iPop = 1:ddat.nCompare
                    set(findobj('Tag', ['crosshairVal' num2str(iPop)]),'String',...
                        sprintf('Value at Voxel: %4.2f',...
                        ddat.img{iPop}(ddat.sag, ddat.cor, ddat.axi)));
                end
            else
                for iPop = 1:ddat.nCompare
                    set(findobj('Tag', ['crosshairVal' num2str(iPop)]),'String',...
                        sprintf('Z = %4.2f',...
                        ddat.img{iPop}(ddat.sag, ddat.cor, ddat.axi)));
                end
            end
        end


guidata(gcbf, ddat);
return