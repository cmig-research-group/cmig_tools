function save_showVol_images(coords, outdir, fstem)

    if isempty(coords) 
        coords = [99 99 129]; 
    end
    if ~exist(outdir, 'dir') 
        mkdir(outdir); 
    end
    if ~exist('fstem', 'var') 
        fstem = 'ScreenShot'; 
    end
    
    % get the hidden data of the showVol window
    showVolWins = findall(groot, 'Type', 'figure', 'Tag', 'showVol');
    
    myWin = showVolWins(1); %just take the first one
    
    handles = guidata(myWin);
    
    
    %orientations = {'coronal' 'sagittal' 'axial'};

    for o = 1%:3 doesnt work
        % handles.ORIENTATION = o;
        % handles = displayNewOrientation(handles);
        % guidata(myWin, handles);
        % updateDisplay_newvol(handles, true);
        % drawnow;
        for v = 1:length(handles.vols)
            % --- Set current volume at the start of the v loop ---
            set(handles.currentVolPB, 'ForegroundColor', 'k');
            tmp = find(handles.hideVol == 0);
            handles.currentVol = tmp(v);
            handles.currentVolPB = handles.volPB(handles.currentVol);
            set(handles.currentVolPB, 'ForegroundColor', 'r');
            handles.Mvxl2lph = handles.vols{handles.currentVol}.Mvxl2lph;
            updateDisplay_newvol(handles, true);
            guidata(myWin, handles);
            drawnow;
       
            % -----------------------------------------------------

            for ci = 1:size(coords,1)
                coord = coords(ci,:);
                rcs = [coord 1];
                M = handles.Mvxl2lph_atlas;
                lph = M * rcs(:);
                clear ud
                ud.lph = lph(1:3);
                ud.M = M;
                ud.sender = myWin;
                set(myWin, 'UserData',ud);

                % hide crosshair
                set(handles.toggleLines1, 'ForegroundColor', 'k');
                set(handles.axes1_X, 'Visible', 'off');

                % hide ROI overlays
                set(handles.anat.aseg.roi_text,'visible','off')
                set(handles.anat.fiber.roi_text,'visible','off')
                set(handles.anat.aparcaseg.roi_text,'visible','off')

                % 2) save screenshot
                drawnow;
                shot = getframe(handles.axes1);
                handles.screenshotnum = v;
                %fname_out = fullfile(outdir, sprintf('%s_orient%d_vol%04d_coord%d.png', fstem, o, v, ci));
                fname_out = fullfile(outdir, sprintf('%s_vol%04d_coord%d.png', fstem, v, ci)); % no point including orientation if 
                imwrite(shot.cdata,fname_out); 
                fprintf(1,'File %s written.\n',fname_out);
            end
        end
    end 
end

