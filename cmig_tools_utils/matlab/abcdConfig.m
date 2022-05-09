function cfg = abcdConfig(arg)
      % abcdConfig  get configuration for ABCD analysis and visualization
      %
      %  This function helps you set up tools to work with ABCD imaging data:
      %
      %   showVol  anatomical volume viewer with ROIs and annotation        (requires showVolData)
      %   FEMA     calculate model statistics on diffusion and other images (requires abcd_sync)
      %
      %  cfg = abcdConfig
      %
      %    reurns config structure, with paths to data in cfg.data.*. If it is the first time it
      %     will create a default config file (~/abcdConfig.json) and prompt you to setup paths.
      %
      %
      %  cfg = abcdConfig(toolname)
      %
      %    tests that data needed by the tool are present
      %
      %
      % This software is Copyright (c) 2022 The Regents of the University of California. All Rights Reserved.
      % See LICENSE.
      
      % hidden usage to write the config file--this could be useful for a tool to automatically download and set up the users data/repos
      %  cfg = abcdConfig(cfg)
      %
      %    writes a cfg structure to the configuaration file at ~/abcdConfig.json. Could use this to add
      %     your own custom fields
      
      configFile = fullfile(userhome,'abcdConfig.json');
      
      if nargin < 1
        arg = '';
      end
      
      %handle case of writing a cfg structure directly
      if isstruct(arg)
        disp('Writing abcdConfig.json')
        cfg=arg;
        if ~isfield(cfg,'data') || ~isfield(cfg.data,'showVolData') || ~isfield(cfg.data,'abcd_sync')
          error('Your cfg structure is incomplete (needs data.showVolData and data.abcd_sync at minimum).')
        end
        writejson(configFile, cfg)
        return
      end
      
      
      %defaults
      default.showVolData = '/edit/path/to/showVolData';
      default.abcd_sync = '/edit/path/to/abcd-sync';
      
      % data resource paths needed for each tool
      paths.FEMA = {'abcd_sync'};
      paths.showVol = {'showVolData'};
      
      % read/initialize the config file
      if exist(configFile,'file') 
        cfg = jsondecode(fileread(configFile));
      else
        cfg = initCFG(default, paths); %write config with placeholders
        writejson(configFile, cfg)
        try
          edit(configFile)
        catch
          disp('Could not open abcdConfig.json for editing')
        end
        error(['%s: We have created the abcdConfig.json configuration file in your home directory. ' ...
          'Please edit it now, adding the paths to abcd_sync and showVolData as needed.'],mfilename)
      end
      
      %if toolname was specified, make sure the needed data exists
      if nargin
        if ~isfield(paths,arg), error('Invalid tool name (%s)', arg), end
        if ~exist(cfg.data.(paths.(arg){1}),'dir')
          try
            edit(configFile)
          catch
            disp('Could not open abcdConfig.json for editing')
          end
          error('%s: Data required by tool %s does not exist. Make sure you have properly set up ~/abcdConfig.json to point to %s.', mfilename, arg, paths.(arg){1})
        end
      end
      
      % ------------------------------------------------------------
      
      function cfg = initCFG(default, paths)
      allTools = fieldnames(paths);
      for iT = 1:length(allTools)
        tool = allTools{iT};
        
        for iD = 1:length(paths.(tool))
          dataName = paths.(tool){iD};
          cfg.data.(dataName) = default.(dataName);
        end
      end
        
      % ------------------------------------------------------------
      
      function writejson(fname, cfg)
      text = jsonencode(cfg);
      text = ['{\n' text(2:end)]; %make it cheaply human readable
      text = strrep(text,'},', '_XX_'); 
      text = strrep(text,':{', '_YY_');
      text = strrep(text,':', ':\n\t\t');
      text = strrep(text,',', ',\n\t');
      text = strrep(text,'}}', '}\n}');
      text = strrep(text,'_XX_', '},\n');
      text = strrep(text,'_YY_', ':\n\t{');
      text = sprintf(text);
      
      fid = fopen(fname,'w');
      fprintf(fid,'%s\n',text);
      fclose(fid);
      