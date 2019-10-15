% MAIN SCRIPT 
% Modification test
clear all; %#ok<CLALL>

addpath('/Users/maxime/Documents/ORBs project/Dynamic model/v122/cprintf');

disp('                                                                  ');
disp('------------------------------------------------------------------');
cprintf('comment', '                    WELCOME TO NavTOP v1.0        \n');
disp('                                                                  ');
disp('Navlab Tool for Orbit Propagation');
disp('Copyright (C) Maxime Largeaud');
disp('Code written by Maxime Largeaud');
disp('                                                                  ');
disp('------------------------------------------------------------------');
disp('                                                                  ');

%% --------------------------- READING INPUTS --------------------------- %

% Prompting user to give an input file
promptInput = 'Please provide an input file :';
file = input(promptInput);

% reading input data from file
disp('Reading inputs ...');
[case_name, mb, sc_mass, sc_refl_area, sc_refl_coeff, rasc, dec, stime, ...
    mb_init_TrueAnom, sc_init_OrbElems, max_deg, time_step, final_time] = ...
    file_reading_functions.read_inputs(file);

% reading main gravitating body data
disp('Reading main gravitating body data ...');
bodyPath = strcat('bodies_data/', mb, '/phys_orb_param.txt');
spHarmPath = strcat('bodies_data/', mb, '/sph_harm.txt');
spHarm = load(spHarmPath); 
[mu, r0, RotVecBF, InMomt, PA2BF, OrbElemsMB] = ...
    file_reading_functions.read_body_phys_prop(bodyPath);


%% ---------------------------- INITIALIZING ---------------------------- %
disp('                                                                  ');
disp('Initializing ...');

% generating main body object
disp('Generating main body object ...');
mb_obj = main_body(mu, r0, RotVecBF, InMomt, PA2BF, rasc, dec, stime, ...
    OrbElemsMB);

% generating spacecraft
disp('Generating spacecraft object ...');
sc_obj = spacecraft(mb_obj, sc_mass, sc_refl_area, sc_refl_coeff, ...
    sc_init_OrbElems);


%% ---------------------------- PROPAGATING ----------------------------- %
disp('                                                                  ');
disp('Propagating ...');

% propagating orbit
[sc_obj, prop_res] = spacecraft.propagate_orbit_orb_elems_RK4 ...
    (sc_obj, 0, final_time, time_step, spHarm, max_deg);
disp('Propagation complete.');


%% -------------------------- PLOTTING RESULTS -------------------------- %
disp('                                                                  ');
disp('Plotting results ...');

% plotting orbits on map
pos = [];
pos(1, :, :) = [prop_res(:, 2), prop_res(:, 3)*180/pi, ...
    prop_res(:, 4)*180/pi];

fin_orb_index = size(prop_res, 1);
plotting_functions.plot_wdata(pos(1,:,2)', pos(1,:,3)', ...
    prop_res(fin_orb_index,5), prop_res(fin_orb_index,6), ...
    prop_res(fin_orb_index,7), prop_res(fin_orb_index, 18:23)', 'Hours', ...
    prop_res(:, 1)/3600, ...
    prop_res(:, 17), '', '~');

orbs_fig = gcf;

% plotting orbital elements and perturbating accelerations
[prop_fig, acc_fig] = plotting_functions.plotPropRes(prop_res(:, 1)/3600, ...
    prop_res(:, 18:23), prop_res(:, 8:16));


%% ------------------------- EXPORTING RESULTS -------------------------- %

expIndex = 'd';
while (expIndex ~= 'y' || expIndex ~= 'n')

    % asking user if results need to be exported
    disp(' ');
    promptExport = 'Export results (''y'' or ''n'') ?';
    expIndex = input(promptExport);

    switch(expIndex)

        case 'y'
            
            addpath ...
            ('/Users/maxime/Documents/ORBs project/Dynamic model/v122/export_fig');
        
            addpath ...
            ('/Users/maxime/Documents/ORBs project/Dynamic model/v122'); 
            
            % test if the results folder does not exist, 
            % create one if it is the case 
            if exist('results', 'dir') == 0
                mkdir results;
            end
            
            % test if case_name folder already exists
            if exist(strcat('results/', case_name), 'dir') == 7
                disp(' ');
                promptOW = (strcat(case_name, ...
                ' result folder already exists, overwrite (''y'' or ''n'') ?'));
                owIndex = input(promptOW); 
                
                switch(owIndex)
                    
                    % user wants to overwrite existing folder
                    case 'y'
                        disp(' ');
                        disp('Exporting results ...');
                        
                        % moving to results directory 
                        cd(strcat('results/', case_name));
                        
                        % copying input file into results folder 
                        copyfile(strcat('../../', file));
                        
                        % export propagation data table
                        file_writing_functions.PropData2File(case_name, prop_res);
                                              
                        % export orbits plot
                        set(0, 'CurrentFigure', orbs_fig)
                        export_fig(strcat('orbs_', case_name));

                        % export orbital elements and perturbating accel. plots
                        set(0, 'CurrentFigure', prop_fig)
                        export_fig(strcat('orbElems_', case_name));
                        set(0, 'CurrentFigure', acc_fig)
                        export_fig(strcat('pertAcc_', case_name));

                        cd ../..
                        disp(' ');
                        disp(strcat('Results were saved in /results/', case_name, '.'));
                        disp('End of program.');
                        return
                     
                    % no overwrite, autosave   
                    case 'n'
                        
                        cd results
                        
                        % test if the autosave folder does not exist, 
                        % create one if it is the case
                        if exist('autosave', 'dir') == 0
                            mkdir autosave;
                        end
                        
                        cd autosave;
                        
                        disp(' ');
                        disp('Exporting results ...');
                        
                        % copying input file into results folder 
                        copyfile(strcat('../../', file));
                        
                        % export propagation data table
                        file_writing_functions.PropData2File(case_name, prop_res);
                        
                        % export orbits plot
                        set(0, 'CurrentFigure', orbs_fig)
                        export_fig(strcat('orbs_', case_name));

                        % export orbital elements and perturbating accel. plots
                        set(0, 'CurrentFigure', prop_fig)
                        export_fig(strcat('orbElems_', case_name));
                        set(0, 'CurrentFigure', acc_fig)
                        export_fig(strcat('pertAcc_', case_name));

                        cd ../.. 
                        disp(' ');
                        disp('Results were saved in /results/autosave.');
                        disp('Please choose a different case name.');
                        disp('End of program.');
                        return
                        
                    % default case, autosave
                    otherwise
                        
                        cd results
                        
                        % test if the autosave folder does not exist, 
                        % create one otherwise
                        if exist('autosave', 'dir') == 0                            
                            mkdir autosave;
                        end
                        
                        cd autosave
                        
                        disp(' ');
                        disp('Exporting results ...');
                        
                        % copying input file into results folder 
                        copyfile(strcat('../../', file));
                        
                        % export propagation data table
                        file_writing_functions.PropData2File(case_name, prop_res);
                        
                        % export orbits plot
                        set(0, 'CurrentFigure', orbs_fig)
                        export_fig(strcat('orbs_', case_name));

                        % export orbital elements and perturbating accel. plots
                        set(0, 'CurrentFigure', prop_fig)
                        export_fig(strcat('orbElems_', case_name));
                        set(0, 'CurrentFigure', acc_fig)
                        export_fig(strcat('pertAcc_', case_name));

                        cd ../.. 
                        disp(' ');
                        warning ...
                        ('No export option could be applied, results were saved in /results/autosave.');
                        disp('End of program.');
                        return
                    
                end 
             
            % if the case_name folder does not exist, save    
            else 
                
                cd results;
                mkdir(case_name);
                cd(case_name);

                disp(' ');
                disp('Exporting results ...');
                
                % copying input file into results folder 
                copyfile(strcat('../../', file));
                
                % export propagation data table
                file_writing_functions.PropData2File(case_name, prop_res);
                
                % export orbits plot
                set(0, 'CurrentFigure', orbs_fig)
                export_fig(strcat('orbs_', case_name));
                
                % export orbital elements and perturbating accel. plots
                set(0, 'CurrentFigure', prop_fig)
                export_fig(strcat('orbElems_', case_name));
                set(0, 'CurrentFigure', acc_fig)
                export_fig(strcat('pertAcc_', case_name));
                
                cd ../..
                disp(' ');
                disp(strcat('Results were saved in /results/', case_name, '.'));
                disp('End of program.');
                return
                
            end 

        case 'n'
            disp(' ');
            disp('No export.');
            disp('End of program.');
            return
        
        otherwise
            disp(' ');
            disp('The program did not understand your entry, please try again.');

    end 
    
end

% if no export options were understood by the program, autosave
addpath ...
('/Users/maxime/Documents/ORBs project/Dynamic model/v122/export_fig');

% test if the results and autosave folders do not exist, create them 
if exist('results', 'dir') == 0

    mkdir results;
    cd results
    mkdir autosave;

else 

    cd results 

    % test if the autosave folder does not exist, create one
    if exist('autosave', 'dir') == 0
        mkdir autosave;
    end

end

cd autosave;

disp(' ');
disp('Exporting results ...');

% copying input file into results folder 
copyfile(strcat('../../', file));

% export propagation data table
file_writing_functions.PropData2File(case_name, prop_res);

% export orbits plot
set(0, 'CurrentFigure', orbs_fig)
export_fig(strcat('orbs_', case_name));

% export orbital elements and perturbating accel. plots
set(0, 'CurrentFigure', prop_fig)
export_fig(strcat('orbElems_', case_name));
set(0, 'CurrentFigure', acc_fig)
export_fig(strcat('pertAcc_', case_name));

cd ../.. 
disp(' ');
warning('No export option could be applied, results were saved in results/autosave.');
disp('End of program.');

% END OF MAIN PROGRAM
