%% BMEN E3810: Biomedical Engineering Laboratory I (INSTRON Data Analysis)
%  Written by: Dr. Lauren N. Heckelman, Jeannine Turgeon, and Anjali Parande
%  Fall 2023

%% Initialize the Workspace:
clear; clc; close all;

%% Define Data Parameters:
materials = {'Silicone_Thin', 'Silicone_Thick', ...
             'Skin_Raw', 'Skin_Treated'};
teamnames = {'Wed01',  'Wed02',  'Wed03',  'Wed04', ...
             'Wed05',  'Wed06', 'Wed07', ...
             'Thurs01','Thurs02','Thurs03', 'Thurs04', ...
             'Thurs05','Thurs06','Thurs07'};

filename = 'III.2 - Biomechanical Testing - Class Data.xlsx';

%% Import Raw Data:
for k = 1:length(materials)
    rawdata.(materials{k}) = readcell(filename, 'Sheet', materials{k});
    rawdata.(materials{k})(cellfun(@(x) any(ismissing(x)), rawdata.(materials{k}))) = {NaN};
end

%% Organize Data:
for k = 1:length(materials)
    for n = 1:length(teamnames)
        data.(materials{k}).L0(1,n)   = (rawdata.(materials{k})(3,4*n-1));
        data.(materials{k}).Width(1,n) = (rawdata.(materials{k})(4,4*n-1));
        data.(materials{k}).Thickness(1,n)  = (rawdata.(materials{k})(5,4*n-1));
        data.(materials{k}).Area(1,n)   = (rawdata.(materials{k})(6,4*n-1));
        
        data.(materials{k}).Pos{1,n}    = (rawdata.(materials{k})(9:end,4*n-2));
        data.(materials{k}).Force{1,n}  = (rawdata.(materials{k})(9:end,4*n-1));
    end
end

%% DISPLAY DATA: Load vs. Extension
% yValues = [];
for k = 1:length(materials)

    yValues = [];  % empty array to place position values into
     for n = 1:length(teamnames)
         Force = cell2mat(data.(materials{k}).Force{1,n});
         Pos  = cell2mat(data.(materials{k}).Pos{1,n});
         
         yValues = [yValues,Pos];
     
         if k == 1
             figure (1);
             hold on
             plot(Pos,Force,'Color',[0,0,1,0.3]); % plot all the data across n with faint line opacity 
         elseif k == 2
             plot(Pos,Force,'Color',[1,0,0,0.3]);
             title("Load vs Extension for Silicone", 'Interpreter', 'None');
             xlabel('Extension (mm)');
             ylabel('Load (N)')
         elseif k == 3
             figure(2);
             hold on
             plot(Pos, Force,'Color',[0,0,1,0.3]);
         elseif k == 4
             plot (Pos,Force, 'Color',[1,0,0,0.3]);
             title("Load vs Extension for Chicken Skin", 'Interpreter', 'None');
             xlabel('Extension (mm)');
             ylabel('Load (N)')
         end
     end
     % plot average position lines for better visualization
     averagePos = nanmean(yValues,2);
     plot(averagePos,Force,'LineWidth',2.5);
     % best way I knew how to create a legend
     if k == 2
         legend('Thin','','','','','','','','','','','','','', ...
             'Thin Average','','','','','','','','','','','','','', ...
             'Thick','Thick Average','Location','Best');
         xlim([0,90]);
     end
     if k == 4
         legend('Untreated','','','','','','','','','','','','','', ...
             'Untreated Average','','','','','','','','','','','','','', ...
             'Treated','Treated Average','Location','Best');
         xlim([0,80]);
     end
end
%% Useful Commands:
    % nanmean(cell2mat(data.Skin_Raw.Height)) --> Calculates the mean when
    % there are NaN entries in a matrix

    % gradient() --> computes the derivative of a dataset - use this for
    % linear region
   
%% DETERMINING VARIABLES AND DISPLAY DATA: Stress vs Strain
Stiffness = zeros(length(teamnames),length(materials));    % Matrix to hold the all the possible values of Stiffness
E = zeros(length(teamnames),length(materials));     % "..." Young's Modulus
uts = zeros(length(teamnames),length(materials));   % "..." ultimate tensile strength
extensibility = zeros(length(teamnames),length(materials)); % "..." extensiblilty

for k = 1:length(materials)
    Stress_values = []; % empty array to place stress values into
    Strain_values = []; % empty array to place strain values into
    max_length = 0; % intializing variable
    best_n = 0; % intializing variable
    for n = 1:length(teamnames)
        Force = cell2mat(data.(materials{k}).Force{1,n});
        Pos  = cell2mat(data.(materials{k}).Pos{1,n});

        % Zero force/position data
        Force = Force - Force(1);
        Pos = Pos - Pos(1);
        
        % Determine if and where NaN entries exist in force/position data
        idx = find(isnan(Force), 1, 'first');
        if isempty(idx) == 0            % safegaurd for if there are no NaN entries
            Force = Force(1:idx-1);
            Pos = Pos(1:idx-1);
        end
        % Crop off all values outside of linear range 
        if k == 1 | k == 2
            Force_linear = Force(:);    % all of the data is within linear region
            Pos_linear = Pos(:);
        elseif k == 3 | k == 4
            dy_dx = gradient(Force, Pos);
            gradient_threshold = 3;     % gradient threshold determined through testing
            linear_region = abs(dy_dx) < gradient_threshold;
            linear_region2 = diff(linear_region);
            idx2 = find(linear_region2 == -1, 1, 'first'); % finds the first -1 value
            if nnz(linear_region) == 0   %if there are no zero values (all linear)
                Force_linear = Force(:);
                Pos_linear = Pos(:);
            elseif idx2 == 1
                idx3 = find((linear_region2) == -1, 2, 'first');
                idx3 = idx3(2);     % finds the second occurance of -1 values if the first -1 value is the first value ever
                Force_linear = Force(1:idx3);
                Pos_linear = Force(1:idx3);
            else
                Force_linear = Force(1:idx2);
                Pos_linear = Pos(1:idx2);
            end
        end

        % Calculate linear regression of force/position data to determine
        % stiffness
            p = polyfit(Pos_linear, Force_linear, 1);
            Stiffness(n,k) = p(1)*1000; % to ensure N/m units
            
        % Stress = F/A
            area = cell2mat(data.(materials{k}).Area(1,n));
            Stress = Force./area;
        % Strain
            initial_length = cell2mat(data.(materials{k}).L0(1,n));
            Strain = Pos./initial_length;
        
        % Young's Modulus
            Stress_linear = Force_linear./area;
            Strain_linear = Pos_linear./initial_length;
            raw_E = Stress_linear./Strain_linear;
            inf_E = isinf(raw_E);      % identifies any infinite entries
            raw_E(inf_E) = NaN;        % makes infinite entries NaN
            nan_E = isnan(raw_E);      % identifies NaN entries
            add_E = mean(raw_E(~nan_E));  % removes NaN entries         
            E(n,k) = add_E;

        % Ultimate Tensile Strength
            uts(n,k) = max(Stress);

        % Extensibilty - end of strain on linear region
            extensibility(n,k) = max(Strain_linear)*100; % to make percentage

        if k == 1
            figure(3);
            hold on
            plot(Strain_linear,Stress_linear,'Color',[0,0,1,0.3]); % plot all the data across n with faint line opacity
        elseif k == 2
            plot(Strain_linear,Stress_linear,'Color',[1,0,0,0.3]);
            title(strcat("Stress vs Strain for Silicone"), 'Interpreter', 'None');
            xlabel('Strain');
            ylabel('Stress (MPa)')
        elseif k == 3
            figure(4)
            hold on
            plot(Strain_linear,Stress_linear,'Color',[0,0,1,0.3]);
        elseif k == 4
            plot(Strain_linear,Stress_linear,'Color',[1,0,0,0.3]);
            title(strcat("Stress vs Strain for Chicken Skin"), 'Interpreter', 'None');
            xlabel('Strain');
            ylabel('Stress (MPa)')
        end
        % create Strain matrix to pull from, determine which n produces longest linear strain 
        current_lengthStrain = length(Strain_linear);
        Strain_values(1:current_lengthStrain,n) = Strain_linear;
        if length(Strain_linear) > max_length
            max_length = length(Strain_linear);
            best_n = n;
        end
        % create Stress matrix, make any zero values NaN
        current_length = length(Stress_linear);
        Stress_values(1:current_length,n) = Stress_linear;
        Stress_values(Stress_values == 0) = NaN;
    end
    % plot average linear stress lines for better visualization
        averageStress = nanmean(Stress_values,2);
    plot(Strain_values(:,best_n),averageStress,'LineWidth',2);
    % best way I knew how to create a legend
    if k == 2
       legend('Thin','','','','','','','','','','','','','', ...
          'Thin Average','','','','','','','','','','','','','', ...
          'Thick','Thick Average','Location','Best');
    end
    if k == 4
       legend('Untreated','','','','','','','','','','','','','', ...
           'Untreated Average','','','','','','','','','','','','','', ...
           'Treated','Treated Average','Location','Best');
    end
end

%% ANOVA
% Stiffness
    Stiff_data = [Stiffness(:,1),Stiffness(:,2),Stiffness(:,3),Stiffness(:,4)];
    pValue_Stiff = anova1(Stiff_data, [],'off'); % off stops table output
    if pValue_Stiff < 0.05
      [h, pStiff_Si, ci, stats] = ttest2(Stiffness(:,1),Stiffness(:,2));
      [h, pStiff_Ch, ci, stats] = ttest2(Stiffness(:,3),Stiffness(:,4));
    end
% Young's Modulus
    E_data = [E(:,1),E(:,2),E(:,3),E(:,4)];
    pValue_E = anova1(E_data, [],'off'); % off stops table output
    if pValue_E < 0.05
      [h, pE_Si, ci, stats] = ttest2(E(:,1),E(:,2));
      [h, pE_Ch, ci, stats] = ttest2(E(:,3),E(:,4));
    end
% Ultimate Tensile Strength
    UTS_data = [uts(:,3),uts(:,4)]; % no uts for silicone because no ruptures
    pValue_UTS = anova1(UTS_data, [],'off'); % off stops table output
    if pValue_UTS < 0.05
      [h, pUTS_Ch, ci, stats] = ttest2(uts(:,3),uts(:,4));
    end
% Extensibility
    Extens_data = [extensibility(:,3),extensibility(:,4)]; % no extensibility for silicone because no ruptures
    pValue_Extens = anova1(Extens_data, [],'off'); % off stops table output
    if pValue_Extens < 0.05
      % [h, pExtens_Si, ci, stats] = ttest2(extensibility(:,1),extensibility(:,2));
      [h, pExtens_Ch, ci, stats] = ttest2(extensibility(:,3),extensibility(:,4));
    end

%% Plotting

% Stiffness
meanStiff = mean(Stiffness);
stdvStiff = std(Stiffness);
figure(5);
bar(meanStiff);
hold on
errorbar(1:size(meanStiff,2),meanStiff,stdvStiff,'k.','LineWidth',0.5);
title('Mean Stiffness of Each Sample ');
xlabel('Sample Type');
ylabel('Stiffness (N/m)');
xticks(1:4);
xticklabels({'Silicone Thin', 'Silicone Thick', 'Chicken Untreated', 'Chicken Treated'});
% include significance if applicable
if pStiff_Si < 0.05
    sigstar({[1,2]},pStiff_Si);
end
if pStiff_Ch < 0.05
    sigstar({[3,4]},pStiff_Ch);
end
% display the mean and std values on the bar plot
for i = 1:length(meanStiff)
    text(i, meanStiff(i) + stdvStiff(i) + 0.1, sprintf('%0.1f ± %0.1f', meanStiff(i), stdvStiff(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% Young's Modulus
meanMod = mean(E);
stdvMod = std(E);
figure(6);
bar(meanMod);
hold on
errorbar(1:size(meanMod,2),meanMod, stdvMod, 'k.','LineWidth', 0.5);
title("Mean Young's Modulus of Each Sample ");
xlabel('Sample Type');
ylabel("Young's Modulus (MPa)");
xticks(1:4);
xticklabels({'Silicone Thin', 'Silicone Thick', 'Chicken Untreated', 'Chicken Treated'});
if pE_Si < 0.05
    sigstar({[1,2]},pE_Si);
end
if pE_Ch < 0.05
    sigstar({[3,4]},pE_Ch);
end
for i = 1:length(meanMod)
    text(i, meanMod(i) + stdvMod(i) + 0.1, sprintf('%0.1f ± %0.1f', meanMod(i), stdvMod(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% Ultimate Tensile Strength
meanUTS = mean(uts(:,3:4)); % doesn't make sense because UTS should be stronger for silicone
stdvUTS = std(uts(:,3:4));
figure (7);
bar(meanUTS);
hold on
errorbar(1:size(meanUTS,2),meanUTS,stdvUTS, 'k.','LineWidth',0.5);
title("Mean Ulitimate Tensile Strength of Each Sample ");
xlabel('Sample Type');
ylabel("Ultimate Tensile Strength (MPa)");
xticks(1:4);
xticklabels({'Silicone Thin', 'Silicone Thick', 'Chicken Untreated', 'Chicken Treated'});
if pUTS_Si < 0.05
    sigstar({[1,2]},pUTS_Si);
end
if pUTS_Ch < 0.05
    sigstar({[3,4]},pUTS_Ch);
end
for i = 1:length(meanUTS)
    text(i, meanUTS(i) + stdvUTS(i) + 0.1, sprintf('%0.1f ± %0.1f', meanUTS(i), stdvUTS(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% Extensibility
meanExtens = mean(extensibility(:,3:4));
stdvExtens = std(extensibility(:,3:4));
figure(8);
bar(meanExtens);
hold on
errorbar(1:size(meanExtens,2),meanExtens,stdvExtens,'k.','LineWidth',0.5);
title("Mean Extensiblity of the Chicken Samples ");
xlabel('Sample Type');
ylabel("Extensibility (%)");
xticks(1:2);
xticklabels({'Chicken Untreated', 'Chicken Treated'});
if pExtens_Ch < 0.05
    sigstar({[1,2]},pExtens_Ch);
end
for i = 1:length(meanExtens)
    text(i, meanExtens(i) + stdvExtens(i) + 0.1, sprintf('%0.1f ± %0.1f', meanExtens(i), stdvExtens(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end





