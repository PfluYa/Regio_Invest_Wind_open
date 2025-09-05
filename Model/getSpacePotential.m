function outputBaseTable = getSpacePotential(inputBaseTable)
%GETSPACEPOTENTUAL Add Data for Free Area in each Nuts3 Region for Solar Open Space and WindOnshore to the Basetable
%   Input: BaseTableRegioInvest -> NutsID's in the first Row are relevant. Other Data is added automatically from
%   DataGeneral/Geographical/ 
%   Filenames: areaSolarOpenSpace.xlsx | areaWindOnshore.xlsx

%% Init BaseTableRegioInvest and add SpaceData
 % Filenames and relevant Variables can be found in local Functions!!!! 
localBaseTable = inputBaseTable;

%% Get you Data
% Add SolarSpaceData 
%localBaseTable = getSolarSpacePotential(localBaseTable);
% Add WindOnshoreData
localBaseTable = getWindSpacePotential(localBaseTable);

%% Write Output
outputBaseTable = localBaseTable;

end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%      Local Functions      %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%% FUNCTION GETSOLARSPACEPOTENTIAL
function outputSolarSpace = getSolarSpacePotential(inputBaseTable)
% Init Options, Filename and Relevant Data
filenameData  = 'areaSolarOpenSpace.xlsx';
relevantVariables = 'relativeAvailableSolarSpace';
localTable = inputBaseTable;

% Read Data Table
localSolarSpace = readtable(filenameData);
% Merge Data with Basetable
localTable = outerjoin(localTable,localSolarSpace,'Type','left','Keys','nutsID','RightVariables',relevantVariables);
% Write Output
outputSolarSpace = localTable;
end
%% FUNCTION GETWINDONSHORESPACEPOTENTIAL
function outputWindSpace = getWindSpacePotential(inputTable)
% Init Options, Filename and Relevant Data
filenameData     = 'areaWindOnshore.xlsx';
relevantVariables = 'relativeAvailableWindSpace';
localTable = inputTable;

% Read Data Table
localWindSpace = readtable(filenameData);
% Merge Data with Basetable
localTable = outerjoin(localTable,localWindSpace,'Type','left','Keys','nutsID','RightVariables',relevantVariables);
% Write Output
outputWindSpace = localTable;
end
