function outputGeoData = getGeoShapeData()
%% Function GETREGIODATA - loads your Shapefile for given GeoScope
%Input: none -> is given by global Variable paraRegioInvest
global paraRegioInvest
%Output: Table with relevant Informations,GeoID's, CountryCodes 

%% Init filenames of Shapefiles and which nuts needs to be
% generally excluded

excludeNutsRegion = {'FRM','FRY'}; % This should be Corsika and something in the Caribbean Sea
filenamePolyShape = 'nuts2021Poly.shp';
filenamePointShape = 'nuts2021Point.shp';

%% Init Global Options, Parameter and Filenames
geoScope = paraRegioInvest.geoScope; % This is only the GeoScope (which Countries?) 

%% Get Shapefile and Data of GeoScope
localGeoData = getShapefile(filenamePolyShape,filenamePointShape);

% Remove irrelevant Countries And Nuts not relevant Nuts Europe
localGeoData = removeCountries(localGeoData,geoScope);
localGeoData = removeNotInEurope(localGeoData);
localGeoData = handleDenmark(localGeoData,geoScope); % This Function handles Denmark West and Denmak East - The nuts regions which are not needed will be removed
localGeoData(contains(localGeoData.NUTS_ID,excludeNutsRegion),:) = [];

% Get Area and Population of each nuts and
localGeoData = [localGeoData, rowfun(@getAreaFromShape,localGeoData,'InputVariables',{'latPoly','lonPoly'},'OutputVariableNames','totalArea')];
% localGeoData = getPopulationForShape(localGeoData,populationYear);

% Change Columnnames, Delete Columns and sort the Table

localGeoData = renamevars(localGeoData,{'NUTS_ID','CNTR_CODE','NAME_LATN'},{'nutsID','countryCode','cityName'});
localGeoData(:,{'LEVL_CODE','MOUNT_TYPE','URBN_TYPE','COAST_TYPE','FID','NUTS_NAME'}) = [];

localGeoData = movevars(localGeoData,{'nutsID','countryCode','cityName'},'Before',{'Geometry'});
localGeoData = movevars(localGeoData, 'latPoly', 'Before', 'lonPoly');
outputGeoData = localGeoData;

%% Local Function GETSHAPEFILE
    function outputShape = getShapefile(inpFilenamePoly,inpFilenamePoint)
        filenamePoly = inpFilenamePoly;
        filenamePoint = inpFilenamePoint;
        
        polyShape = struct2table(shaperead(filenamePoly));
        pointShape = struct2table(shaperead(filenamePoint));
        
        polyShape = renamevars(polyShape,{'X','Y'},{'lonPoly','latPoly'});
        pointShape = renamevars(pointShape,{'X','Y'},{'lonPoint','latPoint'});
        
        outputShape = innerjoin(polyShape,pointShape(:,{'NUTS_ID','latPoint','lonPoint'}));
        
    end
%% Local Function REMOVECOUNTRIES
    function outputShape = removeCountries(inputShape,geoScope)
        localshape = inputShape;
        localGeoScope = cellfun(@(x) x(1:2),geoScope,'UniformOutput',false); % Because of DKW or DKE
        
        idxRelevantScope = contains(localshape.NUTS_ID,localGeoScope);
        localshape = localshape(idxRelevantScope,:);
        outputShape = localshape;
    end
%% Local Function GetAreaFromShape
    function areaOutput = getAreaFromShape(latPoly,lonPoly)
        
        localLatPoly = cell2mat(latPoly);
        localLonPoly = cell2mat(lonPoly);
        
        % area in km^2 (sum over all partial areas of that region)
        localArea = ...
            sum( 1/10^6 * areaint( localLatPoly, localLonPoly, referenceEllipsoid('World Geodetic System 1984'))) ;
        areaOutput = localArea;
        
    end
%% Local Function HANDLEDENMARK
    function outputShape = handleDenmark(inputShape,localGeoScope)
        localShape = inputShape;
        if any(contains(localGeoScope,'DK')) % Check if DK is considered
            if any(ismember(localGeoScope,'DK')) % Check if DK DKW and DKE was selected -> in this Case nothing happens and DK is written out
                %do nothing 
            elseif ~any(ismember(localGeoScope,'DKE'))&& any(ismember(localGeoScope,'DKW')) % Delete DKE - the common Case
                nutsDKE = {'DK011','DK012','DK013','DK014','DK021','DK022'};
                idxToDelete = ismember(localShape.NUTS_ID,nutsDKE);
                localShape(idxToDelete,:) = [];
                idxRenameCountryCode = contains(localShape.NUTS_ID,'DK');
                localShape.CNTR_CODE(idxRenameCountryCode) = {'DKW'};
            elseif ~any(ismember(localGeoScope,'DKW'))&& any(ismember(localGeoScope,'DKE')) % If this happens, ask CW if you really want to consider DKE and not DKW
                nutsDKW = {'DK031','DK032','DK041','DK042','DK050'};
                idxToDelete = ismember(localShape.NUTS_ID,nutsDKW);
                localShape(idxToDelete,:) = [];
                idxRenameCountryCode = contains(localShape.NUTS_ID,'DK');
                localShape.CNTR_CODE(idxRenameCountryCode) = {'DKE'};
            elseif any(ismember(localGeoScope,'DKW'))&& any(ismember(localGeoScope,'DKE'))  % Both is available but we want to separate them -> DKE and DKW
                nutsDKE = {'DK011','DK012','DK013','DK014','DK021','DK022'};
                idxToRenameDKE = ismember(localShape.NUTS_ID,nutsDKE);
                localShape.CNTR_CODE(idxToRenameDKE) = {'DKE'};
                nutsDKW = {'DK031','DK032','DK041','DK042','DK050'};
                idxToRenameDKW = ismember(localShape.NUTS_ID,nutsDKW);
                localShape.CNTR_CODE(idxToRenameDKW) = {'DKW'};
            end
        end
        outputShape = localShape;
    end
%% Local Function REMOVENOTINEUROPE
    function outputTable = removeNotInEurope(inputTable)
        localTable = inputTable;
        boundLat = [35 72];
        boundLon = [-13 40];
        
        idxNotInBoundLat = localTable.latPoint < min(boundLat) | localTable.latPoint > max(boundLat);
        idxNotInBoundLon = localTable.lonPoint < min(boundLon) | localTable.lonPoint > max(boundLon);
        
        idxNotInBound = idxNotInBoundLat|idxNotInBoundLon;
        
        localTable(idxNotInBound,:) = [];
        outputTable = localTable;
    end
end




