%% Automatically generated BEARX Toolbox script
%
% This script was generated based on the user input from the BEARX Toolbox
% Graphical User Interface. Feel free to edit and adapt it further to your
% needs.
%
% Generated 11-Nov-2025 20:31:26
%


%% Clear workspace

% Clear all variables
clear

% Close all figures
close all

% Rehash Matlab search path
rehash path

% Import the correct module
import base.*


%% Define percentile function for summarizing the results

% User choice of percentiles
percentiles = [10 50 90];

% Create a percentiles function used to condense and report some results
percentilesFunc = @(x) prctile(x, percentiles, 2);

% Create a median function used to condense and report some results
medianFunc = @(x) median(x, 2);

% Create a legend for the percentiles
percentilesLegend = compose("%d%%", percentiles);


%% Prepare the output folder

outputFolder = fullfile(".", "output");
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end


%% Prepare an empty array of dummies

dummyObjects = {};


%% Prepare meta information 

% Create a meta information object
meta = Meta( ...
    EndogenousNames=["GDP", "CPI", "STN"] ...
    , ExogenousNames="Oil" ...
    , ShockNames=["DEM", "SUP", "POL"] ...
    , Order=2 ...
    , Intercept=true ...
    , EstimationSpan=datex.span("1977-Q1", "2014-Q4") ...
    , IdentificationHorizon=4 ...
);


%% Load input data table 

% Load the input data table
inputTbl = tablex.fromFile("/Users/myself/Documents/ogr-external-projects/ecb-bear/gui_poc/exampleDataX.csv");
display(inputTbl);


%% Create DataHolder object 

dataHolderObject = DataHolder(meta, inputTbl);


%% Prepare reduced-form estimator 

% Create a reduced-form estimator object
estimatorObject = estimator.NormalWishart( ...
    meta ...
    , Sigma="ar" ...
    , Burnin=0 ...
    , StabilityThreshold=NaN ...
    , MaxNumUnstableAttempts=1000 ...
    , Exogenous=false ...
    , BlockExogenous=false ...
    , Autoregression=0.8 ...
    , Lambda1=0.1 ...
    , Lambda2=0.5 ...
    , Lambda3=1 ...
    , Lambda4=100 ...
    , Lambda5=0.001 ...
);


%% Create reduced-form model 

% Assemble a reduced-form model from the components
redModel = ReducedForm( ...
    Meta=meta ...
    , DataHolder=dataHolderObject ...
    , Estimator=estimatorObject ...
    , Dummies=dummyObjects ...
);

display(redModel);


%% Initialize and presample the reduced-form model 

redModel.initialize();
info = redModel.presample(1000);
display(info);


%% Create Cholesky identification object 

ident = identifier.Cholesky( ...
    Order="" ...
);

display(ident);


%% Create a structural model 

structModel = Structural( ...
    reducedForm=redModel ...
    , identifier=ident ...
);
display(structModel);


%% Initialize and presample the structural model 

structModel.initialize();
info = structModel.presample(1000);
display(info);

% save(fullfile(outputFolder, "structuralModel.mat"), "structModel");


%% Run unconditional forecast using reduced-form model

% Run an unconditional forecast using the structural model
redForecastTbl = redModel.forecast( ...
    datex.span("2015-Q1", "2016-Q4") ...
    , stochasticResiduals=true ...
    , includeInitial=true ...
);

% Condense the forecast to percentiles
redForecastPercentilesTbl = tablex.apply(redForecastTbl, percentilesFunc);
display(redForecastPercentilesTbl);

% Define the output path for saving the results
outputPath = fullfile(outputFolder, "redForecastPercentiles");

% Save the forecast results in percentiles as MAT and/or CSV and/or
% and XLS files
% save(outputPath + ".mat", "redForecastPercentilesTbl");
tablex.writetimetable(redForecastPercentilesTbl, outputPath + ".csv");
tablex.writetimetable(redForecastPercentilesTbl, outputPath + ".xlsx");

% Plot the forecast results in percentiles
figureHandles = chartpack.forecastPercentiles( ...
    redForecastPercentilesTbl, redModel ...
    , "figureTitle", "Reduced-form model forecast percentiles" ...
    , "figureLegend", percentilesLegend ...
);

% Save the figures as a PDF
chartpack.printFiguresPDF(figureHandles, outputPath);


%% Run unconditional forecast using structural model 

% Run an unconditional forecast using the structural model
[structForecastTbl, structForecastContribsTbl] = structModel.forecast( ...
    datex.span("2015-Q1", "2016-Q4") ...
    , stochasticResiduals=true ...
    , includeInitial=true...
    , contributions=false ...
);

% Condense the forecast to percentiles
structForecastPercentilesTbl = tablex.apply(structForecastTbl, percentilesFunc);
display(structForecastPercentilesTbl);

% Define the output path for saving the results
outputPath = fullfile(outputFolder, "structForecastPercentiles");

% Save the forecast results in percentiles as MAT and/or CSV and/or
% and XLS files
% save(outputPath + ".mat", "structForecastPercentilesTbl");
tablex.writetimetable(structForecastPercentilesTbl, outputPath + ".csv");
tablex.writetimetable(structForecastPercentilesTbl, outputPath + ".xlsx");

if ~isempty(structForecastContribsTbl)
    % Condense the forecast contributions to percentiles
    structForecastContribsPercentilesTbl = tablex.apply(structForecastContribsTbl, percentilesFunc);

    % Flatten the 3D contributions table to 2D contributions table
    structForecastContribsPercentilesTbl = tablex.flatten(structForecastContribsPercentilesTbl);

    % Define the output path for saving the results
    outputPath = fullfile(outputFolder, "structForecastContribsPercentiles");

    % Save the percentiles of the forecast contributions
    % save(outputPath + ".mat", "structForecastContribsPercentilesTbl");
    tablex.writetimetable(structForecastContribsPercentilesTbl, outputPath + ".csv");
    tablex.writetimetable(structForecastContribsPercentilesTbl, outputPath + ".xlsx");
end

% Plot the model forecast results in percentiles
figureHandles = chartpack.forecastPercentiles( ...
    structForecastPercentilesTbl, structModel ...
    , "figureTitle", "Structural model forecast percentiles" ...
    , "figureLegend", percentilesLegend ...
);

% Save the figures
chartpack.printFiguresPDF(figureHandles, outputPath);


%% Simulate shock responses

% Simulate the responses to shocks over the shock response horizon
responseTbl = structModel.simulateResponses( ...
    includeInitial=true ...
);

% Condense the results to percentiles and flatten the 3D table to 2D table
responsePercentilesTbl = tablex.apply(responseTbl, percentilesFunc);
responsePercentilesTbl = tablex.flatten(responsePercentilesTbl);
display(responsePercentilesTbl);

% Define the output path for saving the results
outputPath = fullfile(outputFolder, "responsePercentiles");

% Save the shock response results in percentiles as MAT and/or CSV and/or XLSX
% files
% save(outputPath + ".mat", "responsePercentilesTbl");
tablex.writetimetable(responsePercentilesTbl, outputPath + ".csv");
tablex.writetimetable(responsePercentilesTbl, outputPath + ".xlsx");

% Plot the shock response results in percentiles
figureHandles = chartpack.responsePercentiles( ...
    responsePercentilesTbl, structModel ...
    , "figureTitle", "Shock response percentiles" ...
    , "figureLegend", percentilesLegend ...
);

% Save the figures as a PDF
chartpack.printFiguresPDF(figureHandles, outputPath);


