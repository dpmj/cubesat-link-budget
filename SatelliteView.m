
close all;
clearvars;

%% Today's date and stop time.

startTime = datetime(2022, 4, 1, 4, 0, 0);  % 2022-04-01 04:00:00
stopTime = startTime + hours(72);  % after 72 hours
sampleTime = 60;  % seconds


%% Create scenario

sc = satelliteScenario(startTime, stopTime, sampleTime);


%% Satellite

sat = satellite(sc, "EstigiaTLE.tle");

show(sat)
groundTrack(sat,"LeadTime",1200);


%% Ground Stations

name = "Valencia UPV";
lat = 39.47943; 
lon = -0.34230;

gs = groundStation(sc, Name=name, Latitude=lat, Longitude=lon);



%% Visualization

play(sc)
