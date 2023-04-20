
close all;
clearvars;

%% Constants.

R_E = 6371;  % [km] Earth's average radius


%% Today's date and stop time.

startTime = datetime(2022, 4, 1, 4, 0, 0);  % [date, time] 2022-04-01 04:00:00
stopTime = startTime + hours(24);  % Simulation time, stop after 24 hours
sampleTime = 60;  % [s] Simulation step


%% Create scenario

scenario = satelliteScenario(startTime, stopTime, sampleTime);  % Satellite scenario


%% Satellite

sat = satellite(scenario, "EstigiaTLE.tle");  % Import satellite TLE file

groundTrack(sat, LeadTime=1200);  % Show satellite ground tracks with a lead of 20 min

% sat_comms_fov = 179;  % [degrees] Satellite's field of view of the TT&C module 
% If 179º: Access will be constrained by the GS minimum elevation angle.

% sat_comms = conicalSensor(sat, name="TT&C", MaxViewAngle=sat_comms_fov);  % Setup sensor
% sat_fov_track = fieldOfView(sat_comms);



%% Ground Station

name = "Valencia";
lat = 39.47943;  % North latitude
lon = -0.34230;  % East longitude
altitude = 50;  % [m] Altitude of the GS
elevation_min_angle = 12;  % [degrees] Minimum angle of elevation in GS

% Setup Ground Station
gs = groundStation(scenario, Name=name, Latitude=lat, Longitude=lon, Altitude=altitude, ...
                   MinElevationAngle=elevation_min_angle);

% Create a conical sensor (area over the globe) using the elevation. 
% This is equivalent to the satellite coverage with the minimum elevation. 
% Inexact, use only for graphical representation.
% Assuming LEO, low orbit eccentricity

sat_pos_start = states(sat, startTime, CoordinateFrame="geographic");  % [lat,long,height]
sat_height = sat_pos_start(3) / 1e3;  % [km] Height of the satellite over Earth's surface

% [degrees] Field-of-view angle 
fov_angle = 2 * (180/pi) * asin(R_E/(R_E + sat_height) * cos(elevation_min_angle * (pi/180)));

% FOV angle calculus from:
% "The Coverage Analysis for Low Earth Orbiting Satellites at Low Elevation"
% Shkelzen Cakaj, Bexhet Kamo, Algenti Lala, Alban Rakipi
% (IJACSA) International Journal of Advanced Computer Science and Applications, 214

% Conical sensor
sat_coverage = conicalSensor(sat, MaxViewAngle=fov_angle);
sat_fov = fieldOfView(sat_coverage);



%% Compute access to the satellite from ground station

ac_sat_gs = access(sat, gs);  % Compute access to the satellite
ac_sat_intervals = accessIntervals(ac_sat_gs);  % Access intervals from GS

[ac_sat_plot_data, ac_sat_time] = accessStatus(ac_sat_gs);  % Plot the access intervals

ac_sat_plot_data = double(ac_sat_plot_data);  % Casting for the operation below
ac_sat_plot_data(ac_sat_plot_data == false) = NaN;  % Filter non-visibility intervals

% Plot points in time where there is access to the satellite.
figure();
plot(ac_sat_time, ac_sat_plot_data, ' .');
title("coverage intervals");


%% Simulation visualization

v = satelliteScenarioViewer(scenario, ShowDetails=true);
sat.ShowLabel = true;
gs.ShowLabel = true;

show(sat);  % Show satellite

play(scenario);



