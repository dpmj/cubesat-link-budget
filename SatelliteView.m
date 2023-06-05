%% SATELLITE LINK BUDGET ANALYSIS
% Juan Del Pino Mena - 2023
% PLUTON UPV


close all;
clear all;  % To be sure



% ----------------------------------------------------------------------------------------
%% Constants.

R_E = 6371;  % [km] Earth's average radius


% ----------------------------------------------------------------------------------------
%% Configuration

% Transmission and reception parameters

freq_Hz = 433e6;  % [Hz] Transmission frequency 
bitrate = 300;  % [b/s] Transmission bitrate
bandwidth = 250e3;  % [Hz] Transmission bandwidth

% rx_SNR = -20;  % Required SNR


% Ground station
tx_power_gs_dBm = 30;  % [dBm] Transmission power in the Ground Station tx
tx_ant_gain_gs_dB = 0;  % [dBi] Antenna Gain in the GS tx
tx_loss_gs_dB = 2;  % [dB] Power loss in the TX in the GS tx

rx_GNTR_gs_dB = 5;  % [dB/K] Gain-to-Noise-Temperature Ratio in the GS rx
rx_EbNo_gs_dB = 0;  % [dB] Required Eb/No (Bit Energy over Noise Ratio) in GS rx 

% Satellite
tx_power_sat_dBm = 30;  % [dBm] Transmission power in the satellite tx
tx_ant_gain_sat_dB = 0;  % [dBi] Antenna Gain in the sat tx
tx_loss_sat_dB = 2;  % [dB] Power loss in the TX in the sat tx

rx_GNTR_sat_dB = 3;  % [dB] Gain-to-Noise-Temperature Ratio in the sat rx
rx_EbNo_sat_dB = 0;  % [dB] Required Eb/No (Bit Energy over Noise Ratio) in sat rx 


% Simulation time, stop time and step

startTime = datetime(2023, 4, 1, 7, 30, 0);  % [date, time] 2023-04-01 04:00:00
stopTime = startTime + seconds(1200);  % Simulation time, stop after X seconds
sampleTime = 5;  % [s] Simulation step



% ----------------------------------------------------------------------------------------
%% Create scenario

scenario = satelliteScenario(startTime, stopTime, sampleTime);  % Satellite scenario

% Satellite

sat = satellite(scenario, "EstigiaTLE.tle");  % Import satellite TLE file

groundTrack(sat, LeadTime=1200);  % Show satellite ground tracks with a lead of 20 min

% sat_comms_fov = 179;  % [degrees] Satellite's field of view of the TT&C module 
% If 179º: Access will be constrained by the GS minimum elevation angle.
% sat_comms = conicalSensor(sat, name="TT&C", MaxViewAngle=sat_comms_fov);  % Setup sensor
% sat_fov_track = fieldOfView(sat_comms);

% Ground Station

name = "Valencia";
lat = 39.47943;  % North latitude
lon = -0.34230;  % East longitude
altitude = 50;  % [m] Altitude of the GS
elevation_min_angle = 12;  % [degrees] Minimum angle of elevation in GS

% Setup Ground Station
gs = groundStation(scenario, Name=name, Latitude=lat, Longitude=lon, ...
                   Altitude=altitude, MinElevationAngle=elevation_min_angle);

% Create a conical sensor (area over the globe) using the elevation. 
% This is equivalent to the satellite coverage with the minimum elevation. 
% Inexact, use only for graphical representation.
% Assuming LEO, low orbit eccentricity

sat_pos_start = states(sat, startTime, CoordinateFrame="geographic");  % [lat,long,height]
sat_height = sat_pos_start(3) / 1e3;  % [km] Height of the satellite over Earth's surface

% [degrees] Field-of-view angle 
fov_angle = 2 * (180/pi) * asin(R_E/(R_E+sat_height) * cos(elevation_min_angle*(pi/180)));

% FOV angle calculus from:
% "The Coverage Analysis for Low Earth Orbiting Satellites at Low Elevation"
% Shkelzen Cakaj, Bexhet Kamo, Algenti Lala, Alban Rakipi
% International Journal of Advanced Computer Science and Applications (IJACSA)

% Add conical sensor to represent the satellite footprint
sat_coverage = conicalSensor(sat, MaxViewAngle=fov_angle);
sat_fov = fieldOfView(sat_coverage);



% ----------------------------------------------------------------------------------------
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
xlabel("Simulation Time");
ylabel("Access status");
grid on; grid minor;



% ----------------------------------------------------------------------------------------
%% Add TX and RX to the satellite and the Ground Station

% ----------------------------------------------------------------------------------------
% Ground Station

gs_antenna = arrayConfig(Size=[1, 1]);  % Isotropic antenna element from phasedArray T.Box

% Transmitter
gs_transmitter = transmitter(gs, ...  % parent element
                             Antenna=gs_antenna, ...  % Antenna element
                             Frequency=freq_Hz, ...  % Center frequency (Hz)
                             Power=tx_power_gs_dBm - 30, ...  % TX power (dBW)
                             BitRate=bitrate * 1e-6, ...  % Bitrate (Mbps)
                             SystemLoss=tx_loss_gs_dB, ...  % Losses in the equipment (dB)
                             MountingAngles=[0; 0; 0], ...  % Orientation (degrees)
                             Name="GS TX");  % Name

% Receiver
gs_receiver = receiver(gs, ...  % parent element
                       Antenna=gs_antenna, ...  % Antenna element
                       MountingAngles=[0; 0; 0], ...  % Orientation (degrees)
                       Name="GS RX");  % Name

pattern(gs_receiver, freq_Hz, Size=1e5);  % Show radiation pattern in 3D view


% ----------------------------------------------------------------------------------------
% Satellite

sat_antenna = arrayConfig(Size=[1, 1]);  % Isotropic antenna 

% Transmitter
sat_transmitter = transmitter(sat, ...
                              Antenna=sat_antenna, ...  % Antenna element
                              Frequency=freq_Hz, ...  % Center frequency (Hz)
                              Power=tx_power_sat_dBm - 30, ...  % TX power (dBW)
                              BitRate=bitrate * 1e-6, ...  % Bitrate (Mbps)
                              SystemLoss=tx_loss_gs_dB, ...  % Losses in the equipment, dB
                              MountingAngles=[0; 0; 0], ...  % Orientation (degrees)
                              Name="SAT TX");  % Name

% Receiver
sat_receiver = receiver(sat, ...  % parent element
                        Antenna=sat_antenna, ...  % Antenna element
                        MountingAngles=[0; 0; 0], ...  % Orientation (degrees)
                        Name="SAT RX");  % Name

pattern(sat_transmitter, Size=1e5);  % Show radiation pattern in 3D view



% ----------------------------------------------------------------------------------------
%% Link budget analysis - EbNo

% Uplink
link_UL = link(gs_transmitter, sat_receiver);  % Perform uplink analysis
[link_UL_EbNo_dB, link_UL_time] = ebno(link_UL);
link_UL_SNR_dB = link_UL_EbNo_dB + 10 * log10(bitrate / bandwidth);

% Downlink
link_DL = link(sat_transmitter, gs_receiver);  % Perform downlink analysis
[link_DL_EbNo_dB, link_DL_time] = ebno(link_DL);
link_DL_SNR_dB = link_DL_EbNo_dB + 10 * log10(bitrate / bandwidth);

% Plot
figure;

subplot(1, 2, 1);
plot(link_UL_time, link_UL_EbNo_dB); 
hold on;
plot(link_UL_time, link_UL_SNR_dB);
grid on; grid minor;
xlabel("Simulation time");
ylabel("Received Eb/No or SNR (dB)");
title("Uplink Eb/No");

subplot(1, 2, 2);
plot(link_DL_time, link_DL_EbNo_dB); 
hold on;
plot(link_DL_time, link_DL_SNR_dB);
grid on; grid minor;
xlabel("Simulation time");
ylabel("Received Eb/No or SNR (dB)");
title("Downlink Eb/No");



% ----------------------------------------------------------------------------------------
%% Simulation visualization

v = satelliteScenarioViewer(scenario, ShowDetails=true);
sat.ShowLabel = true;
gs.ShowLabel = true;

show(sat);  % Show satellite

play(scenario);



% ----------------------------------------------------------------------------------------
%% Latency
% This function is only available on MatLab R2023a onwards
% https://www.mathworks.com/help/satcom/ref/matlabshared.satellitescenario.satellite.latency.html 

% [lat_delay, lat_time] = latency(sat, gs);  % Compute latency (builtin func in R2023a)
% 
% % Plot latency 
% figure();
% plot(lat_time, lat_delay(1, :) .* 1e3);  % in milliseconds
% 
% xlim([lat_time(1), lat_time(end)]);
% title("First Satellite's Latency vs. Time");
% xlabel("Simulation Time");
% ylabel("Latency (ms)");
% grid on; grid minor;



% ----------------------------------------------------------------------------------------
%% P618 Propagation model

% cfg = p618Config;
% cfg.Frequency = 25e9;
%  % Signal frequency in Hz
% cfg.ElevationAngle = 45;
% cfg.Latitude = 30;
%  % North direction
% cfg.Longitude = 120;
%  % East direction
% cfg.TotalAnnualExceedance = 0.001;
%  % Time percentage of excess for the total
% % Attenuation per annum
% cfg.AntennaEfficiency = 0.65;
% disp(cfg);


