%% SATELLITE LINK BUDGET ANALYSIS - sweep
% Juan Del Pino Mena
% Version 1, December 2023
% 
% Script that sweeps a parameter using link_budget_sim_func()
% 
% CHANGELOG
%
% Version 1: First iteration


close all;
clearvars;



%% CONFIGURATION

% simulation -----------------------------------------------------------------------------

config.sim.startTime = datetime(2023, 4, 1, 7, 30, 0);  % [datetime object] Simulation start time
config.sim.stopTime = config.sim.startTime + hours(4);  % [datetime object] Simulation stop time
config.sim.sampleTime = 10;  % [s] Simulation step, sample time in seconds

% general communication parameters -------------------------------------------------------

config.comms.freq_Hz = 868e6;   % [Hz] Transmission center frequency
config.comms.bitrate_bps = 293;  % [bps] Transmission bitrate
config.comms.symbolrate_Sps = 1/0.032768;  % [Sps] Symbol rate
config.comms.bandwidth_Hz = 125e3;  % [Hz] Transmission bandwidth

% ground station -------------------------------------------------------------------------

config.gs.name = "UPV GS";  % [string] Name of the ground station
config.gs.lat_N_deg = 39.47943;  % [deg] latitude
config.gs.lon_E_deg = -0.34230;  % [deg] longitude
config.gs.altitude_m = 50;  % [m] altitude
config.gs.elevation_min_angle_deg = 12;  % [deg] minimum elevation angle

config.gs.ant.gain_dBi = 0;  % [dBi] Antenna Gain in the ground station
config.gs.ant.ambient_temp_K = 300;  % [dBi] Antenna Gain in the ground station
config.gs.ant.noise_temp_K = 290;  % [dBi] Antenna Gain in the ground station

config.gs.rx.loss_dB = 1;  % [dB] RX system power loss in the ground station
config.gs.rx.nf_dB = 6;  % [dB] RX system noise figure in the ground station

config.gs.tx.power_dBm = 22;  % [dBm] Transmission power in the Ground Station tx
config.gs.tx.loss_dB = 1;  % [dB] TX system power loss in the ground station

% satellite ------------------------------------------------------------------------------

config.sat.tle_file = "EstigiaTLE.tle";  % [string] path to the satellite's TLE file

config.sat.ant.gain_dBi = 0;  % [dBi] Antenna Gain in the ground station
config.sat.ant.ambient_temp_K = 350;  % [dBi] Antenna Gain in the ground station
config.sat.ant.noise_temp_K = 2340;  % [dBi] Antenna Gain in the ground station

config.sat.rx.loss_dB = 1;  % [dB] RX system power loss in the ground station
config.sat.rx.nf_dB = 3;  % [dB] RX system noise figure in the ground station

config.sat.tx.power_dBm = 33;  % [dBm] Transmission power in the Ground Station tx
config.sat.tx.loss_dB = 1;  % [dB] TX system power loss in the ground station

% p618 model -----------------------------------------------------------------------------

config.p618.GasAnnualExceedance = 0.1;  % Avg annual time percentage of excess for gas att
config.p618.CloudAnnualExceedance = 0.1;  % Avg annual time percentage of excess for cloud att
config.p618.RainAnnualExceedance = 0.001;  % Avg annual time percentage of excess for rain att
config.p618.ScintillationAnnualExceedance = 0.01;  % Avg annual time percentage of excess scintillation
config.p618.TotalAnnualExceedance = 0.001;  % Avg annual time % of excess for total att
config.p618.PolarizationTiltAngle = 0;  % [degrees, -90 to 90]  Polarization tilt angle
config.p618.AntennaDiameter = 1;  % [m] Physical diameter of the antenna. Default = 1 m
config.p618.AntennaEfficiency = 0.5;  % [0-1] Antenna efficiency
config.p618.PolLoss_dB = 3;  % [dB] Polarization loss. Worst-case: 3 dB



% ----------------------------------------------------------------------------------------
%% EXECUTION

warning('off','all');  % allows faster code execution - warnings are known
addpath("functions/");  % add function to path

% tic;
output = link_budget_sim_func(config);
% toc;  % approx 3.8 seconds per run


fprintf("\n\n-----------------------------------------\nATMOSPHERIC LOSSES:\n");
fprintf("\nGaseous attenuation = %.3f dB", output.p618.Ag_dB);
fprintf("\nCloud and fog attenuation = %.3f dB", output.p618.Ac_dB);
fprintf("\nRain attenuation = %.3f dB", output.p618.Ar_dB);
fprintf("\nAttenuation due to tropospheric scintillation = %.3f dB", output.p618.As_dB);
fprintf("\nTotal atmospheric attenuation = %.3f dB", output.p618.At_dB);
fprintf("\nCross-polarization discrimination = %.3f dB", output.p618.xpol_discr_dB);
fprintf("\nSky noise temperature = %.3f K", output.p618.temp_sky_K);
fprintf("\nPolarization loss (MANUAL PICK) = %.3f dB", output.p618.pol_loss_dB);
fprintf("\nTOTAL LOSSES = Total att. + Pol. loss = %.3f dB\n\n", output.p618.loss_total_dB);





figure();
subplot(2, 1, 1);
plot(output.latency.time, output.latency.delay_s(1, :) .* 1e3);  % in milliseconds
xlim([output.latency.time(1), output.latency.time(end) .* 1e3]);
title("Satellite Latency in observation time", interpreter="latex");
xlabel("Simulation Time", interpreter="latex");
ylabel("Latency (ms)", interpreter="latex");
grid on; grid minor;

subplot(2, 1, 2);
plot(doppler_time, doppler_fshift * 1e-3);  % kHz
xlim([doppler_time(1), doppler_time(end)]);
title("Doppler Shift in observation time", interpreter="latex");
xlabel("Simulation Time", interpreter="latex");
ylabel("Doppler Frequency Shift (kHz)", interpreter="latex");
grid on; grid minor;




figure;

subplot(2, 1, 1);
plot(link_UL_time, link_UL_rxpower_dBm); 
grid on; grid minor;
xlabel("Simulation time", interpreter="latex");
ylabel("Received signal strength (dBm)", interpreter="latex");
title("Uplink", interpreter="latex");
yline(rx_limit_sensitivity_sat, Color="#5e5e5e", LineStyle="--", Label="Sat sensitivity");

subplot(2, 1, 2);
plot(link_DL_time, link_DL_rxpower_dBm); 
grid on; grid minor;
xlabel("Simulation time", interpreter="latex");
ylabel("Received signal strength (dBm)", interpreter="latex");
title("Downlink", interpreter="latex");
yline(rx_limit_sensitivity_gs, Color="#5e5e5e", LineStyle="--", Label="GS sensitivity");

sgtitle("Received signal strength in uplink and downlink", interpreter="latex");



figure;

subplot(2, 1, 1);
plot(cnr_time, ul_CNR); 
grid on; grid minor;
xlabel("Simulation time", interpreter="latex");
ylabel("Carrier-to-Noise Ratio (dB)", interpreter="latex");
title("Uplink", interpreter="latex");
yline(rx_limit_cnr_sat, Color="#5e5e5e", LineStyle="--", Label="Sat CNR limit");

subplot(2, 1, 2);
plot(cnr_time, dl_CNR); 
grid on; grid minor;
xlabel("Simulation time", interpreter="latex");
ylabel("Carrier-to-Noise Ratio (dB)", interpreter="latex");
title("Downlink", interpreter="latex");
yline(rx_limit_cnr_gs, Color="#5e5e5e", LineStyle="--", Label="GS CNR limit");

sgtitle("Carrier-to-Noise Ratio (CNR) in uplink and downlink", interpreter="latex");







