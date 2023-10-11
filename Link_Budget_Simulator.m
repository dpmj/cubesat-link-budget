%% SATELLITE LINK BUDGET ANALYSIS
% Juan Del Pino Mena - 2023
% PLUTON UPV
% NOTE: This program requires MatLab >= R2023a

% Link budget estimator for LEO satellites. This script computes access intervals in a
% given simulation time, the latency and Doppler frequency shift. Also estimates the
% propagation and atmospheric losses, and provides a complete link budget.


close all;
clearvars;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants

k_boltzmann = 1.380649e-23;  % [J/K] Boltzmann's constant
T_0 = 290; % [K] room temperature for noise calculus, usually 290 K

R_E = 6371;  % [km] Earth's average radius
% Not for orbit propagation, but for auxiliary calculations.


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Configuration parameters

% ----------------------------------------------------------------------------------------
% Transmission and reception parameters

freq_Hz = 433e6;  % [Hz] Transmission frequency 
bitrate_bps = 300;  % [b/s] Transmission bitrate
symbolrate_Sps = 1/0.032768;  % [Sps] Symbol rate (LoRa: Ts = 2^SF / BW)
bandwidth_Hz = 125e3;  % [Hz] Transmission bandwidth


% ----------------------------------------------------------------------------------------
% Ground Station

gs_name = "UPV GS";
gs_lat_N = 39.47943;  % North latitude
gs_lon_E = -0.34230;  % East longitude
gs_altitude = 50;  % [m] Altitude of the GS
elevation_min_angle = 12;  % [degrees] Minimum angle of elevation in GS

% Power, gain and losses

tx_power_gs_dBm = 22;  % [dBm] Transmission power in the Ground Station tx
tx_ant_gain_gs_dB = 0;  % [dBi] Antenna Gain in the GS tx
tx_loss_gs_dB = 3;  % [dB] TX system power loss in the ground station
rx_loss_gs_dB = 3;  % [dB] RX system power loss in the ground station

% Noise 

rx_nf_gs_dB = 3;  % [dB] System noise figure in RX in the ground station
% rx_noise_temp_gs_K = (10^(rx_nf_gs_dB/10) - 1) * T_0;  % [K] Equiv. noise temp of GS rx

ant_ambient_temp_gs_K = 290;  % [K] Antenna ambient temperature in the ground station
ant_noise_temp_gs_K = 290;  % [K] Antenna noise temperature in the ground station
% The equivalent noise temperature seen at the receive output of the antenna. Describes 
% how much noise an antenna produces in a given environment. For uplink 290 K is 
% considered as the worst case, and for downlink is 2340 K.

% sys_temp_gs_K = rx_noise_temp_gs_K + ant_noise_temp_gs_K;  % [K] system equiv temp GS rx

% rx_PN_gs_dB = 10 * log10(k_boltzmann * ant_noise_temp_K * bandwidth);
% [dB] noise power at GS antenna

% G/T (Gain-to-Noise-Temperature Ratio, GNTR) relates the receive antenna gain and the 
% system noise temperature. system_noise_temp = antenna_noise_temp + receiver_noise_temp
% If we are not measuring with an receiver then system_noise_temp = antenna_noise_temp
% This is not a representative value for calculating G/T, since it relates to the receive 
% performance of both antenna and receiver. 

% G/T = GR − Nf − 10 * log10( T0 + (Ta − T0) * 10^(−0.1 * Nf) )
% GR is the receive antenna gain in dBi
% Nf is the noise figure in dB
% T0 is the ambient temperature in degrees Kelvin
% Ta is the antenna temperature in degrees Kelvin
% https://www.mathworks.com/help/satcom/ug/nb-iot-ntn-link-budget-analysis.html

rx_GNTR_gs_dB_K = tx_ant_gain_gs_dB - rx_nf_gs_dB - 10 * log10(ant_ambient_temp_gs_K + ...
    (ant_noise_temp_gs_K - ant_ambient_temp_gs_K) * 10^(-0.1 * rx_nf_gs_dB));  % [dB/K]

% rx_EbNo_gs_dB = 0;  % [dB] Required Eb/No (Bit Energy over Noise Ratio) in GS rx 


% ----------------------------------------------------------------------------------------
% Satellite

% Power, gain and losses

tx_power_sat_dBm = 22;  % [dBm] Transmission power in the satellite tx
tx_ant_gain_sat_dB = 0;  % [dBi] Antenna Gain in the sat tx
tx_loss_sat_dB = 3;  % [dB] System power loss in TX in the satellite

rx_loss_sat_dB = 3;  % [dB] System power loss in RX in the satellite
rx_GNTR_sat_dB = 3;  % [dB] Gain-to-Noise-Temperature Ratio in the sat rx

% Noise 

rx_nf_sat_dB = 3;  % [dB] System noise figure in RX in the ground station
% rx_noise_temp_sat_K = (10^(rx_nf_sat_dB/10) - 1) * T_0;  % [K] Equiv. noise temp of sat rx

ant_ambient_temp_sat_K = 100;  % [K] Antenna ambient temperature in the satellite
ant_noise_temp_sat_K = 2340;  % [K] Antenna noise temperature in the satellite
% sys_temp_sat_K = rx_noise_temp_sat_K + ant_noise_temp_sat_K;  % [K] syst equiv temp sat rx

% [dB/K] G/T in the sat rx
rx_GNTR_sat_dB_K = tx_ant_gain_sat_dB - rx_nf_sat_dB ... 
    - 10 * log10(ant_ambient_temp_sat_K ...
    + (ant_noise_temp_sat_K - ant_ambient_temp_sat_K) * 10^(-0.1 * rx_nf_sat_dB)); 

% rx_EbNo_sat_dB = 0;  % [dB] Required Eb/No (Bit Energy over Noise Ratio) in sat rx 


% ----------------------------------------------------------------------------------------
% Simulation time, stop time and step

startTime = datetime(2023, 4, 1, 7, 30, 0);  % [date, time] 2023-04-01 07:30:00
stopTime = startTime + hours(4);  % Simulation time, stop after 4 hours
sampleTime = 10;  % [s] Simulation step



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create scenario

scenario = satelliteScenario(startTime, stopTime, sampleTime);  % Satellite scenario

% ----------------------------------------------------------------------------------------
% Satellite

sat = satellite(scenario, "EstigiaTLE.tle");  % Import satellite TLE file
groundTrack(sat, LeadTime=1200);  % Show satellite ground tracks with a lead of 20 min

% ----------------------------------------------------------------------------------------
% Ground Station
gs = groundStation(scenario, Name=gs_name, Latitude=gs_lat_N, Longitude=gs_lon_E, ...
                   Altitude=gs_altitude, MinElevationAngle=elevation_min_angle);

% Create a conical sensor (area over the globe) using the elevation. This is equivalent to 
% the satellite coverage with the minimum elevation. May be inexact, use only for 
% graphical representation. Assuming LEO, low orbit eccentricity.

% FOV angle calculus from: "The Coverage Analysis for Low Earth Orbiting Satellites at Low
% Elevation" S. Cakaj, B. Kamo, A. Lala, A. Rakipi. International Journal of Advanced 
% Computer Science and Applications (IJACSA)

sat_pos_start = states(sat, startTime, CoordinateFrame="geographic");  % [lat,long,height]
sat_height = sat_pos_start(3) / 1e3;  % [km] Height of the satellite over Earth's surface

% [degrees] Field-of-view angle 
fov_angle = 2 * (180/pi) * asin(R_E/(R_E+sat_height) * cos(elevation_min_angle*(pi/180)));

% Add conical sensor to represent the satellite footprint
sat_coverage = conicalSensor(sat, MaxViewAngle=fov_angle);
sat_fov = fieldOfView(sat_coverage);




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute access to the satellite from ground station

ac_sat_gs = access(sat, gs);  % Compute access to the satellite
ac_sat_intervals = accessIntervals(ac_sat_gs);  % Access intervals from GS

[ac_sat_plot_data, ac_sat_time] = accessStatus(ac_sat_gs);  % Plot the access intervals

ac_sat_plot_data = double(ac_sat_plot_data);  % Casting for the operation below
% ac_sat_plot_data(ac_sat_plot_data == 0) = NaN;  % Filter non-visibility intervals

% ----------------------------------------------------------------------------------------
% Plot points in time where there is access to the satellite.

figure();
plot(ac_sat_time, ac_sat_plot_data, '.-');
title("Observation intervals: line-of-sight with $\varepsilon >= \varepsilon_{min}$" + ...
      ". [1: access; 0: no access]", interpreter="latex");
xlabel("Simulation Time", interpreter="latex");
ylabel("Access status", interpreter="latex");
ylim([-0.5, 1.5]);
grid on;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add TX and RX to the satellite and the Ground Station

% ----------------------------------------------------------------------------------------
% Ground Station

gs_antenna = arrayConfig(Size=[1, 1]);  % Isotropic antenna element from phasedArray T.Box

% Transmitter
gs_transmitter = transmitter( ...
    gs, ...  % parent element
    Antenna = gs_antenna, ...  % Antenna element
    Frequency = freq_Hz, ...  % Center frequency (Hz)
    Power = tx_power_gs_dBm - 30, ...  % TX power (dBW)
    BitRate = bitrate_bps * 1e-6, ...  % Bitrate (Mbps)
    SystemLoss = tx_loss_gs_dB, ...  % Equipment losses (dB)
    MountingAngles = [0; 0; 0], ...  % Orientation (degrees)
    Name = "GS TX" ...  % Name
);

% Receiver
gs_receiver = receiver( ...
    gs, ...  % parent element
    Antenna = gs_antenna, ...  % Antenna element
    MountingAngles = [0; 0; 0], ...  % Orientation (degrees)
    SystemLoss = rx_loss_gs_dB, ...  % Equipment losses (dB)
    Name = "GS RX" ...  % Name
);

pattern(gs_receiver, freq_Hz, Size=1e5);  % Show radiation pattern in 3D viewer


% ----------------------------------------------------------------------------------------
% Satellite

sat_antenna = arrayConfig(Size=[1, 1]);  % Isotropic antenna 

% Transmitter
sat_transmitter = transmitter( ...
    sat, ...  % parent element
    Antenna = sat_antenna, ...  % Antenna element
    Frequency = freq_Hz, ...  % Center frequency (Hz)
    Power = tx_power_sat_dBm - 30, ...  % TX power (dBW)
    BitRate = bitrate_bps * 1e-6, ...  % Bitrate (Mbps)
    SystemLoss = tx_loss_sat_dB, ...  % Equipment losses (dB)
    MountingAngles = [0; 0; 0], ...  % Orientation (degrees)
    Name = "SAT TX" ...  % Name
);

% Receiver
sat_receiver = receiver( ...
    sat, ...  % parent element
    Antenna = sat_antenna, ...  % Antenna element
    MountingAngles = [0; 0; 0], ...  % Orientation (degrees)
    SystemLoss = rx_loss_sat_dB, ...  % Equipment losses (dB)
    Name = "SAT RX" ...  % Name
);

pattern(sat_transmitter, Size=1e5);  % Show radiation pattern in 3D viewer



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Latency and Doppler Shift
% These functions are only available on MatLab R2023a onwards
% Dynamic - computes latency and doppler shift over the complete simulation time

[latency_delay, latency_time] = latency(sat, gs);
[doppler_fshift, doppler_time, doppler_info] = dopplershift(sat, gs, Frequency=freq_Hz);


figure();
subplot(1, 2, 1);
plot(latency_time, latency_delay(1, :) .* 1e3);  % in milliseconds
xlim([latency_time(1), latency_time(end)]);
title("Satellite Latency in observation time", interpreter="latex");
xlabel("Simulation Time", interpreter="latex");
ylabel("Latency (ms)", interpreter="latex");
grid on; grid minor;

subplot(1, 2, 2);
plot(doppler_time, doppler_fshift * 1e-3);  % kHz
xlim([doppler_time(1), doppler_time(end)]);
title("Doppler Shift in observation time", interpreter="latex");
xlabel("Simulation Time", interpreter="latex");
ylabel("Doppler Frequency Shift (kHz)", interpreter="latex");
grid on;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% P618 propagation losses - atmospheric losses
% Estimate the atmospheric losses expected near Valencia.
% Static - does not recalculate for every simulation point

% ----------------------------------------------------------------------------------------
% Download and extract the digital maps, if not available on path

maps = exist("maps.mat","file");
p836 = exist("p836.mat","file");
p837 = exist("p837.mat","file");
p840 = exist("p840.mat","file");
matFiles = [maps p836 p837 p840];

if ~all(matFiles)
    if ~exist("ITURDigitalMaps.tar.gz","file")
        url = "https://www.mathworks.com/supportfiles/spc/P618/ITURDigitalMaps.tar.gz";
        websave("ITURDigitalMaps.tar.gz",url);
        untar("ITURDigitalMaps.tar.gz");
    else
        untar("ITURDigitalMaps.tar.gz");
    end
    addpath(cd)
end


% ----------------------------------------------------------------------------------------
% Configuration object

% KNOWN PROBLEMS:
% - The lowest frequency the model admits is 1 GHz.
% - Where is the polarization loss?
% - The cross-polarization discrimination prediction method is valid for the frequency 
%   values in the range 4 to 55 GHz. For less than 4 GHz, the 4 GHz point will be used.

p618cfg = p618Config(...
    Frequency = 1e9, ...  % [Hz] Carrier frequency
    ElevationAngle = elevation_min_angle, ...  % [degrees] Elevation Angle, worst-case
    Latitude = gs_lat_N, ...  % Degrees North
    Longitude = gs_lon_E, ...  % Degrees East
    GasAnnualExceedance = 0.1, ...  % Average annual time percentage of excess for gas att
    CloudAnnualExceedance = 0.1, ...  % Avg annual time percentage of excess for cloud att
    RainAnnualExceedance = 0.001, ...  % Avg annual time percentage of excess for rain att
    ScintillationAnnualExceedance = 0.01, ...  % Avg annual time % of excess scintillation
    TotalAnnualExceedance = 0.001, ...  % Avg annual time % of excess for total att
    PolarizationTiltAngle = 0, ...  % [degrees, -90 to 90]  Polarization tilt angle
    AntennaDiameter = 1, ...  % [m] Physical diameter of the antenna. Default = 1 m
    AntennaEfficiency = 0.5 ...  % [0-1] Antenna efficiency
);

% ----------------------------------------------------------------------------------------
% Calculate Earth-space propagation losses, cross-polarization discrimination and sky 
% noise temperature with the above configuration
[p618_atm_loss, p618_xpol_discr, p618_temp_sky] = p618PropagationLosses(p618cfg);

% p618_atm_loss.Ag - Gaseous attenuation (dB)
% p618_atm_loss.Ac - Cloud and fog attenuation (dB)
% p618_atm_loss.Ar - Rain attenuation (dB)
% p618_atm_loss.As - Attenuation due to tropospheric scintillation (dB)
% p618_atm_loss.At - Total atmospheric attenuation (dB)
% p618_xpol_discr  - Cross-polarization discrimination (dB) not exceeded for the 
%                    percentage of the RainAnnualExceedance.
% p618_temp_sky    - Sky noise temperature (K) at the ground station antenna.

% ----------------------------------------------------------------------------------------
% Polarization loss

p618_pol_loss = 3;  % [dB] Polarization loss. Worst-case: 3 dB

% ----------------------------------------------------------------------------------------
% Total losses
p618_loss_total = p618_atm_loss.At + p618_pol_loss;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CNR: Carrier-to-noise ratio for configured satellite link budget parameters
% Static by default, but here is forced to be dynamic: loop computes CNR in every access=1
% simulation point

% ----------------------------------------------------------------------------------------
% UPLINK: Ground Station -> Satellite 

ul_CNR_cfg = satelliteCNRConfig(...
    TransmitterPower = tx_power_gs_dBm - 30, ...  % [dBW] Transmit power 
    TransmitterSystemLoss = tx_loss_gs_dB, ...  % [dB] Losses in the transmitter
    TransmitterAntennaGain = tx_ant_gain_gs_dB, ...  % [dBi] Antenna gain
    Distance = 1700, ...  % [km] Distance. Worst-case
    Frequency = freq_Hz * 1e-9, ...  % [GHz] Center frequency
    MiscellaneousLoss = p618_loss_total, ...  % [dB] Misc attenuation
    GainToNoiseTemperatureRatio = rx_GNTR_sat_dB_K, ...  % [dB/K] G/T ratio
    ReceiverSystemLoss = rx_loss_sat_dB, ...  % [dB] Receiver system loss
    BitRate = bitrate_bps * 1e-6, ...  % [Mbps] Bit rate
    SymbolRate = symbolrate_Sps * 1e-6, ...  % [MSps] Symbol rate
    Bandwidth = bandwidth_Hz * 1e-6 ...  [MHz] Bandwidth
);

[ul_CNR, ul_CNR_info] = satelliteCNR(ul_CNR_cfg);


% ----------------------------------------------------------------------------------------
% DOWNLINK: Satellite -> Ground Station

dl_CNR_cfg = satelliteCNRConfig(...
    TransmitterPower = tx_power_sat_dBm - 30, ...  % [dBW] Transmit power 
    TransmitterSystemLoss = tx_loss_sat_dB, ...  % [dB] Losses in the transmitter
    TransmitterAntennaGain = tx_ant_gain_sat_dB, ...  % [dBi] Antenna gain
    Distance = 1200, ...  % [km] Distance. Worst-case
    Frequency = freq_Hz * 1e-9, ...  % [GHz] Center frequency
    MiscellaneousLoss = p618_loss_total, ...  % [dB] Misc attenuation
    GainToNoiseTemperatureRatio = rx_GNTR_gs_dB_K, ...  % [dB/K] G/T ratio
    ReceiverSystemLoss = rx_loss_gs_dB, ...  % [dB] Receiver system loss
    BitRate = bitrate_bps * 1e-6, ...  % [Mbps] Bit rate
    SymbolRate = symbolrate_Sps * 1e-6, ...  % [MSps] Symbol rate
    Bandwidth = bandwidth_Hz * 1e-6 ...  [MHz] Bandwidth
);

[dl_CNR, dl_CNR_info] = satelliteCNR(dl_CNR_cfg);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Link budget analysis: Received Power

% EIRP: Equivalent Isotropic Radiated Power: EIRP = P_tx + G_ant_tx
% PISO: Received isotropic power in the antenna: PISO_rx = EIRP_tx - FSPL 
% PRI:  Received power at receiver input: PRI_rx = EIRP_tx - FSPL + G_ant_rx - Pre_rx_loss
%       This value is more interesting than PISO
% Pre_rx_loss: Total loss before the receiver input in the receiver system. Sum of feeder 
%              loss, radome loss, loss due to polarization mismatch, etc.


% ----------------------------------------------------------------------------------------
% Uplink

link_UL = link(gs_transmitter, sat_receiver);  % Perform uplink analysis

[link_UL_PISO_dBW, link_UL_PRI_dBW, link_UL_time] = sigstrength(link_UL);  % Signal power
% [link_UL_EbNo_dB, link_UL_time] = ebno(link_UL);  % EbNo
% link_UL_SNR_dB = link_UL_EbNo_dB + 10 * log10(bitrate / bandwidth);  % SNR?


% ----------------------------------------------------------------------------------------
% Downlink

link_DL = link(sat_transmitter, gs_receiver);  % Perform downlink analysis

[link_DL_PISO_dBW, link_DL_PRI_dBW, link_DL_time] = sigstrength(link_UL);  % Signal power 
% [link_DL_EbNo_dB, link_DL_time] = ebno(link_DL);    % EbNo
% link_DL_SNR_dB = link_DL_EbNo_dB + 10 * log10(bitrate / bandwidth);  % SNR?


% ----------------------------------------------------------------------------------------
% Plot
figure;

subplot(1, 2, 1);
plot(link_UL_time, link_UL_PRI_dBW + 30); 
grid on; grid minor;
xlabel("Simulation time", interpreter="latex");
ylabel("Received signal strength (dBm)", interpreter="latex");
title("Uplink", interpreter="latex");

subplot(1, 2, 2);
plot(link_DL_time, link_DL_PRI_dBW + 30); 
grid on; grid minor;
xlabel("Simulation time", interpreter="latex");
ylabel("Received signal strength (dBm)", interpreter="latex");
title("Downlink", interpreter="latex");

sgtitle("Received signal strength in uplink and downlink", interpreter="latex");



% Plot
% figure;
% 
% subplot(1, 2, 1);
% plot(link_UL_time, link_UL_EbNo_dB); 
% hold on;
% plot(link_UL_time, link_UL_SNR_dB);
% grid on; grid minor;
% xlabel("Simulation time");
% ylabel("Received Eb/No or SNR (dB)");
% title("Uplink");
% legend("$E_b/N_o$", "SNR", interpreter="latex");
% 
% subplot(1, 2, 2);
% plot(link_DL_time, link_DL_EbNo_dB); 
% hold on;
% plot(link_DL_time, link_DL_SNR_dB);
% grid on; grid minor;
% xlabel("Simulation time");
% ylabel("Received Eb/No or SNR (dB)");
% title("Downlink");
% legend("$E_b/N_o$", "SNR", interpreter="latex");
% 
% sgtitle("$E_b/N_o$ and SNR in uplink and downlink", interpreter="latex");





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation 3D visualization

% v = satelliteScenarioViewer(scenario, ShowDetails=true);
% sat.ShowLabel = true;
% gs.ShowLabel = true;
% 
% show(sat);  % Show satellite
% 
% play(scenario);  % Show scenario, play it








