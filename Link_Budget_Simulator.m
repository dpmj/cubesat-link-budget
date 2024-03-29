%% SATELLITE LINK BUDGET ANALYSIS - STANDALONE SCRIPT
% Juan Del Pino Mena
% Version 2, December 2023
% 
% Link budget estimator for LEO satellites. Provides a complete link budget.
% This script computes access intervals in a given simulation time, azimuth, elevation, 
% range, latency, Doppler frequency shift, FSPL losses, atmospheric losses, received power
% and CNR.
%
% REQUIREMENTS
% This program requires MatLab >= R2023a and the satellite communications toolbox.
%
% CHANGELOG
%
% Version 2: Renamed to "standalone" to differenciate between script and function. Solved 
% a bug in the downlink signal strength calculus.
%
% Version 1: First iteration. High-level link budget calculations, Configuration, 3D 
% viewer, access, azimuth, elevation, range, latency, doppler, p618 losses, signal 
% strength, CNR


close all;
clearvars;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script config

show3Dviewer = false;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants

k_boltzmann = 1.380649e-23;  % [J/K] Boltzmann's constant
T_0 = 290; % [K] room temperature for noise calculus, usually 290 K
R_E = 6371;  % [km] Earth's average radius. (not for orbit propagation, only aux calc.)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Configuration parameters

% ----------------------------------------------------------------------------------------
% Common transmission parameters

freq_Hz = 434e6;  % [Hz] Transmission frequency 
bitrate_bps = 293;  % [b/s] Transmission bitrate
symbolrate_Sps = 1/0.032768;  % [Sps] Symbol rate (LoRa: Ts = 2^SF / BW)
bandwidth_Hz = 125e3;  % [Hz] Transmission bandwidth


% ----------------------------------------------------------------------------------------
% Transceiver limits (only used for plots)

rx_limit_sensitivity_gs = -137;  % [dBm] The minimum admitted signal power in the gs
rx_limit_cnr_gs = -20;  % [dB] The minimum SNR / CNR the gs transceiver admits 

rx_limit_sensitivity_sat = -137;  % [dBm] The minimum admitted signal power in the sat
rx_limit_cnr_sat = -20;  % [dB] The minimum SNR / CNR the sat transceiver admits


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

% KNOWN PROBLEM:
% Which System Loss (insertion, connectors, etc) to consider?

% Noise 

rx_nf_gs_dB = 6;  % [dB] System noise figure in RX in the ground station

% KNOWN PROBLEM:
% Receivers noise figure? Should consider the complete receiver chain NF, not only the LNA
% (i.e.: from antenna to the transceiver input pin)

ant_ambient_temp_gs_K = 300;  % [K] Antenna ambient temperature in the ground station
ant_noise_temp_gs_K = 290;  % [K] Antenna noise temperature in the ground station
% The equivalent noise temperature seen at the receive output of the antenna. Describes 
% how much noise an antenna produces in a given environment.

% KNOWN PROBLEM:
% Which antennae ambient (physical) and noise temperature to consider? 
% Physical temperature may vary a lot, especially in space. Need to consider the worst 
% case. For uplink 290 K is considered as the worst case (antenna points upwards), and for
% downlink is 2340 K (antenna points to Earth). The temperature depends on directivity

% G/T (Gain-to-Noise-Temperature Ratio, GNTR) relates the receive antenna gain and the 
% system noise temperature.
% G/T = GR − Nf − 10 * log10( T0 + (Ta − T0) * 10 ^ (−0.1 * Nf) )
% GR is the receive antenna gain in dBi
% Nf is the noise figure in dB
% T0 is the ambient temperature in degrees Kelvin
% Ta is the antenna temperature in degrees Kelvin
% https://www.mathworks.com/help/satcom/ug/nb-iot-ntn-link-budget-analysis.html

rx_GNTR_gs_dB_K = tx_ant_gain_gs_dB - rx_nf_gs_dB - 10 * log10(ant_ambient_temp_gs_K + ...
    (ant_noise_temp_gs_K - ant_ambient_temp_gs_K) * 10^(-0.1 * rx_nf_gs_dB));  % [dB/K]


% ----------------------------------------------------------------------------------------
% Satellite

% Power, gain and losses

tx_power_sat_dBm = 33;  % [dBm] Transmission power in the satellite tx
tx_ant_gain_sat_dB = 0;  % [dBi] Antenna Gain in the sat tx
tx_loss_sat_dB = 1;  % [dB] System power loss in TX in the satellite
rx_loss_sat_dB = 1;  % [dB] System power loss in RX in the satellite

% Noise 

rx_nf_sat_dB = 3;  % [dB] System noise figure in RX in the ground station

ant_ambient_temp_sat_K = 350;  % [K] Antenna ambient temperature in the satellite
ant_noise_temp_sat_K = 2340;  % [K] Antenna noise temperature in the satellite

% [dB/K] G/T in the sat rx
rx_GNTR_sat_dB_K = tx_ant_gain_sat_dB - rx_nf_sat_dB ... 
    - 10 * log10(ant_ambient_temp_sat_K ...
    + (ant_noise_temp_sat_K - ant_ambient_temp_sat_K) * 10^(-0.1 * rx_nf_sat_dB)); 


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
%% Compute access to the satellite from the ground station

ac_sat_gs = access(sat, gs);  % Compute access to the satellite
ac_sat_intervals = accessIntervals(ac_sat_gs);  % Access intervals from GS

[ac_sat_plot_data, ac_sat_time] = accessStatus(ac_sat_gs);  % Plot the access intervals

ac_sat_plot_data = double(ac_sat_plot_data);  % Casting for the operation below
% ac_sat_plot_data(ac_sat_plot_data == 0) = NaN;  % Filter out of sight intervals


% ----------------------------------------------------------------------------------------
% Output

fprintf("\n-----------------------------------------\n" + ...
    "ACCESS INTERVALS:\n\n");
disp(ac_sat_intervals);


% ----------------------------------------------------------------------------------------
% Plot points in time where there is access to the satellite.

figure();
plot(ac_sat_time, ac_sat_plot_data, '.-');
title("Observation intervals: line-of-sight with $\varepsilon >= \varepsilon_{min}$" + ...
      ". [1: access; 0: no access]", interpreter="latex");
xlabel("Simulation Time", interpreter="latex");
ylabel("Access status", interpreter="latex");
ylim([-0.5, 1.5]);
grid on; grid minor;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add TX and RX to the satellite and the Ground Station

% ----------------------------------------------------------------------------------------
% Ground Station transceivers

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
% Satellite transceivers

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
%% Azimuth, elevation and range between satellite and ground station

% Compute azimuth, elevation, range
[azimuth_sat_gs_deg, elevation_sat_gs_deg, range_sat_gs_m] = aer(gs, sat);  % [deg,deg,m]

range_sat_gs_full_m = range_sat_gs_m;  % Save a copy of the complete vector for CNR calc.
aer_sat_gs_time = ac_sat_time;  % time vector for plotting

% Discard data when not in line of sight
azimuth_sat_gs_deg(~logical(ac_sat_plot_data)) = NaN;  
elevation_sat_gs_deg(~logical(ac_sat_plot_data)) = NaN;
range_sat_gs_m(~logical(ac_sat_plot_data)) = NaN; 


% ----------------------------------------------------------------------------------------
% Plot azimuth, elevation, range

figure;

subplot(3, 1, 1);
plot(aer_sat_gs_time, azimuth_sat_gs_deg); 
grid on; grid minor;
xlabel("Simulation time", interpreter="latex");
ylabel("Azimuth (degrees)", interpreter="latex");
title("Ground Station Azimuth", interpreter="latex");

subplot(3, 1, 2);
plot(aer_sat_gs_time, elevation_sat_gs_deg); 
grid on; grid minor;
xlabel("Simulation time", interpreter="latex");
ylabel("Elevation (degrees)", interpreter="latex");
title("Ground Station Elevation", interpreter="latex");

subplot(3, 1, 3);
plot(aer_sat_gs_time, range_sat_gs_m); 
grid on; grid minor;
xlabel("Simulation time", interpreter="latex");
ylabel("Distance (m)", interpreter="latex");
title("Range between GS and Sat", interpreter="latex");



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Latency and Doppler Shift
% These functions are only available on MatLab R2023a onwards
% Dynamic - computes latency and doppler shift over the complete simulation time

[latency_delay, latency_time] = latency(sat, gs);
[doppler_fshift, doppler_time, doppler_info] = dopplershift(sat, gs, Frequency=freq_Hz);


% ----------------------------------------------------------------------------------------
% Plot

figure();
subplot(2, 1, 1);
plot(latency_time, latency_delay(1, :) .* 1e3);  % in milliseconds
xlim([latency_time(1), latency_time(end)]);
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

[p618_atm_loss_dB, p618_xpol_discr_dB, p618_temp_sky_K] = p618PropagationLosses(p618cfg);

% p618_atm_loss_dB.Ag - Gaseous attenuation (dB)
% p618_atm_loss_dB.Ac - Cloud and fog attenuation (dB)
% p618_atm_loss_dB.Ar - Rain attenuation (dB)
% p618_atm_loss_dB.As - Attenuation due to tropospheric scintillation (dB)
% p618_atm_loss_dB.At - Total atmospheric attenuation (dB)
% p618_xpol_discr_dB  - Cross-polarization discrimination (dB) not exceeded for the 
%                    percentage of the RainAnnualExceedance.
% p618_temp_sky_K     - Sky noise temperature (K) at the ground station antenna.


% ----------------------------------------------------------------------------------------
% Polarization loss

p618_pol_loss = 3;  % [dB] Polarization loss. Worst-case: 3 dB


% ----------------------------------------------------------------------------------------
% Total losses

p618_loss_total = p618_atm_loss_dB.At + p618_pol_loss;


% ----------------------------------------------------------------------------------------
% Output

fprintf("\n\n-----------------------------------------\nATMOSPHERIC LOSSES:\n");
fprintf("\nGaseous attenuation = %.3f dB", p618_atm_loss_dB.Ag);
fprintf("\nCloud and fog attenuation = %.3f dB", p618_atm_loss_dB.Ac);
fprintf("\nRain attenuation = %.3f dB", p618_atm_loss_dB.Ar);
fprintf("\nAttenuation due to tropospheric scintillation = %.3f dB", p618_atm_loss_dB.As);
fprintf("\nTotal atmospheric attenuation = %.3f dB", p618_atm_loss_dB.At);
fprintf("\nCross-polarization discrimination = %.3f dB", p618_xpol_discr_dB);
fprintf("\nSky noise temperature = %.3f K", p618_temp_sky_K);
fprintf("\nPolarization loss (MANUAL PICK) = %.3f dB", p618_pol_loss);
fprintf("\nTOTAL LOSSES = Total att. + Pol. loss = %.3f dB\n\n", p618_loss_total);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Link budget analysis: Received Power

% EIRP: Equivalent Isotropic Radiated Power: EIRP = P_tx + G_ant_tx
% Pre_rx_loss: Total loss before the receiver input in the receiver system. Sum of feeder 
%              loss, radome loss, loss due to polarization mismatch, etc.

% PISO: Received isotropic power in the antenna: PISO_rx = EIRP - FSPL 
% PRI:  Received power at receiver input: PRI_rx = EIRP_tx - FSPL + G_ant_rx - Pre_rx_loss

% The PRI result is the most interesting. But since the PRI does not account for the p618
% losses, we will substract them from the received power


% ----------------------------------------------------------------------------------------
% Uplink

link_UL = link(gs_transmitter, sat_receiver);  % Create link
[link_UL_PISO_dBW, link_UL_PRI_dBW, link_UL_time] = sigstrength(link_UL);  % Signal power

link_UL_rxpower_dBm = link_UL_PRI_dBW + 30 - p618_loss_total;  % [dBm] UL received power


% ----------------------------------------------------------------------------------------
% Downlink

link_DL = link(sat_transmitter, gs_receiver);  % Create link
[link_DL_PISO_dBW, link_DL_PRI_dBW, link_DL_time] = sigstrength(link_DL);  % Signal power 

link_DL_rxpower_dBm = link_DL_PRI_dBW + 30 - p618_loss_total;  % [dBm] DL received power


% ----------------------------------------------------------------------------------------
% Plot

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



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CNR: Carrier-to-noise ratio for configured satellite link budget parameters (SNR)
% Static by default, forced to be dynamic (SNR computed in every simulation point)

% ----------------------------------------------------------------------------------------
% Computes CNR in every simulation point

N_points_sim = length(ac_sat_plot_data);

ul_CNR = zeros(N_points_sim, 1);
dl_CNR = zeros(N_points_sim, 1);
% ul_CNR_info = zeros(N_points_sim, 1);  % Cannot allocate due to type not being double
% dl_CNR_info = zeros(N_points_sim, 1);

for i = 1:1:N_points_sim

    % UPLINK: Ground Station -> Satellite
    ul_CNR_cfg = satelliteCNRConfig(...
        TransmitterPower = tx_power_gs_dBm - 30, ...  % [dBW] Transmit power 
        TransmitterSystemLoss = tx_loss_gs_dB, ...  % [dB] Losses in the transmitter
        TransmitterAntennaGain = tx_ant_gain_gs_dB, ...  % [dBi] Antenna gain
        Distance = range_sat_gs_full_m(i) * 1e-3, ...  % [km] Distance. Varies.
        Frequency = freq_Hz * 1e-9, ...  % [GHz] Center frequency
        MiscellaneousLoss = p618_loss_total, ...  % [dB] Misc attenuation
        GainToNoiseTemperatureRatio = rx_GNTR_sat_dB_K, ...  % [dB/K] Receiver G/T ratio
        ReceiverSystemLoss = rx_loss_sat_dB, ...  % [dB] Receiver system loss
        BitRate = bitrate_bps * 1e-6, ...  % [Mbps] Bit rate
        SymbolRate = symbolrate_Sps * 1e-6, ...  % [MSps] Symbol rate
        Bandwidth = bandwidth_Hz * 1e-6 ...  % [MHz] Bandwidth
    );
    % [ul_CNR(i), ul_CNR_info(i)] = satelliteCNR(ul_CNR_cfg);
    [ul_CNR(i), ~] = satelliteCNR(ul_CNR_cfg);
    
    % DOWNLINK: Satellite -> Ground Station
    dl_CNR_cfg = satelliteCNRConfig(...
        TransmitterPower = tx_power_sat_dBm - 30, ...  % [dBW] Transmit power 
        TransmitterSystemLoss = tx_loss_sat_dB, ...  % [dB] Losses in the transmitter
        TransmitterAntennaGain = tx_ant_gain_sat_dB, ...  % [dBi] Antenna gain
        Distance = range_sat_gs_full_m(i) * 1e-3, ...  % [km] Distance. Varies.
        Frequency = freq_Hz * 1e-9, ...  % [GHz] Center frequency
        MiscellaneousLoss = p618_loss_total, ...  % [dB] Misc attenuation
        GainToNoiseTemperatureRatio = rx_GNTR_gs_dB_K, ...  % [dB/K] Receiver G/T ratio
        ReceiverSystemLoss = rx_loss_gs_dB, ...  % [dB] Receiver system loss
        BitRate = bitrate_bps * 1e-6, ...  % [Mbps] Bit rate
        SymbolRate = symbolrate_Sps * 1e-6, ...  % [MSps] Symbol rate
        Bandwidth = bandwidth_Hz * 1e-6 ...  % [MHz] Bandwidth
    );
    % [dl_CNR(i), dl_CNR_info(i)] = satelliteCNR(dl_CNR_cfg);
    [dl_CNR(i), ~] = satelliteCNR(dl_CNR_cfg);

end

% Filter data out of line of sight
ul_CNR(~logical(ac_sat_plot_data)) = NaN;
dl_CNR(~logical(ac_sat_plot_data)) = NaN;

cnr_time = ac_sat_time;  % time vector for plotting


% ----------------------------------------------------------------------------------------
% Plot

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



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation 3D visualization
% This representation consumes time and their windows can be an annoyance.

if show3Dviewer == true  % only if requested

    v = satelliteScenarioViewer(scenario, ShowDetails=true);
    sat.ShowLabel = true;
    gs.ShowLabel = true;
    
    show(sat);  % Show satellite
    
    play(scenario);  % Show scenario, play it

end


