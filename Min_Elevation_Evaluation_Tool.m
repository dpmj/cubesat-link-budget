%% Minimum Elevation Angle calculator
% Juan Del Pino Mena - 2023
% PLUTON UPV

% Graphically calculates the received power and SNR by minimum elevation angle. Considers
% only the worst case scenario (minimum elevation, maximum distances).
% This uses a simple FSPL model and a circular orbit

close all;
clearvars;



%% Constants

RT = 6371e3; % Radius of the Earth (m)

epsilon_0 = 8.854188e-12;  % [F/m] Vacuum electrical permitivity
mu_0 = 4 * pi * 1e-7;  % [T x m/A] Vacuum magnetic permitivity
c_0 = 1 / sqrt(epsilon_0 * mu_0);  % [m/s] Wave propagation speed in vacuum

k_B = 1.381e-23;  %  [J/K] Boltzmann's constant



%% Configuration

% Orbit parameters

h = 520e3;  % [m] Height of the satellite orbit avobe the Earth's surface


% Radio parameters

freq = 434e6;  % [Hz] Carrier frequency
BW = 125e3;  % [Hz] Bandwidth
P_tx_dBm = 20;  % [dBm] TX Power
S_rx_dBm = -137;  % [dBm] RX Sensitivity
SNR_limit_dBm = -20;  % [dB] Minimum SNR required

G_sat = 0;  % [dBi] Satellite antenna gain
G_gs = 0;  % [dBi] Ground Station antenna gain

% Noise in the antennas
noise_temp_ant_ul_K = 290;  % [K] For uplink 290 K is considered as the worst case
noise_temp_ant_dl_K = 2340; % [K] For downlink 2340 K for an electromagnetically busy area


% Losses

% These losses must be calculated prior execution. The values specified here are suitable 
% for a 434 MHz, circularly polarized carrier. Pessimistic approximation (complete loss of
% a polarization component due to Faraday rotation, and half power lost due to insertion 
% losses). Atmosphetic losses are low at this frequency.

loss_atmos_dB = 0.5;  % [dB] Atmospheric losses
loss_faraday_dB = 3;  % [dB] Faraday rotation losses
loss_insertion_dB = 1;  % [dB] Insertion losses


% Elevation angle sweep setup
min_elevation_degree = 0:0.5:89;  % [degrees]
min_elevation_rad = (pi / 180) .* min_elevation_degree;  % [rads]



%% Calculus of losses

% [sr] Solid angle of the satellite coverage.
solid_angle_sr = pi/2 - min_elevation_rad ...
                 - asin((RT .* cos(min_elevation_rad)) ./ (RT + h));

% [m] Slant range at minimum elevation (worst case)
max_distance_m = (RT + h) .* (sin(solid_angle_sr) ./ cos(min_elevation_rad));

% FSPL: Free-Space path losses
fspl_linear = ((4 .* pi .* max_distance_m .* freq) ./ (c_0)) .^ 2;
fspl_dB = 10 .* log10(fspl_linear);

% Total losses
loss_total_dB = fspl_dB + loss_atmos_dB + loss_faraday_dB + loss_insertion_dB;


figure();
plot(min_elevation_degree, fspl_dB, '--k');
hold on;
plot(min_elevation_degree, loss_total_dB);
title("Losses as a function of elevation angle")
legend("FSPL Reference", "Total Losses");
xlabel("Minimum elevation degree (º)");
ylabel("Losses (dB)");
grid on; grid minor;



%% Link budget

P_rx_dBm = P_tx_dBm + G_sat + G_gs - loss_total_dB;  % [dB] Received power


% Distinguish between the angles that comply and the ones which do not

P_rx_good = P_rx_dBm(P_rx_dBm >= S_rx_dBm);
P_rx_min_elevation_good = min_elevation_degree(P_rx_dBm >= S_rx_dBm);
P_rx_bad = P_rx_dBm(P_rx_dBm < S_rx_dBm);
P_rx_min_elevation_bad = min_elevation_degree(P_rx_dBm < S_rx_dBm);


figure();
subplot(1,2,1);
plot(P_rx_min_elevation_bad, P_rx_bad, '-r');
hold on;
plot(P_rx_min_elevation_good, P_rx_good, '-g');
yline(S_rx_dBm, "--", sprintf("Sensitivity: %d dBm", S_rx_dBm));
title("Received power as a function of elevation angle");
xlabel("Minimum elevation degree (º)");
ylabel("P_{RX} (dBm)");
grid on; grid minor;

subplot(1,2,2);
plot(P_rx_min_elevation_bad, P_rx_bad - S_rx_dBm, '-r');
hold on;
plot(P_rx_min_elevation_good, P_rx_good - S_rx_dBm, '-g');
yline(0, "--", "Limit of power compliance");
title("Received power margin over Sensitivity, as a function of elevation angle");
xlabel("Minimum elevation degree (º)");
ylabel("P_{RX} - S_{RX} (dB)");
grid on; grid minor;



%% SNR

N_ul_dBm = 10*log10(BW * k_B * noise_temp_ant_ul_K) + 30;  % [dBm] Noise power in uplink
N_dl_dBm = 10*log10(BW * k_B * noise_temp_ant_dl_K) + 30;  % [dBm] Noise power in downlink

SNR_ul_dB = P_rx_dBm - N_ul_dBm;  % [dB] Uplink SNR
SNR_dl_dB = P_rx_dBm - N_dl_dBm;  % [dB] Downlink SNR


% Distinguish between the angles that comply and the ones which do not

SNR_ul_good = SNR_ul_dB(SNR_ul_dB >= SNR_limit_dBm);
SNR_ul_min_elevation_good = min_elevation_degree(SNR_ul_dB >= SNR_limit_dBm);
SNR_ul_bad = SNR_ul_dB(SNR_ul_dB < SNR_limit_dBm);
SNR_ul_min_elevation_bad = min_elevation_degree(SNR_ul_dB < SNR_limit_dBm);

SNR_dl_good = SNR_dl_dB(SNR_dl_dB >= SNR_limit_dBm);
SNR_dl_min_elevation_good = min_elevation_degree(SNR_dl_dB >= SNR_limit_dBm);
SNR_dl_bad = SNR_dl_dB(SNR_dl_dB < SNR_limit_dBm);
SNR_dl_min_elevation_bad = min_elevation_degree(SNR_dl_dB < SNR_limit_dBm);


figure();
subplot(1,2,1);
plot(SNR_ul_min_elevation_bad, SNR_ul_bad, '-r');
hold on;
plot(SNR_ul_min_elevation_good, SNR_ul_good, '-g');
yline(SNR_limit_dBm, "--", sprintf("SNR limit: %d dB", SNR_limit_dBm));
title("Uplink SNR as a function of elevation angle");
xlabel("Minimum elevation degree (º)");
ylabel("SNR at receiver (dB)");
grid on; grid minor;

subplot(1,2,2);
plot(SNR_dl_min_elevation_bad, SNR_dl_bad, '-r');
hold on;
plot(SNR_dl_min_elevation_good, SNR_dl_good, '-g');
yline(SNR_limit_dBm, "--", sprintf("SNR limit: %d dB", SNR_limit_dBm));
title("Downlink SNR as a function of elevation angle");
xlabel("Minimum elevation degree (º)");
ylabel("SNR at receiver (dB)");
grid on; grid minor;














