%% 
%% OPT 212: Optical Communication Pulse Trains Coding Exercises
% Student version
% 
% These exercises accompany the Optical Communication Pulse Trains active script. 
% They focus on sampled data sequences, additive Gaussian receiver noise, first-order 
% electrical filtering, and eye diagrams.
% 
% Complete the TODO blocks in each work area. The companion solution file is 
% Week11_optical_communication_pulse_trains_solutions.m.
% 
% Setup — clear removes workspace variables, close all closes prior figures, 
% clc clears the Command Window, and format compact reduces printed whitespace. 
% The groot settings set readable default font sizes for the plots.

clear
close all
clc
format compact
set(groot,'DefaultAxesFontSize',14,'DefaultTextFontSize',14, ...
    'DefaultLegendFontSize',11,'DefaultColorbarFontSize',14);
%% Exercise 1: NRZ Optical Power Waveform
% A baseband optical communication simulation usually represents the slowly 
% varying optical power envelope rather than the optical carrier oscillation itself. 
% For this exercise, use a wavelength of 1550 nm as the system context, but simulate 
% the normalized optical power envelope.
% 
% Use a bit rate of |40 Gb/s|, a |48| bit pseudorandom sequence, and |128| samples 
% per bit. Since |1 Gb/s| has the same numerical period as |1 GHz|, the bit period 
% in picoseconds is
% 
% $$T_b={1000\over R_b},$$
% 
% when $R_b$ is in Gb/s and $T_b$ is in ps. Use |round(rand(\[1,48\]))| or a 
% similar command to generate the binary data.
% 
% Tasks:
% 
% 1. Convert the bit rate to a bit period in ps.
% 
% 2. Create a uniformly sampled time vector for all 48 bits.
% 
% 3. Generate an NRZ waveform where each bit slot is completely filled by either 
% the zero level or the one level.
% 
% 4. Plot the first 8 bits of the NRZ waveform.
%% Exercise 1 Work Area
% TODO: Add your code for Exercise 1 here.
%% Worked solution 1
% Step 1 — Specify the physical context and sampling. The wavelength labels 
% the optical system; this baseband calculation models normalized power, so it 
% does not resolve the optical carrier. NRZ means non-return-to-zero. The waveform 
% uses one-level power 1 and zero-level power 0, holding each level for the full 
% bit slot.

wavelength_nm = 1550 ;  % wavelength in nm 
bit_rate_Gbps =  40;
number_of_bits = 48 ;  % total number of bits in waveform
samples_per_bit = 128 ; 

%% 
% Step 2 — Convert units and construct the grid. At 40 Gb/s, each bit lasts 
% 25 ps. Dividing by 128 gives dt = 0.1953125 ps. There are 48 × 128 = 6144 samples. 
% Starting sample numbers at zero places the first sample at t = 0; ending at 
% N−1 avoids duplicating the endpoint of the 1200 ps record.

bit_period_ps = 1000/bit_rate_Gbps;
dt_ps = bit_period_ps/samples_per_bit;
time_ps = [0:number_of_bits*samples_per_bit]*dt_ps;
%% 
% Step 3 — Generate repeatable random data and form the NRZ waveform. rng fixes 
% the random-number seed so a rerun gives the same example. rand generates uniform 
% numbers, and round maps them to equally probable zeros and ones. The make_nrz_waveform 
% helper holds each bit value for its entire slot. Consecutive one bits therefore 
% form an uninterrupted high level; consecutive zero bits form an uninterrupted 
% low level.

rng(2129)  % Optional -- fixing the random number seed guarantees the
% same waveform every time
bits = round(rand([number_of_bits 1])) ;  % create a random sequence of bits
bit_period_ps = 1000/bit_rate_Gbps;  % width (in ps) of every bit 
% % floor maps time to a zero-based bit slot; +1 converts to MATLAB indexing.
bit_index = floor(time_ps/bit_period_ps) + 1;
% % Guard against time samples outside the supplied bit sequence.
valid = bit_index >= 1 & bit_index <= numel(bits);
nrz_power = zeros(size(time_ps));
% % Repeat each bit value at every sample assigned to its slot.
nrz_power(valid) = bits(bit_index(valid))   ;

%% 
% Step 4 — Select the display interval without shortening the simulation. The 
% logical mask selects t < 200 ps, the first eight bit slots. 

first_bits_to_show = 8;
display_mask = time_ps < first_bits_to_show*bit_period_ps;
%% 
% Step 5 — Inspect the NRZ waveform. stairs emphasizes the piecewise-constant 
% samples, and dotted vertical lines mark bit boundaries. The vertical limits 
% leave space around the zero and one levels. Transitions occur only when adjacent 
% bits differ; the signal stays high across consecutive one bits.

figure('Name','Solution 1: NRZ waveform')
plot(time_ps(display_mask),nrz_power(display_mask),'k-', ...
     'LineWidth',1.4)
% 
grid on
xlabel('time (ps)','FontSize',14)
ylabel('power','FontSize',14)
ylim([-0.05 1.15])
title('NRZ waveform: each bit fills the full bit slot');
%% 
% Step 6 — Report the parameters and transmitted sequence. table labels each 
% quantity with its units. Converting bits to strings and joining them produces 
% a readable 48-character sequence. Omitting the final semicolon displays the 
% table in the Live Editor.

waveform_summary = table(wavelength_nm,bit_rate_Gbps,bit_period_ps,dt_ps, ...
     number_of_bits,samples_per_bit, ...
     string(join(string(bits),'')), ...
     'VariableNames',{'Wavelength_nm','BitRate_Gbps','BitPeriod_ps', ...
     'SampleSpacing_ps','NumberOfBits','SamplesPerBit', ...
     'BitSequence'})
%% Exercise 2: Additive Gaussian Noise and Decisions
% A simple receiver model adds Gaussian noise to the normalized power trace:
% 
% $$P_{\it{meas}}(t)=P(t)+\sigma n(t),$$
% 
% where $n(t)$ is a zero-mean random sequence with unit standard deviation. 
% Use the NRZ waveform from Exercise 1, noise RMS |sigma = 0.18|, and a decision 
% threshold of |0.5|. Sample the waveform once per bit at the bit center.
% 
% For equal zero and one levels at 0 and 1, equal Gaussian noise on both levels, 
% and a threshold halfway between them, the approximate analytical bit-error probability 
% is
% 
% $$P_e=Q\!\left({0.5\over \sigma}\right) ={1\over 2}\it{erfc}\!\left({0.5\over 
% \sqrt{2}\sigma}\right).$$
% 
% Tasks:
% 
% 1. Add Gaussian noise to the NRZ waveform.
% 
% 2. Sample the noisy waveform at the center of each bit.
% 
% 3.  Plot the clean waveform, noisy waveform, and bit-center samples for the 
% first 8 bits.
%% Exercise 2 Work Area
% TODO: Add your code for Exercise 2 here.
%% Worked solution 2
% Step 1 — Set the receiver model. The noise standard deviation is 0.18 in normalized 
% power units, and the decision threshold is halfway between 0 and 1. Run the 
% script from the top: this exercise reuses the NRZ waveform and sampling grid 
% from Exercise 1.

noise_rms =  0.18   ;
decision_threshold = 0.5  ;
%% 
% Step 2 — Add noise, sample, and decide. randn(size(nrz_power)) produces one 
% independent unit-variance Gaussian noise sample per waveform sample; multiplying 
% by noise_rms sets its standard deviation. The seed makes this example reproducible. 
% Negative measured values are allowed because the model represents an additive 
% receiver signal, not a constraint on physical optical power. 

rng(2130)
noisy_nrz_power = nrz_power + noise_rms*randn(size(nrz_power));

%% 
% Step 3 — Inspect the first eight decisions. plot shows the noisy trace, stairs 
% overlays the clean reference, and circular markers show only the samples actually 
% used by the detector. hold on retains the existing curves. The horizontal line 
% is the threshold; padded limits keep noise excursions visible. A sample above 
% the threshold is classified as one, regardless of the waveform between sampling 
% instants.

figure
hold on
plot(time_ps(display_mask),nrz_power(display_mask), '--r','LineWidth',1.4)
plot(time_ps(display_mask),noisy_nrz_power(display_mask), ...
     'Color',[0.55 0.55 0.55],'LineWidth',1.0)
 hold off
 grid on
 xlabel('time (ps)','FontSize',14)
 ylabel('normalized power','FontSize',14)
 
 legend('noisy NRZ trace','noise-free NRZ reference')
 title(sprintf(['Noisy NRZ trace, sigma = %.2f, threshold = ', ...
     '%.2f'],noise_rms,decision_threshold));
%% 
% 


%%