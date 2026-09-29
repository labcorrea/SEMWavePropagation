function waveform = sem_hann_sine(frequency, cycles, amplitude, startTime, samplesPerCycle)
%SEM_HANN_SINE Finite sine burst multiplied by a Hann window.
% The returned structure has sampled time and values arrays for a source.
if nargin < 3 || isempty(amplitude), amplitude = 1; end
if nargin < 4 || isempty(startTime), startTime = 0; end
if nargin < 5 || isempty(samplesPerCycle), samplesPerCycle = 40; end
validateattributes(frequency, {'numeric'}, {'scalar','real','positive','finite'});
validateattributes(cycles, {'numeric'}, {'scalar','real','positive','finite'});
validateattributes(amplitude, {'numeric'}, {'scalar','real','finite'});
validateattributes(startTime, {'numeric'}, {'scalar','real','finite'});
validateattributes(samplesPerCycle, {'numeric'}, {'scalar','integer','>=',2});
duration = cycles/frequency;
tau = linspace(0,duration,ceil(cycles*samplesPerCycle)+1);
window = 0.5*(1-cos(2*pi*tau/duration));
values = amplitude*sin(2*pi*frequency*tau).*window;
values([1,end]) = 0;
waveform = struct('time',startTime+tau,'values',values);
end
