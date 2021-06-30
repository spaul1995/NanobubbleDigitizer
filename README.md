# NanobubbleDigitizer
A MATLAB based algorithm for detection of nanobubble induced current blockage signals observed during resistive pulse sensing in solid-state nanopores.

M1,M2,...M20 are the current vs time graphs for 20 consecutive voltage pulses of 9V.
detectingsignals_array.M is the MATLAB code file which reads M1 to M20 in a loop. Within each loop it detects the bubble blockages and calcultes the blockage duration, waiting time, nucleation current and waiting time.
