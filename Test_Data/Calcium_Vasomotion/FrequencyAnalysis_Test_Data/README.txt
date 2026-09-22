Calcium Vasomotion Frequency Analysis - Test Data

test_01_0p06Hz.csv
Tests calcium and vessel diameter oscillating together at 0.06 Hz, near the lower end of the vasomotor analysis range.

test_02_0p10Hz.csv
Tests calcium and vessel diameter oscillating together at 0.10 Hz, the approximate 10-second vasomotor rhythm relevant to this analysis.

test_03_0p18Hz.csv
Tests calcium and vessel diameter oscillating together at 0.18 Hz, near the upper end of the analysis range.

test_04_mixed_frequencies.csv
Tests signals containing both 0.08 Hz and 0.14 Hz components to confirm the analysis handles signals with more than one frequency.

test_05_different_amplitudes_0p10Hz.csv
Tests signals at the same 0.10 Hz frequency but with different amplitudes to check frequency, bandpower, and Hilbert amplitude measurements.

test_06_short_recording.csv
Tests a 70-second 0.10 Hz recording to confirm the script uses its shorter Welch-window setting for recordings under 90 seconds.

Notes
The files use a sampling frequency of 1.727 Hz. Numeric data begin on row 12, with calcium in column F and vessel diameter in column G. Diameter is percent of baseline 100.
