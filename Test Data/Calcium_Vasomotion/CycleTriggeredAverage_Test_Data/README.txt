Cycle-Triggered Average Test Data

test_01_calcium_peak_at_diameter_peak.csv
Tests a calcium oscillation whose peak occurs at approximately the same time as the vessel diameter peak.

test_02_calcium_peak_before_diameter_peak.csv
Tests a calcium oscillation whose peak occurs approximately 2 seconds before the vessel diameter peak.

test_03_calcium_peak_after_diameter_peak.csv
Tests a calcium oscillation whose peak occurs approximately 2 seconds after the vessel diameter peak.

test_04_calcium_peak_at_diameter_trough.csv
Tests an inverse-phase relationship in which the calcium peak occurs approximately at the vessel diameter trough.

test_05_weak_calcium_oscillation.csv
Tests a lower-amplitude calcium oscillation while retaining a clear vasomotor oscillation in vessel diameter.

test_06_rejected_short_recording.csv
Tests the minimum-duration check and should be skipped because the recording is shorter than 50 seconds.

test_07_poor_vasomotor_band_content.csv
Tests a diameter signal dominated by a frequency below the specified vasomotor analysis band.

test_08_flat_diameter.csv
Tests an essentially flat diameter signal that should not provide reliable vasomotor cycles for cycle-triggered averaging.

Notes
The full-length files are 180 seconds long and use a sampling frequency of 1.951 Hz. Numeric data begin on row 12, with time in column E, calcium in column F, and vessel diameter in column G. Diameter is percent of baseline 100.
