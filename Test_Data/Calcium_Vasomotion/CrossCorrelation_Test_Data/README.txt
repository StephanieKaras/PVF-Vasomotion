Cross-Correlation Test Data

test_01_calcium_leads_positive.csv
Tests a strong positive relationship in which calcium leads vessel diameter by approximately 1 second.

test_02_diameter_leads_positive.csv
Tests a strong positive relationship in which vessel diameter leads calcium by approximately 1 second.

test_03_synchronized_positive.csv
Tests a strong positive relationship with calcium and vessel diameter occurring at approximately the same time.

test_04_synchronized_negative.csv
Tests a strong inverse relationship with calcium and vessel diameter occurring at approximately the same time.

test_05_uncoupled_noise.csv
Tests signals with no intentionally imposed temporal relationship between calcium and vessel diameter.

test_06_rejected_flat_diameter.csv
Tests the flat-diameter gating condition and should be rejected because vessel diameter has essentially no variation.

test_07_rejected_no_diameter_event.csv
Tests the event-based gating condition and should be rejected because diameter varies but never reaches the required 5-unit event threshold.

Notes
All files use a sampling frequency of 1.727 Hz. Numeric data begin on row 12, with calcium in column F and vessel diameter in column G. Diameter is percent of baseline 100.
