# Note on Number of Realizations

The Adams, Larsen, and Pomraning paper cites a resolution of $10^5$ resolutions for each benchmarking scheme. I have found that this is not enough resolutions to actually correctly model the results for the low-thickness elements. Simply changing the initiating random number seed changes the calculated parameters for only $10^5$ resolutions.

Origin | Reflection | Transmission
--- | --- | ---
Adams Published | 0.0491 | 0.9331
Calculated Seed 123456 | 0.0049970 | 0.93048
Calculated Seed 1234567 | 0.0041057 | 0.93157
Calculated Seed 12345678 | 0.0038978 | 0.93238

It can be seen that the marked difference in total cross-section between the two materials is the cause of some instability in the initial reflection term of the source on the left wall. The material averaging tends to converge much faster for the transmission term on the right wall.
