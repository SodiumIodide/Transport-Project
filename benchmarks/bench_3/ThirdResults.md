# Results of Third Benchmark

## Parameters Given

Parameter | Material 1 | Material 2
--- | --- | ---
$\sigma_t$ | 10/99 | 100/11
$c$ | 0.90 | 0.90
$\lambda$ | 99/100 | 11/100

## Results

The exact values are those resulting from $1 \times 10^8$ realizations, with the exception of the 10.0 thickness data, which is calculated for $1 \times 10^6$ realizations for expediency.

The model values are those resulting from a convergence to $1 \times 10^{-9}$ percent error between iterations.

### Thickness = 0.1

Parameter | Published Exact | Published Model | Calculated Exact | Calculated Model
--- | --- | --- | --- | ---
Reflection | 0.0480 | 0.0473 | 0.0488849 | 0.0475541
Transmission | 0.9341 | 0.9344 | 0.9328820 | 0.9372570

### Thickness = 1.0

Parameter | Published Exact | Published Model | Calculated Exact | Calculated Model
--- | --- | --- | --- | ---
Reflection | 0.2563 | 0.2178 | 0.2520195 | 0.2184895
Transmission | 0.5985 | 0.6267 | 0.5938456 | 0.6286463

### Thickness = 10.0

Alpha = 1.7392527

Parameter | Published Exact | Published Model | Calculated Exact | Calculated Model | Alpha Closure | Atomic Mix
--- | --- | --- | --- | --- | --- | ---
Reflection | 0.4785 | 0.3707 | 0.4799376 | 0.3717843 | 0.4036980 | 0.4807157
Transmission | 0.0159 | 0.0237 | 0.0159978 | 0.0238105 | 0.0140277 | 0.0038574

Paramter | $10^6$ Realizations | $10^7$ Realizations
--- | --- | ---
Reflection | 0.4799704 | 0.4799719
Transmission | 0.0157653 | 0.0158591

Monte Carlo Data:
Parameter | Atomic Mix (1E7 p) | Closure (1E7 p) | Alpha (1E7 p) | Exact (1E4 p, 5E5 r)
--- | --- | --- | --- | ---
Reflection | 0.4785959 | 0.3700725 | 0.4015677 | 0.4772735
Transmission | 0.0038917 | 0.0238894 | 0.0139971 | 0.0160364
