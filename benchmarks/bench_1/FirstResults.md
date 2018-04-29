# Results of First Benchmark

## Parameters Given

Parameter | Material 1 | Material 2
--- | --- | ---
$\sigma_t$ | 10/99 | 100/11
$c$ | 0.00 | 1.00
$\lambda$ | 99/100 | 11/100

## Results

The exact values are those resulting from $1 \times 10^8$ realizations, with the exception of the 10.0 thickness data, which is calculated for $1 \times 10^6$ realizations for expediency.

The model values are those resulting from a convergence to $1 \times 10^{-9}$ percent error between iterations.

### Thickness = 0.1

Parameter | Published Exact | Published Model | Calculated Exact | Calculated Model
--- | --- | --- | --- | ---
Reflection | 0.0491 | 0.0479 | 0.0049543 | 0.0048039
Transmission | 0.9331 | 0.9343 | 0.9325667 | 0.9371273

### Thickness = 1.0

Parameter | Published Exact | Published Model | Calculated Exact | Calculated Model
--- | --- | --- | --- | ---
Reflection | 0.2495 | 0.2187 | 0.2426435 | 0.2193373
Transmission | 0.5950 | 0.6254 | 0.5887067 | 0.6273042

### Thickness = 10.0

Alpha = 3.2712594

Parameter | Published Exact | Published Model | Calculated Exact | Calculated Model | Alpha Closure | Atomic Mix
--- | --- | --- | --- | --- | --- | ---
Reflection | 0.4342 | 0.3760 | 0.4365407 | 0.3771022 | 0.4482795 | 0.4981935
Transmission | 0.0146 | 0.0259 | 0.0146235 | 0.0259773 | 0.0104524 | 0.0047435

Monte Carlo Data:
Parameter | Atomic Mix (1E7 p) | Closure (1E7 p) | Alpha (1E7 p) | Exact (1E4 p, 5E5 r)
--- | --- | --- | --- | ---
Reflection | 0.4958776 | 0.3783326 | 0.4485059 | 0.4367145
Transmission | 0.0047801 | 0.0264183 | 0.0106775 | 0.0148474
