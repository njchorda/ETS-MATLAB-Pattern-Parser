# ETSDataFormat3D

A MATLAB function for parsing and loading 3D antenna pattern data from ETS-Lindgren CSV exports into structured arrays.

## Overview

`ETSDataFormat3D` reads a 3D antenna measurement CSV file (as exported by ETS-Lindgren systems) and returns antenna pattern data organized by polarization (phi, theta, and/or total power), along with the corresponding angle and frequency vectors. To speed up repeated loads, the function automatically saves a `.mat` cache file on first use and loads from it on subsequent calls.

## Usage

```matlab
% Four-output form — returns whichever single polarization is available
[D, phi, theta, freqs] = ETSDataFormat3D(filename)

% Six-output form — returns both polarizations AND vector sum (total gain)
[D_total, D_phi, D_theta, phi, theta, freqs] = ETSDataFormat3D(filename)
```

## Caching Behavior

The function uses 'readcell()' to fetch the data from the csv file. For large datasets this operation can be very slow and use a lot of RAM.

On the first call, the function parses the CSV and saves a `.mat` file in the same directory with the same base filename. On subsequent calls, if the `.mat` file is detected it is loaded directly, significantly reducing load time for large datasets. To force a re-parse, simply delete the `.mat` file.

### Output Arguments

| Variable | Size | Description |
|---|---|---|
| `D` / `D_total` | `numPhi × numTheta × numFreqs` | Total power directivity array |
| `D_phi` | `numPhi × numTheta × numFreqs` | Phi-polarization directivity array |
| `D_theta` | `numPhi × numTheta × numFreqs` | Theta-polarization directivity array |
| `phi` | `numPhi × 1` | Phi angles (degrees) |
| `theta` | `1 × numTheta` | Theta angles (degrees) |
| `freqs` | `1 × numFreqs` | Frequencies (Hz) |

### Input Arguments

| Argument | Description |
|---|---|
| `filename` | Full or relative path to the ETS-Lindgren `.csv` measurement file |

## Examples

```matlab
% Load total power only
[D, phi, theta, freqs] = ETSDataFormat3D('myAntenna_3D.csv');

% Plot pattern at first frequency
figure;
imagesc(theta, phi, D(:,:,1));
xlabel('Theta (deg)'); ylabel('Phi (deg)');
title('Total Power Pattern');
colorbar;

% Load all polarizations
[D_total, D_phi, D_theta, phi, theta, freqs] = ETSDataFormat3D('myAntenna_3D.csv');
```

## Notes

- Frequencies in the raw CSV are assumed to be in MHz and are automatically converted to Hz on output.
- The six-output form requires that all three polarizations (total, phi, theta) are present in the CSV. An error is thrown if any are missing.
- The four-output form returns whichever single polarization is available (`D_total` takes priority, followed by `D_vert`, then `D_horiz`).

## Requirements

- MATLAB R2019b or later (uses `readcell`)
- No additional toolboxes required
