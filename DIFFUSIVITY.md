# Diffusivity Calculation Guide

This guide describes how the standalone diffusivity utility calculates the value printed by:

```bash
./cal_diffusivity example.input out.dat
```

The result is printed as one number, for example:

```text
3.04640e-08
```

## Reported Unit

The reported diffusivity has units of:

```text
m^2/s
```

This follows from the code:

```text
D = 0.5 * k * L^2
```

where `k` has units of `1/s` and `L` is a hopping length in meters.

## Files Used

The utility uses two files:

```bash
./cal_diffusivity input_file output_path_file
```

For the example:

```bash
./cal_diffusivity example.input out.dat
```

It reads:

1. Cell size, cell angles, molecular mass, temperature, direction, and number of string points from `example.input`
2. The optimized string path and external potential values from `out.dat`

## Required Output File Format

The path file must have one row per string point. Each row has seven columns:

```text
frac_a frac_b frac_c alpha beta gamma V
```

| Column | Meaning |
| --- | --- |
| `frac_a frac_b frac_c` | Adsorbate center-of-mass position in fractional coordinates |
| `alpha beta gamma` | Euler angles in radians |
| `V` | External potential along the path, in Kelvin energy units |

The number of rows must match `#_of_points` in the input file.

## Step 1: Convert the Path to Meters

The path positions in `out.dat` are fractional coordinates. The code converts each point to Cartesian coordinates using the unit-cell lengths and angles:

```text
frac_a frac_b frac_c  ->  x y z
```

The cell lengths in the input are in Angstrom, so the code converts Cartesian coordinates to meters:

```text
x_m = x_A * 1e-10
y_m = y_A * 1e-10
z_m = z_A * 1e-10
```

## Step 2: Calculate the String Arc Length

The code builds a cumulative distance coordinate `s` along the string:

```text
s[0] = 0
s[i] = s[i-1] + distance(point_i, point_{i-1})
```

Each segment distance is calculated in 3D Cartesian space:

```text
distance = sqrt((dx)^2 + (dy)^2 + (dz)^2)
```

Because the Cartesian coordinates were converted to meters, `s` is also in meters.

## Step 3: Choose the Hopping Length

The hopping length `L` is chosen from the input `Direction`:

| Direction | Hopping length used |
| --- | --- |
| `1` | `La * 1e-10` meters |
| `2` | `Lb * 1e-10` meters |
| `3` | `Lc * 1e-10` meters |

So the direction affects the final diffusivity through `L^2`.

## Step 4: Convert Molecular Mass

The molecular mass is read from the input line:

```text
cutoff(A) FH_signal Mass(g/mol) Tempearture(K)  Running_steps
```

The code converts mass from `g/mol` to `kg/molecule`:

```text
m = Mass(g/mol) / 1000 / N_A
```

where:

```text
N_A = 6.02214076e23
```

## Step 5: Calculate the Thermal Velocity Factor

The velocity-like prefactor is:

```text
sqrt(k_B * T / (2 * pi * m))
```

where:

| Symbol | Meaning | Unit |
| --- | --- | --- |
| `k_B` | Boltzmann constant | J/K |
| `T` | Temperature from the input file | K |
| `m` | Adsorbate molecular mass | kg/molecule |

This term has units of `m/s`.

## Step 6: Use the Potential Profile Along the String

For each string point, the utility calculates:

```text
exp(-V[i] / T)
```

This works because `V` is treated as an energy in Kelvin units. Therefore `V/T` is dimensionless.

Then it integrates this Boltzmann factor along the path using the trapezoidal rule:

```text
integral = trapz(s, exp(-V/T))
```

Since `s` is in meters, the integral has units of meters.

## Step 7: Calculate the Rate Constant

The rate-like quantity is:

```text
k = sqrt(k_B * T / (2 * pi * m)) * exp(-max(V) / T) / trapz(s, exp(-V/T))
```

Unit check:

```text
(m/s) / m = 1/s
```

So `k` has units of inverse seconds.

## Step 8: Calculate Diffusivity

The final diffusivity is:

```text
D = 0.5 * k * L^2
```

where:

| Symbol | Meaning | Unit |
| --- | --- | --- |
| `D` | Diffusivity | `m^2/s` |
| `k` | Rate-like crossing factor | `1/s` |
| `L` | Hopping length from the selected cell direction | `m` |

## Full Formula

Combining the terms:

```text
D = 0.5 * L^2
    * sqrt(k_B * T / (2 * pi * m))
    * exp(-max(V) / T)
    / integral[ exp(-V(s) / T) ds ]
```

where the integral is evaluated numerically along the string path using the trapezoidal rule.

## Important Notes

- The standalone utility `cal_diffusivity` reads `T` from the input file.
- The main CUDA program also computes internal `D_1` and `D_2` values to decide whether to write the initial or final string path, but in the current source that internal comparison uses `T = 300` hard-coded.
- The standalone utility does not print labels or units, only the numeric diffusivity value.
- If `Direction` changes, the hopping length changes, so the diffusivity changes even for the same path and energy profile.
- The output potential `V` is assumed to be in Kelvin energy units, consistent with the Lennard-Jones epsilon values used in the input.

## Example

For the current files:

```bash
./cal_diffusivity example.input out.dat
```

the utility prints:

```text
3.04640e-08
```

Interpreted with units:

```text
D = 3.04640e-08 m^2/s
```
