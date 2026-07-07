# GaSSM Input and Output Guide

This guide describes the structure of `example.input` as read by the current program. The input format is fixed-order: many label lines are ignored by the code, but they should remain in place because the parser advances through the file line by line.

## Input File Overview

The input file has four main parts:

1. Cell and run settings
2. String calculation settings
3. Adsorbate molecule definition
4. Adsorbent/framework atom list

## Cell and Run Settings

```text
Nmaxa Nmaxb Nmaxc:
70 70 70
La Lb Lc dL
25.832 25.832 25.832 0.5
Alpha Beta Gamma
90.000000 90.000000 90.000000
cutoff(A) FH_signal Mass(g/mol) Tempearture(K)  Running_steps
12.90000 1 28 300.000000 10000
```

| Field | Meaning | Notes |
| --- | --- | --- |
| `Nmaxa Nmaxb Nmaxc` | Legacy grid-size fields | Currently skipped by the code. Keep both lines in the file. |
| `La Lb Lc` | Unit-cell lengths | Angstrom. |
| `dL` | Grid spacing / legacy spacing field | Read by the code, but not used later in the current source. |
| `Alpha Beta Gamma` | Unit-cell angles | Degrees. Converted internally to radians. |
| `cutoff(A)` | Interaction cutoff | Angstrom. Used to decide how many periodic images of the framework are included. |
| `FH_signal` | Flag field | Read by the code, but not used later in the current source. |
| `Mass(g/mol)` | Total adsorbate molecular mass | Should equal the sum of adsorbate site masses. |
| `Temperature(K)` | Temperature | Kelvin. Used in diffusivity-related calculations. |
| `Running_steps` | Maximum string iterations | If convergence is not reached by this count, the status flag is `0`. |

## String Calculation Settings

```text
---------String Calculation Settings---------
Direction
1
#_of_points delta_frac delta_angle_degree
401 0.000100 1.000000
convergence_setting
default
```

| Field | Meaning | Notes |
| --- | --- | --- |
| `Direction` | Overall transport direction | Options are `1`, `2`, or `3`. See below. |
| `#_of_points` | Number of string images | This is the number of points along the diffusion path. The output file has this many rows. |
| `delta_frac` | Translation update size | Fractional-coordinate step size used during string optimization. |
| `delta_angle_degree` | Rotation update size | Angular step size in degrees. Converted internally to radians. |
| `convergence_setting` | Convergence mode | In the current source, only `default` is handled. It sets both translation and rotation convergence percentages to `30`. |

## Direction Options and Impact

`Direction` defines the periodic direction that the string crosses:

| Direction | Meaning | Initial string behavior | Diffusivity length scale |
| --- | --- | --- | --- |
| `1` | Crosses along fractional `a` direction | Initial string runs from `a = 0` to `a = 1` | Uses `La` |
| `2` | Crosses along fractional `b` direction | Initial string runs from `b = 0` to `b = 1` | Uses `Lb` |
| `3` | Crosses along fractional `c` direction | Initial string runs from `c = 0` to `c = 1` | Uses `Lc` |

The direction does not mean the path is forced to remain a perfectly straight line. It sets the overall cell-crossing direction. During optimization, each string image can relax in fractional position and orientation:

```text
frac_a frac_b frac_c alpha beta gamma
```

For example, with `Direction = 1`, the path crosses the cell in the `a` direction, but the optimized chain can curve through `b` and `c` if that lowers the energy.

## Adsorbate Block

```text
------------------Adsorbate------------------
Number of sites
2
x(A) y(A) z(A) Epsilon(K) Sigma(A) Charge(e) Mass(g/mol)
0.00 0.00 0.00 92.8 3.68 0 14
0.00 0.00 1.33 92.8 3.68 0 14
```

This block defines the rigid adsorbate molecule. Each row is one adsorbate interaction site:

```text
x y z epsilon sigma charge mass
```

| Field | Meaning | Units |
| --- | --- | --- |
| `x y z` | Site coordinates relative to the adsorbate geometry | Angstrom |
| `epsilon` | Lennard-Jones epsilon | Kelvin |
| `sigma` | Lennard-Jones sigma | Angstrom |
| `charge` | Partial charge | Elementary charge |
| `mass` | Site mass | g/mol |

The code computes the adsorbate center of mass from these site coordinates and masses. The molecule is then treated as rigid, with each string image carrying a center-of-mass position and Euler angles.

## Adsorbent / Framework Block

```text
------------------Adsorbent------------------
Number of atoms
424
ID diameter(A) Epsilon(K) Charge(e) mass(g/mol) frac_x frac_y frac_z atom_name
1 2.461553158 62.3988584 1.275 65.38 0.2934 0.2066 0.2066 Zn
...
```

This block defines the framework atoms. The number after `Number of atoms` must match the number of framework atom rows that follow.

Each framework row is:

```text
ID sigma epsilon charge mass frac_x frac_y frac_z atom_name
```

| Field | Meaning | Units / Notes |
| --- | --- | --- |
| `ID` | Atom index | Read numerically, but not used as an identifier later. |
| `sigma` / `diameter(A)` | Lennard-Jones size parameter | Angstrom. The header calls it `diameter(A)`, but the code stores it as `sigma_frame`. |
| `epsilon` | Lennard-Jones epsilon | Kelvin |
| `charge` | Partial charge | Elementary charge |
| `mass` | Atomic mass | g/mol |
| `frac_x frac_y frac_z` | Framework atom fractional coordinates | Fractional unit-cell coordinates |
| `atom_name` | Element/name label | Not read by the calculation; kept for readability. |

The framework is periodically expanded based on the cutoff and cell dimensions. The status line reports the expanded atom count, not just the original framework atom count.

## Output File

The main run command is:

```bash
./GPU_string_polyatmoic example.input out.dat
```

The output file has one row per string point. With `#_of_points = 401`, the output file has 401 rows.

Each row has seven columns:

```text
frac_a frac_b frac_c alpha beta gamma external_potential
```

| Column | Meaning | Notes |
| --- | --- | --- |
| `frac_a frac_b frac_c` | Adsorbate center-of-mass position | Fractional coordinates in the simulation cell. |
| `alpha beta gamma` | Adsorbate Euler angles | Radians in the current output. |
| `external_potential` | External potential energy at that string image | Effectively in Kelvin energy units, consistent with the force-field parameters. |

The output is the path selected by the program after comparing the initial and final path diffusivity estimates. In the current code, if the initial estimate `D_1` is larger than the final estimate `D_2`, it writes the initial string; otherwise it writes the final optimized string.

## Terminal `info` Line

The program also prints a terminal line like:

```text
info: 1 800 424 0.123456
```

or:

```text
info: 0 10000 424 0.123456
```

The fields are:

| Field | Meaning |
| --- | --- |
| `1` or `0` | Convergence status. `1` means converged; `0` means maximum running steps reached without convergence. |
| Iteration count | The iteration number when the run stopped. |
| Expanded framework atoms | `Number of atoms` multiplied by the number of periodic images needed for the cutoff. |
| Runtime | CPU clock time printed by the program, in seconds. |

## Checklist for Creating a New Structure Input

When adapting the example for a new framework, update:

1. `La Lb Lc` and `Alpha Beta Gamma`
2. `cutoff(A)`, if needed
3. `Mass(g/mol)`, if changing the adsorbate
4. `Direction`, depending on which diffusion direction you want
5. Adsorbate site rows, if changing the molecule
6. `Number of atoms`
7. All framework atom rows: `ID sigma epsilon charge mass frac_x frac_y frac_z atom_name`

Keep the section headers and label lines in the same order unless the parser is changed.
