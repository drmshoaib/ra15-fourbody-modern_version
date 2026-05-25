# Four-Body Orbital Integrator

A modern C++17 reimplementation of a four-body coplanar gravitational integrator,
originally written in FORTRAN and C++ (Visual Studio 6.0, 2003).

The core numerical algorithm is Everhart's **RA15** — a 15th-order Gauss-Radau
integrator — which is one of the most accurate methods for long-term orbital
integration. No FORTRAN compiler or legacy DLL is required.

---

## Repository structure

```
FourBodyIntegrator_Modern/
├── FourBodyIntegrator.sln          Visual Studio 2022 solution
│
├── Integrator/                     Static library: the integration engine
│   ├── include/
│   │   ├── Force.h                 Gravitational force function declaration
│   │   └── RA15.h                  RA15 integrator declaration
│   └── src/
│       ├── Force.cpp               Pairwise accelerations for 4 bodies
│       └── RA15.cpp                Everhart RA15 algorithm (C++17)
│
└── FourBody/                       Console application: simulation driver
    ├── src/
    │   └── main.cpp                Reads initial conditions, runs simulation
    └── data/
        └── initialcond.txt         Sample initial conditions (4 equal masses)
```

---

## Building

**Requirements:** Visual Studio 2022 with the "Desktop development with C++" workload.

1. Open `FourBodyIntegrator.sln`
2. Select configuration **Release | x64** (or Debug for debugging)
3. Press **Ctrl+Shift+B** to build both projects

The build order is automatic: `Integrator` (static lib) is built first,
then `FourBody` links against it.

---

## Running

The executable reads an initial conditions file and writes a trajectory file.

**Default paths** (relative to the `FourBody/` project directory):
- Input:  `data/initialcond.txt`
- Output: `result1.txt`

**Optional command-line overrides:**
```
FourBody.exe [input_file] [output_file]
```

**From the command line:**
```
cd FourBody
x64\Release\FourBody.exe data\initialcond.txt result1.txt
```

**From Visual Studio:** The working directory is set to the `FourBody/` project
folder automatically (configured in `FourBody.vcxproj`), so pressing F5 will
find `data/initialcond.txt` without any extra setup.

---

## Initial conditions file format

Plain text, one value per line, no comments:

```
x0          <- position of body 0 (x)
y0          <- position of body 0 (y)
x1          <- position of body 1 (x)
y1          <- position of body 1 (y)
x2
y2
x3
y3
vx0         <- velocity of body 0 (x)
vy0
vx1
vy1
vx2
vy2
vx3
vy3
m0          <- mass of body 0
m1
m2
m3
T           <- total simulation time
```

All values are in natural (dimensionless) units with G = 1.

---

## Output file format

Space-separated, 18 columns per line, 20 decimal places:

```
time  x0 y0 x1 y1 x2 y2 x3 y3  vx0 vy0 vx1 vy1 vx2 vy2 vx3 vy3  energy
```

The first line is a `#` header. If energy conservation is violated beyond
the tolerance (|ΔE/E₀| > 10⁻⁹), the simulation stops early and logs the
offending time step.

---

## Integration parameters

Defined in `main.cpp` and easily adjusted:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `timestep` | `0.001` | Fixed sequence size passed to RA15 |
| `maxSteps` | `1 000 000` | Maximum number of output steps |
| `energyTol` | `1e-9` | Fractional energy drift threshold |
| `ll` | `-1` | Negative → constant step size mode |
| `nclass` | `-2` | Conservative second-order: y'' = F(y, t) |

To switch to **adaptive step size**, change `ll` to a positive integer (6–12).
A value of `ll = 10` gives approximately 10⁻¹⁰ local truncation error per step.

---

## The RA15 algorithm

RA15 (Radau, 15th order) uses eight Gauss-Radau quadrature points per
sequence to build a 15th-degree polynomial approximation to the solution.
Key properties:

- **Order 15** — extremely low truncation error for smooth orbital problems
- **Adaptive step size** — automatically adjusts to maintain a target accuracy
- **Conservative systems** — uses the specialised y'' = F(y, t) code path
  (NCLASS = -2), which exploits time-reversal symmetry for better energy conservation
- **Self-starting** — no separate predictor stage required

### References

- E. Everhart, *"An efficient integrator that uses Gauss-Radau spacings"*,
  in *Dynamics of Comets: Their Origin and Evolution*, A. Carusi & G. B.
  Valsecchi (eds.), Reidel, Dordrecht, 1985, pp. 185–202.

---

## Equations of motion

Four coplanar bodies in barycentric coordinates, G = 1:

$$\ddot{\mathbf{r}}_i = \sum_{j \neq i} \frac{m_j (\mathbf{r}_j - \mathbf{r}_i)}{|\mathbf{r}_j - \mathbf{r}_i|^3}$$

Total energy (conserved quantity monitored at every step):

$$E = \underbrace{\frac{1}{2}\sum_i m_i |\dot{\mathbf{r}}_i|^2}_{K} - \underbrace{\sum_{i < j} \frac{m_i m_j}{|\mathbf{r}_i - \mathbf{r}_j|}}_{U}$$

---

## Differences from the original (2003) code

| Original | This version |
|----------|-------------|
| FORTRAN `.for` files + C++ DLL | Pure C++17, no external dependencies |
| Visual Studio 6.0 `.dsp`/`.dsw` | Visual Studio 2022 `.sln`/`.vcxproj` |
| 32-bit only | x64 target |
| Fixed-size FORTRAN arrays | `std::vector` with bounds-safe access |
| Computed `GO TO` statements | `switch` statements (identical numeric behaviour) |
| Three projects (DLL, console app, MFC GUI) | Two projects (static lib + console app) |
| Output file only | Header line + optional command-line paths |
