# SEMWavePropagation

Standalone MATLAB spectral element method (SEM) examples for isotropic elastic waves. The functions in `src` are extracted from the corresponding mesh, matrix, solver, and ParaView export ideas in the parent repository; they do not depend on its classes.

## Cases

| Dimension | Mode | Boundaries | Main field |
|---|---|---|---|
| 1D | Longitudinal | Free / PML at both ends | Axial displacement |
| 1D | Transverse | Free / PML at both ends | Scalar transverse displacement |
| 2D | SH | Free / PML on all sides | Out-of-plane displacement |
| 2D | P-SV | Free / PML on all sides | Coupled in-plane displacement |
| 3D | Full elastic | Free / PML on all sides | Three displacement components |
| 2D | P-SV half-space | Free top, PML on sides and bottom | Coupled in-plane displacement |
| 3D | Full elastic half-space | Free top, PML on sides and bottom | Three displacement components |

There are 12 runnable scripts in `examples`: five with free boundaries, five with PMLs approximating an infinite domain, and two half-space cases. For example, from the repository root:

```matlab
run('examples/wave_2d_psv_pml.m')
run('examples/wave_3d_halfspace.m')
```

Each script writes a `.pvd` collection and `.vtu` snapshots under `examples/output/`. Open the `.pvd` in ParaView. The `displacement` field is scalar for the 1D and SH modes, and vector-valued for P-SV and 3D. Each SEM element is divided into linear cells between GLL points for visualization.

Each example also animates its saved snapshots in a MATLAB figure: a line plot in 1D, a surface in 2D, or three orthogonal slices in 3D. Scalar modes show signed displacement; P-SV and 3D show displacement magnitude. The 2D surface height is scaled for visibility, while its colors represent the computed displacement. To replay a solution at a different speed, call `sem_animate_wave(model,result,struct('framePause',0.1))` after solving; `framePause` is in seconds.

Every example uses at least five GLL intervals per wavelength at the source's carrier frequency. `sem_points_per_wavelength` checks the shortest bulk wavelength (the shear wave for P-SV and 3D modes) against the largest adjacent GLL-node spacing on any axis. The half-space meshes include additional resolution for slower surface waves. A finite Hann-windowed pulse is not strictly band-limited, so this criterion applies to its carrier frequency.

## API and time step

`sem_mesh(bounds,nelem,p)` creates a Cartesian tensor-product GLL mesh. `bounds` has one `[min max]` row per dimension, `nelem` specifies elements per axis, and `p` is the polynomial order. `sem_assemble(mesh,material,mode)` builds sparse stiffness `model.K` and diagonal sparse mass `model.M` from `material.E`, `material.nu`, and `material.rho`. It also computes the largest generalized eigenvalue of `K*v = lambda*M*v` and the central-difference limit `model.dtCritical = 2/sqrt(model.lambdaMax)`.

The user must specify `options.dt`, for example `options.dt = 0.5*model.dtCritical`. Both solvers reject a value above the computed limit. They may shorten the actual step slightly to divide `tEnd` into an integer number of steps. `options.snapshotEvery` is optional.

Every source must contain `position`, `direction`, and a sampled `waveform` structure with `time` and `values` arrays. The solvers interpolate these samples linearly and apply zero force outside the supplied time interval. All examples construct a finite sine burst with a Hann window:

```matlab
waveform = sem_hann_sine(1e5,2,1); % frequency, cycles, amplitude
source = struct('position',[0.03 0.03],'direction',1,'waveform',waveform);
```

`sem_hann_sine` also accepts optional `startTime` and `samplesPerCycle` arguments. A custom sampled waveform can be assigned directly to `source.waveform`.

Use `sem_solve(model,source,tEnd,options)` for traction-free external boundaries. It integrates displacement with central differences. Use `sem_export_paraview(model,result,filenameBase)` to write a time series.

## PML and half-spaces

Build the mesh over the physical bounds **plus** the PML widths. Then call `sem_pml_profile(model,physicalBounds,thickness,targetReflection,order)` and `sem_solve_pml(model,pml,source,tEnd,options)`. `thickness` is a `dim x 2` array of lower/upper layer widths. The optional `targetReflection` defaults to `1e-3`; the optional polynomial order defaults to `2`.

The PML solver uses a split-field velocity-stress formulation. Directional first derivatives are damped within the layer, and explicit RK4 advances the fields. Its additional damping stability condition is `max(pml.sigma(:))*result.dt <= 2.5`. Setting a face's width to zero leaves it traction-free. The half-space examples use `y=0` (2D) and `z=0` (3D) as free surfaces, with depth in the positive coordinate direction.

The examples use SI units and homogeneous isotropic materials. The 1D and SH models are scalar with unit cross-section or thickness. PMLs approximate an unbounded domain, with residual reflection controlled by layer width, damping, spatial resolution, and wave frequency. The meshes are Cartesian and the code does not implement anisotropy.

## Verification

```matlab
addpath('tests')
test_sem_wave
test_sem_pml
```

The tests check mass, symmetry, rigid translations, a 1D eigenfrequency, the time-step guard, source interpolation and zero extrapolation, and PML reflection reduction. The VTU/PVD files produced by the examples are XML files readable by ParaView.

## License and citation

This project is licensed under the [Apache License 2.0](LICENSE). If you redistribute the code or a derivative work, retain the license and the attribution in [NOTICE](NOTICE) as required by the license.

If you use this software in a publication or research project, please cite Lucio de Abreu Correa and SEMWavePropagation. For questions, contact [labcorrea@gmail.com](mailto:labcorrea@gmail.com). The repository's [CITATION.cff](CITATION.cff) provides citation metadata that GitHub can display through its **Cite this repository** feature. A plain-text reference is:

> Correa, Lucio de Abreu. *SEMWavePropagation*. https://github.com/labcorrea/SEMWavePropagation
