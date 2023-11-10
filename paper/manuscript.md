---
title: 'BoundaryLayerDynamics.jl v1.0: a modern codebase for atmospheric boundary-layer simulations'
correspondence: Manuel F. Schmid, m.schmid@civil.ubc.ca
author:
- name: Schmid
  firstname: Manuel F.
  institute: 1,2
  correspondence: mfs2173@columbia.edu
- name: Giometto
  firstname: Marco G.
  institute: 2
- name: Lawrence
  firstname: Gregory A.
  institute: 1
- name: Parlange
  firstname: Marc B.
  institute: 3,4
shortauthor: Schmid et al.
institute:
- name: Dept. of Civil Engineering, University of British Columbia, Vancouver, Canada
- name: Dept. of Civil Engineering & Engineering Mechanics, Columbia University, New York, NY, USA
- name: Dept. of Civil Engineering, Monash University, Melbourne, VIC, Australia
- name: University of Rhode Island, Kingston, RI, USA
abstract: |
  \noindent
  We present BoundaryLayerDynamics.jl, a new code for turbulence-resolving simulations of atmospheric boundary-layer flows as well as canonical turbulent flows in channel geometries.
  The code performs direct numerical simulation as well as large-eddy simulation using a hybrid (pseudo)spectral and finite-difference approach with explicit time advancement.
  Written in Julia, the code strives to be flexible and adaptable without sacrificing performance, and extensive automated tests aim to ensure that the implementation is and remains correct.
  We show that the simulation results are in agreement with published results and that the performance is on par with an existing Fortran implementation of the same methods.
short-summary: |
  Turbulence-resolving flow models have strict performance requirements, as simulations often run for weeks using hundreds of processes.
  Many flow scenarios also require the flexibility to modify physical and numerical models for problem-specific requirements.
  With a new code written in Julia we hope to make such adaptations easier without compromising on performance.
  In this paper we discuss the modeling approach and present validation and performance results.
acknowledgements: |
  We acknowledge the financial support from the NSERC (Canada), the Swiss National Science Foundation, the Department of Civil Engineering and Engineering Mechanics at Columbia University, Monash University, and the University of Rhode Island.
  This work made use of the Stampede2 system at the Texas Advanced Computing Center (TACC) through the Extreme Science and Engineering Discovery Environment \citep[XSEDE;][]{Towns+2014}, which is supported by National Science Foundation grant number ACI-1548562.
availability: |
  BoundaryLayerDynamics.jl is open source software and available under the MIT License through the official Julia package repository and on GitHub\footnote{\url{https://github.com/efpl-columbia/BoundaryLayerDynamics.jl}}, where the public repository of the package is currently hosted.
  The version described in this article is archived on Zenodo [@Schmid2023a].
  The data and code required to reproduce this paper are also made available on Zenodo [@Schmid+2023a].
keywords:
- large-eddy simulation
- direct numerical simulation
- atmospheric boundary layer
key-points: # Each must be 140 characters or fewer with no special characters or punctuation and must be complete sentences
- BoundaryLayerDynamics.jl is a new code for turbulence-resolving simulations of atmospheric boundary-layer flows.
- The code is validated for direct numerical simulation and large-eddy simulation of turbulent channel flows.
- Performance is shown to be in line with a Fortran implementation of the same numerical approach.
plain-language-summary: |
  When the atmosphere interacts with the built and natural environment at the surface of our planet, the air moves in complex, swirling patterns.
  Simulating these motions in detail requires the equivalent of hundreds or even thousands of typical personal computers.
  As computational models have to make efficient use of these resources, they are complex pieces of software and it can be challenging to make significant changes without jeopardizing the efficiency or correctness of the simulations.
  However, changes are often necessary to represent a wide range of flow scenarios in an optimal manner.
  This paper presents a new computational model that aims to make such changes easier and more robust.
  The paper shows that the model produces the expected results for typical configurations and that its performance is in line with a more classical implementation of the same modeling approach.
  While the adaptability of the new model can only be truly evaluated over time, the code is already being used to explore questions that would have been difficult to answer otherwise.
header-includes:
- '\newcommand{\kx}{\kappa_1}'
- '\newcommand{\ky}{\kappa_2}'
- '\newcommand{\iu}{{i\mkern1mu}}'
- '\newcommand{\euler}{\mathrm{e}}'
- '\newcommand{\ndkern}{\mkern -.2mu}'
- '\newcommand{\re}{R \ndkern e}'
- '\newcommand{\filt}[1]{\widetilde{#1}}'
- '\newcommand{\vfd}[1]{\boldsymbol{\hat{#1}}}'
- '\newcommand{\Dfd }[0]{ \hat{D} }'
- '\newcommand{\Op}[1]{\operatorname{#1}}'
- '\newcommand{\vphi}{\vfd{\varphi}}'
- '\newcommand{\Cpp}{C\hspace{-.02ex}\raisebox{.3ex}{\textbf{\scriptsize ++}}}'
- '\newcommand{\dd}{\mathrm{d}}'
- '\newcommand{\dv}[2]{\frac{\dd #1}{\dd #2}}'
- '\newcommand{\dvinline}[2]{\dd #1 / \dd #2}'
- '\newcommand{\pdvinline}[2]{\partial #1 / \partial #2}'
- '\newcommand{\pdvn}[1]{\frac{\partial}{\partial #1}}'
- '\newcommand{\pdv}[2]{\frac{\partial #1}{\partial #2}}'
- '\newcommand{\pdvd}[3]{\frac{\partial^2 #1}{\partial #2 \partial #3}}'
- '\newcommand{\qq}[1]{\quad \text{#1} \quad}'
- '\newcommand{\qqr}[1]{\text{#1} \quad}'
- '\newcommand{\vb}[1]{\boldsymbol{#1}}'
- '\newcommand{\abs}[1]{\lvert #1 \rvert}'
- '\newcommand{\eval}[1]{\left. #1 \right|}'
- '\newcommand{\order}[1]{\mathcal{O} \mkern -1mu \left( #1 \right)}'
include-before:
draft: true
conclusions: |
  Turbulence-resolving ABL flow simulations are subject to a number of competing requirements that have to be considered when developing simulation code. Availability of physical and numerical models, performance and scalability, ease of use and ease of modification, safeguards against implementation and usage errors, as well as license terms may vary considerably and trade-offs are often inevitable.
  The Julia programming language is promising more favorable trade-offs by offering the ergonomics of a modern high-level language without sacrificing performance.

  In this paper, we have introduced a new code for turbulence-resolving flow simulations, designed for the requirements of atmospheric boundary-layer research and written in Julia.
  The performance is shown to be in line with a Fortran implementation of the same modelling approach.
  In fact, it even appears that easier experimentation with algorithmic approaches and implementation trade-offs might have a stronger impact on performance than the remaining computational overhead compared to highly optimized Fortran compilers.

  The code also places a focus on automated testing and minimizing the chances for errors both during development and usage.
  This is particularly important in exploratory research, where the expected behavior of a new model or flow system is not known *a priori* and it is difficult to discern between inconspicuous errors and novel results.

  The code provides the core functionality for both direct numerical simulation and large-eddy simulation in channel-flow geometries.
  In the future, we expect to expand the scope by adding functionality such as more advanced subgrid-scale models, support for temperature, humidity and transport of passive scalars, and partially resolved complex terrain.
contributions: |
  MFS developed the model code, performed validation and performance testing, and prepared the manuscript.
  MGG, GAL, and MBP acquired project funding, supervised the work, and reviewed the manuscript.
  Access to computational resources was provided by MGG.
conflicts: |
  The authors declare that they have no conflict of interest.
---

some historical perspective: ABL research has relied on computation for a long time
- numerical simulation of ABL flows
- deardorff 1970–1972 early les
––
Since Deardorff’s early studies [@Deardorff1969; @Deardorff1970a; @Deardorff1970b], numerical simulations of the three-dimensional, unsteady flow field have become an integral part of microscale atmospheric boundary-layer (ABL) research.
Direct numerical simulation (DNS) provides an extremely accurate tool to study fundamental properties of turbulent flows and their scaling from low to moderate Reynolds numbers [@MoinMahesh1998].
Large-eddy simulation (LES) provides a numerical model for a wide range of ABL flow phenomena at realistic Reynolds numbers while relying on modest, well-supported modeling assumptions [@MeneveauKatz2000; @Stoll+2020].
Together, DNS and LES constitute the backbone of the computational study of turbulent flow dynamics and have contributed many insights to our current understanding.

available codes
- many codes implement simulation methods
- large community projects, open-source codes: openfoam, palm, wrf
- many research groups have their own projects that are passed between researchers and evolve over time
- also some commercial options such as ansys fluent
- Q: include some detail about Parlange–Albertson code history here?
- choice of code might be guided to some degree by familiarity and access, but also real differences along many dimensions
––
Many different implementation of these methods are in use for ABL research and continue to be actively developed.
Projects such as PALM [@Maronga+2020], OpenFOAM [@Chen+2014], and WRF [@Skamarock+2021] develop open-source codes with broad applicability in a community effort.
In addition, many research groups have their own codes that have been developed and extended over decades and are passed person-to-person.
Some studies also rely on commercial software such as Ansys Fluent, although not having access to implementation details is problematic for scientific reproducibility.

dimensions of code differences
- perhaps most important: physical and mathematical methods
- performance also very important: many studies limited by performance
- ease of use: might require a lot of knowledge to set up & run correctly, and may or may not be documented, may rely to different degrees on knowledge & experience passed down person to person
- many turbulence studies also require changes to the code, so also important how easy it is to change code without affecting performance & correctness)
- all dimensions desirable, but sometimes require large effort & sometimes are in conflict with each other & require trade-offs
––
Codes differ along many dimensions, perhaps most importantly in the physical models that are implemented.
The numerical methods used to compute a solution can also lead to important differences, particularly for LES, where the smallest resolved scales of motion are an integral part of the turbulence dynamics [@KravchenkoMoin1997].
The performance characteristics of a code are also central as most turbulence-resolving simulations continue to be limited by their computational cost.
Furthermore, there can be important differences in the effort and model-specific experience required for setting up simulations and making changes and additions to the source code while ensuring the correctness of the results.
When developing a new code or selecting an existing model for a simulation, these different qualities have to be weighed against each other and trade-offs are inevitable.

new trade-offs enabled by Julia
- trade-offs not constant but shift over time, e.g. with availability of new hardware, software, and mathematical & physical methods
- one such development is Julia (publicly launched in 2012, stable with 1.0 in aug 2018, further matured since)
- promises no or minimal performance overhead compared to statically compiled languages like Fortran, C, C++ commonly used for implementing performance-critical numerics
- promises easier development (automatic memory management, automatic typing, easier to keep orthogonal functionality separate with multiple dispatch, modern package system)
––
The advent of the Julia programming language [@Bezanson+2017] represents a shift in the landscape of possible trade-offs between conflicting goals.
Publicly launched in 2012 and stabilized with version 1.0 in 2018, Julia promises to combine the performance of Fortran, C, and \Cpp with the convenience of Python, Matlab, and R.
Automatic memory management, dynamic typing with type inference, multiple dispatch, and a built-in modern package manager facilitate rapid development of clear, concise code that keeps orthogonal functionality separate.
At the same time, Julia’s “just-ahead-of-time” compilation model allows code to run with no or minimal computational overhead.

goals of this project
- new code for research of turbulent flows, with focus on ABL research
- goals: basic functionality for dns & les of turbulent channel flow, high performance (similar to existing), high expandability, high correctness both now & in future
- methods: Julia, automated tests, orthogonal functionality
––
In this paper, we present and discuss BoundaryLayerDynamics.jl [@Schmid2023a], a new code for turbulence-resolving flow simulation optimized for ABL research.
The code has been written to provide core functionality for DNS and LES of channel-flow configurations with a focus on making it easy to use, adapt, and extend the code without jeopardizing the correctness of the results.
To achieve better trade-offs along these dimensions, the implementation relies on the Julia programming language, on automated testing, and on a modular design.

scope of the project
- none of the functionality is limited to ABL flows
- generally applies to turbulence-resolving simulations of flow in channel-like geometries
- as design is guided by needs of ABL flows, some design decisions may not be appropriate for all turbulent flows
––
The suitability of the new code is of course not limited to simulations of ABL turbulence.
In fact, none of the current functionality is specific to ABL applications.
However, the choice of physical models and numerical methods is guided by the needs of such applications and future developments will similarly prioritize those use cases.

overview over paper
- concise numerics: (pseudo)spectral horizontally, finite differences vertically, with MPI parallelism, explicit time integration, dns & les with simple wall & sgs model
- describe automated validation, provide results for manual validation of turbulent flows
- performance: show that it matches Fortran as used in Vancouver paper & explore performance characteristics
––
The code solves the incompressible Navier–Stokes equations relying on (pseudo\discretionary{-)}{}{)}spectral and finite difference methods for discretization in the horizontal and vertical directions respectively, making use of the Message Passing Interface (MPI) for parallelization in the vertical direction, and performing fully explicit time integration.
The details of the numerical methods are described in section \ref{sec:numerics}.
The implementation is validated via a number of automated tests as described in section \ref{sec:validation}, where we also present a validation against DNS and LES results computed with different codes.
The performance analysis presented in section \ref{sec:performance} shows that the computational cost is comparable to a Fortran implementation of the same numerical approach and that parallel performance scales favorably up to the maximum supported number of parallel processes.

BoundaryLayerDynamics.jl is open source software and available under the MIT License through the official Julia package repository and on GitHub\footnote{\url{https://github.com/efpl-columbia/BoundaryLayerDynamics.jl}}, where the public repository of the package is currently hosted.
The version described in this article is archived on Zenodo [@Schmid2023a].

# Governing equations and numerical methods {#sec:numerics}

## what is the overarching goal? {.outline}

accurate simulation of turbulent flow dynamics for research
- primary goal: simulate atmosphere at turbulence-resolving resolutions at different length scales with focus on ABL
- secondary goal: simulate other turbulent flows in engineering & in natural environment
- trade-offs: numerics optimized for best results rather than widest applicability, initial functionality limited compared to more mature code
––
The choice of governing equations and numerical methods is guided by the goal of studying the turbulent flow dynamics in the atmospheric boundary layer.
Other turbulent flows in engineering and in the natural environment are also considered insofar as their requirements do not conflict with those of atmospheric boundary-layer flows.

## what is the physical problem that is being solved? {.outline}

starting point: incompressible NS are very accurate (not much simplification)
- rotational form of incompressible NS in Cartesian coordinates
- constant density & viscosity
- need to define all introduced variables
- not clear: why rotational form? different forms should be equivalent up to machine precision (is this correct? should re-read CTR paper), but advection form would require 9 iFFTs (ui, dui/dx1, dui/dx2) and 3 FFTs to compute in physical domain, divergence form would require 3 iFFTs (ui) and 6 FFTs (uiuj), while rotational form requires 6 iFFTs (ui, wi) and 3 FFTs
––
It is well-established that the Navier–Stokes equations are an extremely accurate physical model for flows with a Knudsen number $K \ndkern n \ll 1$ and that compressibility effects are minimal for flows with a Mach number $M \ndkern a \lesssim 1/3$ [@Panton2013].
Since these conditions are met for atmospheric boundary-layer flows, the incompressible Navier–Stokes equations,
\begin{equation}
  \pdv{u_i}{t} + u_j \left( \pdv{u_i}{x_j} - \pdv{u_j}{x_i} \right) =
  \frac{1}{\re} \pdvd{u_i}{x_j}{x_j}
  - \pdv{p}{x_i}
  + f_i
  \qq{and}
  \pdv{u_i}{x_i} = 0
  \,,
\label{eq:ns-pd}
\end{equation}
are used as the mathematical model for this work, given here with the rotational form of the advection term [@Orszag1971a], for which the total kinematic pressure $p = p^\mathrm{\mkern2mu gauge}/\rho + \frac{1}{2} u_i u_i$ includes both the static and dynamic pressure.
Quantities are non-dimensionalized with a length scale $\mathcal{L}$ and a velocity $\mathcal{U}$, producing the Reynolds number $\re = \mathcal{U} \mathcal{L}/\nu$.
The Cartesian coordinates $x_i$ denote the primary horizontal ($i=1$, usually streamwise), secondary horizontal ($i=2$, usually cross-stream), and vertical ($i=3$) directions while the corresponding velocity components are given as $u_i$.
The term $f_i$ denotes the components of a body force, usually gravity or a constant pressure gradient, and the kinematic viscosity $\nu$ and fluid density $\rho$ are assumed scalar constants.

other physical processes could be added, but are currently not
- Coriolis, temperature and humidity are the dominant ones
- could be added relatively easily
- for neutral atmosphere, relatively accurate
––
While this formulation can serve as a reasonable representation of a neutrally stratified atmospheric boundary layer, there are many important ABL processes that are not included.
Coriolis forces, temperature, and humidity in particular are of central importance for many applications.
The current code is meant to provide the core functionality necessary for ABL flow simulations and serve as a foundation for a more comprehensive set of physical and numerical models that can be added over time.

high Re: use filtered NS instead
- cannot resolve all scale → filter equations
- define filter with integral
- sgs stresses: have to be modeled as function of filtered velocity field
––
For flows with a moderate to high Reynolds number, current computational capabilities generally do not permit resolving the full range of scale of motions.
In this case, the filtered Navier–Stokes equations,
\begin{equation}
  \pdv{\filt{u}_i}{t} + \filt{u}_j \left( \pdv{\filt{u}_i}{x_j} - \pdv{\filt{u}_j}{x_i} \right) +
  \pdv{\tau_{ij}^\mathrm{sgs}}{x_j} =
  \frac{1}{\re} \pdvd{\filt{u}_i}{x_j}{x_j}
  - \pdv{\filt{p}}{x_i}
  + \filt{f}_i
  \qq{and}
  \pdv{\filt{u}_i}{x_i} = 0
  \,,
\label{eq:ns-les}
\end{equation}
are used as the computational model, where $\filt{u}_i$ represents the spatially filtered velocity field, i.e. $\filt{u}_i = \int G(\vb{r},\vb{x}) u_i(\vb{x}-\vb{r}) \dd{\vb{r}}$ with $G$ defining the filtering operation.
The subgrid-scale stress tensor $\tau_{ij}^\mathrm{sgs} = \tau_{ij}^\mathrm{R} - \frac{1}{3} \tau_{ii}^\mathrm{R} \delta_{ij}$ represents the anisotropic component of the residual stress tensor $\tau_{ij}^\mathrm{R} = \filt{u_i u_j} - \filt{u}_i \filt{u}_j$ and has to be modeled as a function of the resolved velocity $\filt{u}_i$.
The modified pressure $\filt{p} = \filt{p}^\mathrm{\mkern2mu gauge}/\rho + \frac{1}{2} \filt{u}_i \filt{u}_i + \frac{1}{3} \tau_{ii}^\mathrm{R}$ now includes contributions from the filtered gauge pressure $\filt{p}^\mathrm{\mkern2mu gauge}$, the resolved kinetic energy, and the unresolved kinetic energy.
The forcing term $\filt{f}_i$ is simply the spatially filtered $f_i$.
In the following, the same notation is used to represent both the unfiltered (DNS) and filtered (LES) equations to simplify the notation.

The past decades have seen several efforts to develop a suitable model for $\tau_{ij}^\mathrm{sgs}$ [@Smagorinsky1963; @Schumann1975; @BardinaFerzigerReynolds1980; @Germano+1991; @MeneveauLundCabot1996; @PorteAgelMeneveauParlange2000; @BouZeidMeneveauParlange2005].
The current implementation includes the static @Smagorinsky1963 subgrid-scale model
\begin{equation}
  \tau_{ij}^\mathrm{sgs} =
  - 2 l_\mathrm{S}^2 \mathcal{S} S_{ij}
  \,,
\end{equation}
where $S_{ij} = 1/2 \left( \pdv{u_i}{x_j} + \pdv{u_j}{x_i} \right)$ is the resolved strain rate, $\mathcal{S} = \sqrt{2 S_{ij} S_{ij}}$ is the characteristic or total strain rate, and $l_\mathrm{S}$ is the Smagorinsky lengthscale, taken to be the product of the filter width $\Delta$ and a constant Smagorinsky coefficient $C_\mathrm{S}$.

horizontal domain size & boundary conditions
- specifying boundary conditions remains difficult
- common strategy is use periodic boundary conditions
- disadvantage: problem has to be formulated in periodic manner, domain needs to be large enough to decorrelate flow, which can be difficult especially if flow inside domain is not homogeneous
- QUESTION: should we always have repeating elements inside domain?
- NOTE: in a way, we have the same considerations at both ends of the wavenumber range – just like we want to resolve all the small wavenumbers that are relevant, we also want the domain to be large enough that the lowest wavenumber is not very relevant – but it is possible that the lowest wavenumber is due to the repeating nature of the flow/geometry and that repeating it twice would give the same behavior at the second wavenumber
- periodic: horizontal only
- advantage: physical flow everywhere, can express flow field with periodic basis functions
- (put equation here)
––
Specifying boundary conditions for turbulent flows remains a challenging research problem, since chaotic velocity fluctuations have to be prescribed in a physically accurate manner.
To avoid these difficulties, turbulence-resolving simulations are often run with periodic boundary conditions.
While this requires that the problem is formulated such that it can be approximated with a periodic flow field, unphysical border regions are avoided and accurate results can be obtained in the whole domain as long as the domain is large enough to accommodate all relevant scales of motion.
Furthermore, the flow field can then be expressed in terms of periodic basis functions.
For a horizontally-periodic domain of size $L_1 \times L_2$, the velocity field can be written as
\begin{equation}
  u_i(x_1, x_2, x_3) = \mathop{\sum\sum}_{\substack{-\infty < \kx < \infty \\ -\infty < \ky < \infty}} \hat{u}_i^{\kx\ky}(x_3) \, \euler^{\iu \kx 2 \pi x_1/L_1} \, \euler^{\iu \ky 2 \pi x_2/L_2}
  \,.
\end{equation}
With similar expressions for the pressure $p$ and the forcing $f_i$, we can rewrite the governing equations for a single mode with wavenumber $\kx$ in $x_1$-direction and $\ky$ in $x_2$-direction,
\begin{equation}
  \begin{aligned}
&\pdvn{t} \hat{u}_i^{\kx\ky}(x_3)
+ \mathop{\sum\sum}_{\substack{-\infty < \kx^\prime < \infty \\ -\infty < \ky^\prime < \infty}}
\left( \Dfd_j^{\kx^\prime \ky^\prime} \hat{u}_i^{\kx^\prime \ky^\prime}(x_3) - \Dfd_i^{\kx^\prime \ky^\prime} \hat{u}_j^{\kx^\prime \ky^\prime}(x_3) \right)
\hat{u}_j^{(\kx-\kx^\prime)(\ky-\ky^\prime)}(x_3)
\\
&\hspace{1em}
+ \Dfd_j^{\kx\ky} \hat{\tau}_{ij}^{\mathrm{sgs}\,\kx\ky}(x_3)
= \frac{1}{\re} \left( \Dfd_j^{\kx\ky} \right)^2 \hat{u}_i^{\kx\ky}(x_3)
- \Dfd_i^{\kx\ky} \hat{p}^{\kx\ky}(x_3)
+ \hat{f}_i^{\kx\ky}(x_3)
\,,
  \end{aligned}
  \label{eq:mom-fd}
\end{equation}
where $\Dfd_1^{\kx\ky} = \frac{2\pi\iu\kx}{L_1}$, $\Dfd_2^{\kx\ky} = \frac{2\pi\iu\ky}{L_2}$, and $\Dfd_3^{\kx\ky}=\pdvn{x_3}$ are the differential operators, with $\iu$ denoting the imaginary unit.
For direct numerical simulations the subgrid-scale term is omitted.
The continuity equation becomes
\begin{equation}
\Dfd_i^{\kx\ky} \hat{u}_i^{\kx\ky}(x_3) = 0
\,.
\label{eq:cont-fd}
\end{equation}

boundaries in vertical direction
- periodic boundary conditions in vertical less useful if boundary-layer flows are of interest
- for channel geometry, bcs are needed at 0 & l3, with integral constraint
- for engineering flows such as smooth-wall closed & open channel flow, dirichlet & neumann BCs are sufficient
- complex surfaces found in ABL require significant modeling effort, either resolving (simplified) features or representing effect of surface with a wall model
- wall-models are generally formulated for discretized equations
––
In the vertical direction, turbulence is not homogeneous for ABL flows and periodic boundary conditions are not applicable.
For a channel geometry, boundary conditions have to be specified for $u_i(x_3 = 0)$ and $u_i(x_3 = L_3)$ with the constraint that
\begin{equation}
  \int_0^{L_2} \int_0^{L_1} u_3(x_3=0) - u_3(x_3=L_3) \, \dd{x_1} \dd{x_2} = 0\,,
\end{equation}
which can be obtained from integrating the continuity equation over the whole domain.
A number of engineering flows such as smooth-wall open and closed channel flows can be modeled with Dirichlet and Neumann boundary conditions.
The complex boundaries of atmospheric flows can require significant modeling effort and simulations generally have to partially resolve surfaces (e.g. immersed boundary method, terrain-following coordinates) or represent their effect with a wall model for $\tau_{i3}^\mathrm{sgs}$, usually formulated for the discretized equations [@PiomelliBalaras2002].
The current implementation includes an algebraic equilibrium rough-wall model defined similar to @MasonCallen1986 with
\begin{equation}
  \tau_{i3}^\mathrm{sgs}(x_3=0) =
  \frac{-\kappa^2 \sqrt{u_1(x_3^\mathrm{ref})^2 + u_2(x_3^\mathrm{ref})^2}}
  {\log(x_3^\mathrm{ref} / z_0)^2} u_i(x_3^\mathrm{ref})
\end{equation}
for $i=1,2$ and $\tau_{33}^\mathrm{sgs}(x_3=0) = 0$, where $z_0$ is the roughness length, $\kappa \approx 0.4$ is the von Kármán constant, and $x_3^\mathrm{ref} > z_0$ is a reference height at which the (resolved) velocity is obtained, usually chosen as the first grid point.
To improve the near-wall behavior of the subgrid-scale model, the Smagorinsky length scale is adjusted to $l_\mathrm{S}^{-n} = \left(C_\mathrm{S} \Delta\right)^{-n} + \left(\kappa x_3 \right)^{-n}$ as proposed by @MasonThomson1992, with $n=2$ as the default value.

## how is it approximated in space? {.outline}

spatial discretization: wavenumbers and grid points
- limited number of wavenumbers and grid points are solved for
- grid stretching can be applied in vertical direction
- grid is staggered for vertical velocity
––
For the numerical solution of equations \eqref{eq:mom-fd} and \eqref{eq:cont-fd} we limit ourselves to $N_1 \times N_2$ wavenumbers at $N_3$ vertical grid points.
The wavenumbers are selected symmetrically around $\kappa_i = 0$, i.e. $\abs{\kx} \le (N_1-1)/2$ and $\abs{\ky} \le (N_2-1)/2$.
This results in an odd number of wavenumbers in each direction and avoids the need for a special treatment of Nyquist frequencies.
Since $\hat{\phi}^{-\kx-\ky} = \hat{\phi}^{\kx\ky\ast}$ for any real-valued $\phi$, we only need to explicitly solve for half the modes and can obtain the others through complex conjugation.
In vertical direction, equidistant grid points are selected from the interval $[0, 1]$, which is then mapped to the domain with a function $x_3: \left[0, 1\right] \to \left[0, L_3\right]$, $\zeta \mapsto x_3(\zeta)$.
This function can be used for grid stretching in the vertical direction; the choice of $x_3: \zeta \mapsto L_3 \zeta$ defines a uniform grid.
A staggered arrangement of grid points with $\zeta_C$ at the center of the $N_3$ segments and $\zeta_I$ at the $N_3-1$ interfaces between them, i.e.,
\begin{equation}
  \begin{aligned}
  \zeta_C &\in \left\{ \frac{1/2}{N_3}, \frac{3/2}{N_3}, \dots, \frac{N_3-1/2}{N_3} \right\}
  && \qq{for} u_1, u_2, p, f_1, f_2, \tau_{ii}^\mathrm{sgs}, \tau_{12}^\mathrm{sgs},
  \\
  \zeta_I &\in \left\{ \frac{1}{N_3}, \frac{2}{N_3}, \dots, \frac{N_3-1}{N_3} \right\}
  && \qq{for} u_3, f_3, \tau_{13}^\mathrm{sgs}, \tau_{23}^\mathrm{sgs},
  \end{aligned}
\end{equation}
avoids the need to specify boundary conditions for the pressure field, prevents odd–even decoupling, and results in a smaller effective grid spacing [@FerzigerPericStreet2020].
When running large-eddy simulations, the discretization implicitly defines the spatial filter $G$ and the filter width $\Delta$ is taken to be $\Delta = \sqrt[3]{\Delta_1\Delta_2\Delta_3}$ [@ScottiMeneveauLilly1993] with $\Delta_1 = L_1 / N_1$, $\Delta_2 = L_2 / N_2$, and $\Delta_3 = 1 / N_3 \dd x_3 / \dd \zeta$.

spatial discretization: derivatives
- horizontal derivatives are exact
- vertical derivatives are approximated with second-order finite differences
- truncation error for derivative
- mention boundary stencils (give full stencils? in appendix?)
––
The horizontal derivatives $\hat{D}_1^{\kx\ky}$ and $\hat{D}_2^{\kx\ky}$ can be computed exactly.
For the vertical derivative $\hat{D}_3^{\kx\ky}$, we use central second-order finite differences on the staggered $\zeta$-nodes [@MoinVerzicco2016] as well as the analytical derivative of $x_3(\zeta)$, i.e.
\begin{equation}
  \eval{\hat{D}_3^{\kx\ky} \hat{\phi}^{\kx\ky}}_{\zeta}
  = \eval{\dv{\zeta}{x_3}}_\zeta \eval{\pdv{\hat{\phi}^{\kx\ky}}{\zeta}}_\zeta
  = \eval{\dv{\zeta}{x_3}}_\zeta \frac{\hat{\phi}^{\kx\ky}(\zeta + \delta \zeta / 2) - \hat{\phi}^{\kx\ky}(\zeta - \delta \zeta / 2)}{\delta \zeta} + \order{\delta \zeta^2}
\end{equation}
for any field $\phi$, where $\delta\zeta \equiv 1/N_3$ is the grid spacing in the $\zeta$ coordinate.
Vertical derivatives are therefore evaluated at the opposite set of grid points to the ones where $\phi$ is defined, as typical for staggered grids.
At the boundary, one-sided second-order stencils are employed.
This approximation of the vertical derivatives results in a truncation error of order $\order{\delta \zeta^2}$.

spatial discretization: error from non-linear advection term
- interpolations in vertical direction → truncation error
- sum over limited range of wavenumbers → truncation error
––
The non-linear advection term of Eq. \eqref{eq:mom-fd} requires further approximations.
First, some of the terms (e.g. $i=1, j=3$) are evaluated at the opposite set of vertical grid points than where they are required and have to be interpolated.
With the simple interpolation $\hat{\phi}\left(\zeta\right) = 1/2 \bigl( \hat{\phi}\left(\zeta-\delta\zeta/2\right) + \hat{\phi}\left(\zeta+\delta\zeta/2\right) \bigr) + \order{\delta \zeta^2}$, the truncation error generally increases but remains of order $\order{\delta\zeta^2}$.
Furthermore, the double sum can only be computed over the resolved range of wavenumbers, i.e. $\abs{\kx^\prime} \le (N_1-1)/2$ and $\abs{\ky^\prime} \le (N_2-1)/2$, producing another truncation error that decreases exponentially with the number of resolved wavenumbers.
The same applies to the non-linear expressions involved in the evaluation of $\tau_{ij}^\mathrm{sgs}$.
This discretization of the advection term in rotational form conserves kinetic energy in the absence of time-integration errors [@Mansour+1979] as long as the grid is uniform.

spatial discretization: computation of non-linear advection term
- FFTs to avoid expensive convolution
- number of grid points in physical domain is a parameter
- for npd > 3/2 n, product of two terms has no aliasing error (only truncation error)
- order of vertical derivative and FFT
- sgs bwd: 3 velocity + 6 horizontal derivatives
- sgs fwd: 6 symmetric tij (or 5 if we don’t need t33)
- 3/2 rule: Npd = 3*kmax+1
- odd N: kmax=(N-1)/2, Npd = 3/2 (N-1) + 1 = 3/2 N - 1/2
- even N: kmax=(N-2)/2, Npd = 3/2 (N-2) + 1 = 3/2 N - 2
––
To simplify the computation of non-linear terms and avoid evaluating expensive convolutions, those terms are computed on $N_1^\mathrm{PD} \times N_2^\mathrm{PD}$ equidistant grid points in the physical domain, relying on the fast Fourier transform (FFT) algorithm for forward and backward transforms [@Orszag1969; @Orszag1971b].
In principle, $N_1^\mathrm{PD}$ and $N_2^\mathrm{PD}$ are parameters that can be chosen independently of $N_1$ and $N_2$, but the choice of $N_i^\mathrm{PD} \ge 1 + 3 \kappa_i^\mathrm{max}$  avoids introducing aliasing errors for a simple product of two variables such as the resolved advection term [@PattersonOrszag1971].
In this case, the physical-domain evaluation is equivalent to a true spectral Galerkin method computing the convolution of Eq. \eqref{eq:mom-fd} over all resolved wavenumbers.
Contributions from wavenumbers $\abs{\kappa_i} > (N_i-1)/2$ are discarded upon return to the Fourier domain and vertical derivatives can be computed before or after the horizontal Fourier transforms as the two operations commute.

For the more complex non-linear expressions introduced by the SGS model, full dealiasing is generally not feasible and physical-domain evaluations incur aliasing errors in addition to the truncation errors.
This approach, dubbed the pseudospectral method by @Orszag1971b, can achieve similar accuracy to a Galerkin method [@Orszag1972].
While it is common to set $N_i^\mathrm{PD} = N_i$ for pseudospectral approximations and only discard the Nyquist wavenumber, $N_i^\mathrm{PD}$ can in principle be chosen freely for more control over truncation errors.
Furthermore, it can be beneficial to choose different values of $N_i^\mathrm{PD}$ for each non-linear term, since computing the resolved advection requires only nine transforms and is known to be sensitive to aliasing errors [@KravchenkoMoin1997; @Margairaz+2018] while the evaluation of the SGS model requires 15 transforms.

summary of spatial discretization: semi-discrete equations
- collect velocity components into an array, pressure as well
- explain components of semi-discrete equations
- forcing could be dependent on space and on velocity
––
Combining the velocity components into a single vector $\vfd{u}$ of length $N_1 \times N_2 \times (3 N_3 - 1)$, the spatially discretized momentum equation can be written as
\begin{equation}
  \dv{\vfd{u}}{t} = \Op{Adv}(\vfd{u}) + \frac{1}{\re} \Delta \vfd{u} + \frac{1}{\re} \vfd{b}_\Delta - \Op{G} \vfd{p} + \vfd{f}
  \,.
  \label{eq:mom-vec}
\end{equation}
Here, $\Op{Adv}(\cdot)$ is the non-linear advection operator that includes both resolved and subgrid-scale contributions.
$\Delta$ is the linear operator that computes the Laplacian of each velocity component, with $\vfd{b}_\Delta$ denoting the contributions from vertical boundary conditions.
$\Op{G}$ is the linear operator that computes the gradient of the pressure, discretized as a vector $\vfd{p}$ of size $N_1 \times N_2 \times N_3$.
The forcing term $\vfd{f}$, a vector of the same length as $\vfd{u}$, can be constant in time and space (e.g. pressure-driven channel flow), vary in time only (e.g. constant-mass-flux channel flow), vary in space only (e.g. baroclinic flow), or even vary in time and space as a function of $\vfd{u}$, which can be used to model vegetation drag or to represent complex geometry with an immersed-boundary method.
The discrete continuity equation,
\begin{equation}
  \Op{D} \vfd{u} + \vfd{b}_\mathrm{D} = 0
  \,,
  \label{eq:cont-vec}
\end{equation}
contains the linear divergence operator $\Op{D}$ with the contributions $\vfd{b}_\mathrm{D}$ from the vertical boundary conditions of $u_3$.

summary of space discretization: advantages & disadvantages of approach
- (justify with track record of Parlange–Albertson code & papers Beatrice & Weiyi, but also mention drawbacks for complex geometry and mass conservation)
––
This hybrid approach, relying on spectral approximations in horizontal direction (pseudospectral for the evaluation of $\tau_{ij}^\mathrm{sgs}$) and second-order-accurate finite differences in vertical direction, has long been employed for computational studies of turbulent flows in channel geometries [@MoinKim1982; @Moeng1984; @AlbertsonParlange1999a; @AlbertsonParlange1999b].
It combines the fast convergence and low dissipation of spectral methods [@GiacominiGiometto2021] with the ease of parallelization and simple handling of boundary conditions of finite differences.
Conversely, handling complex domains and non-periodic boundaries can be problematic, though still possible [@ChesterMeneveauParlange2007; @Schmid2015; @LiBouZeidAnderson2016].

## how is it approximated in time? {.outline}

approach for time integration
- follow perot 1993
- corresponds to fractional step method, but rigorous derivation from fully discrete equations
- no boundary conditions for pressure required (pressure equation depends on boundary conditions of continuity equation)
- easy to extend without introducing errors (unsteady BCs, constant-flux forcing)
––
Following @Perot1993, we obtain the expressions for time integration through a block LU decomposition of the fully-discretized equations.
This results in expressions in the style of the fractional step method [@Chorin1968; @Temam1969], but avoids the need for boundary conditions for the intermediate velocity and the pressure on a staggered grid and can easily be adapted when new terms are included or different numerical methods are employed.

fully-discretized system for ab methods
––
Adams–Bashforth methods solving ordinary differential equations of the form $\dvinline{\vb{u}}{t} = f(\vb{u}),\; \vb{u}(t_0) = \vb{u}_0$ can be written as $\vb{u}^{(n+1)} = \vb{u}^{(n)} + \Delta t \sum_{i=0}^{s-1} \beta_i f(\vb{u}^{(n-i)})$, where $\beta_i$ are the coefficients of the method, $s$ is the order of accuracy, and superscripts denote the time step [@HairerNorsettWanner1993].
With $\beta_0 = 1$ this corresponds to the forward Euler method ($s=1$) while $\beta_0 = 3/2, \beta_1 = -1/2$ gives second-order accuracy ($s=2$).
Applied to the momentum equation \eqref{eq:cont-vec}, this can be written as
\begin{equation}
\vfd{u}^{(n+1)} = \vfd{u}^{(n)} + \Delta t \sum_{i=0}^{s-1} \beta_i \left( \Op{F}( \vfd{u}^{(n-i)} ) - \Op{G} \vfd{p}^{(n-i)} \right)\,,
\end{equation}
where terms are grouped with the definition $\Op{F}(\vfd{u}) \equiv \Op{Adv}(\vfd{u}) + \frac{1}{\re} \Delta \vfd{u} + \frac{1}{\re} \vfd{b}_\Delta + \vfd{f}$ to simplify the notation.
Together with the continuity equation $\eqref{eq:cont-vec}$, the fully-discretized equations become
\begin{equation}
  \begin{aligned}
  \vfd{u}^{(n+1)} + \Op{G} \vphi^{(n+1)} &=
  \vfd{u}^{(n)} + \Delta t \sum_{i=0}^{s-1} \beta_i
  \Op{F}\left( \vfd{u}^{(n-i)} \right)
  \quad\text{and} \\
  \Op{D} \vfd{u}^{(n+1)} &= - \vfd{b}_{\Op{D}}
  \end{aligned}
\end{equation}
if we group the pressure contributions with $\vphi^{(n+1)} \equiv \Delta t \sum_{i=0}^{s-1} \beta_i\, \vfd{p}^{(n-i)}$.

fully-discretized system for explicit Runge–Kutta methods
––
Similarly, explicit $s$-stage Runge–Kutta methods can be written as $\vb{u}^{(n,i)} = \sum_{k=0}^{i-1} \left( \alpha_{ik} \vb{u}^{(n,k)} + \Delta t \beta_{ik} f(\vb{u}^{(n,k)}) \right)$ for $i=1,\dots,s$, with $\vb{u}^{(n,0)} = \vb{u}^{(n)}$ and $\vb{u}^{(n+1)} = \vb{u}^{(n,s)}$ [@GottliebKetchesonShu2009].
This is referred to as the Shu–Osher form [@ShuOsher1988], where the coefficients $\alpha_{ik}$ and $\beta_{ik}$ are not uniquely determined by the Butcher tableau of the method and can be chosen to minimize storage requirements.
In this case, the fully discretized equations can be written as
\begin{equation}
  \begin{aligned}
  \vfd{u}^{(n,i)} + \Op{G} \vphi^{(n,i)} &=
  \sum_{k=0}^{i-1} \left( \alpha_{ik} \vfd{u}^{(n,k)} + \Delta t \beta_{ik} \Op{F}\left(\vfd{u}^{(n,k)},\, t^{(n,k)}\right) \right)
  \quad\text{and} \\
  \Op{D} \vfd{u}^{(n,i)} &= - \vfd{b}_{\Op{D}}
  \end{aligned}
\end{equation}
with the definition $\vphi^{(n,i)} \equiv \Delta t \sum_{k=0}^{i-1} \beta_{ik}\, \vfd{p}^{(n,k)}$.

solution for each step/stage
- LU-decomposition of system
- sequential steps for solution
––
Each step or stage requires the solution of a system of equations in the form $\vfd{u} + \Op{G} \vphi = \vfd{a}$ and $\Op{D} \vfd{u} = \vfd{b}$.
This can be solved with the LU-decomposition
\begin{equation}
  \begin{pmatrix} \Op{I} & \Op{G} \\ \Op{D} & 0 \end{pmatrix}
  \begin{pmatrix} \vfd{u} \\ \vphi \end{pmatrix}
    =
  \begin{pmatrix} \Op{I} & 0 \\ \Op{D} & -\Op{D}\Op{G} \end{pmatrix}
  \begin{pmatrix} \Op{I} & \Op{G} \\ 0 & \Op{I} \end{pmatrix}
  \begin{pmatrix} \vfd{u} \\ \vphi \end{pmatrix}
    =
  \begin{pmatrix} \Op{I} & 0 \\ \Op{D} & -\Op{D}\Op{G} \end{pmatrix}
  \begin{pmatrix} \vfd{u}^\star \\ \vphi \end{pmatrix}
  = \begin{pmatrix} \vfd{a} \\ \vfd{b}\end{pmatrix}
  \,,
\end{equation}
where $\vfd{u}^\star \equiv \vfd{u} + \Op{G} \vphi$ has been introduced.
The steps to compute the solution are therefore
\begin{equation}
\begin{aligned}
  \vfd{u}^\star &= \vfd{a}
  \,, \\
  \label{eq:projection}
  \Op{D} \Op{G} \vphi &= \Op{D} \vfd{u}^\star - \vfd{b}
  \,, \quad \text{and} \\
  \vfd{u} &= \vfd{u}^\star - \Op{G} \vphi
  \,,
\end{aligned}
\end{equation}
where the second step requires solving a linear system.
Since the operator $\left(\Op{D}\Op{G}\right)$ has no coupling between different wavenumbers, Eq. \eqref{eq:projection} can be decomposed into $N_1 \times N_2$ tridiagonal systems of size $N_3$.
For $\kx = \ky = 0$ the system is singular due to the fact that the governing equations only include the gradient of the pressure and do not place any restrictions on the absolute magnitude of the pressure variable.
This system therefore has to be solved iteratively or the equations have to be regularized, e.g. by specifying an arbitrary value for one element of $\vphi$.
The current implementation relies on the Thomas algorithm [@QuarteroniSaccoSaleri2007] to solve the tridiagonal systems and includes the forward Euler and second-order Adams–Bashforth methods for time integration as well as the strong stability preserving Runge–Kutta methods SSPRK (2,2) and SSPRK (3,3) [@GottliebKetchesonShu2009].
Adding other explicit methods is straightforward, provided they can be formulated in a similar fashion.

## how is it implemented? {.outline}

Julia, MPI (vertical layers), FFTW
––
The simulation code is written in the Julia programming language, relying on the Julia bindings to the FFTW library [@FrigoJohnson2005] for fast Fourier transforms.
For parallelization, the domain is vertically split into up to $N_3$ blocks that are computed by separate processes exchanging information through the Message Passing Interface (MPI).


# Model validation {#sec:validation}

\begin{figure*}[!t]%
  \centering
  \includegraphics[width=\textwidth]{fig01.pdf}
  \caption{Error convergence for transient two-dimensional laminar flows. The left panel shows second-order convergence as the vertical grid resolution is refined for flows set up along the vertical direction and a randomly chosen horizontal direction.
  The other panels show first-, second-, and third-order convergence as the time resolution is refined for a Taylor–Green vortex set up along the two horizontal directions, in which case the spatial discretization is exact and the order of convergence of the time integration methods is measured. Grid lines show the formal order of convergence for each case.}\label{fig:validation-laminar}
\end{figure*}

## overview of validation methods {.outline}

methods
- automatic tests for individual terms
- automatic tests for laminar flows
- manual tests for turbulent flows
- randomize parameters in all cases
––
The validation efforts presented in this section aim to confirm that the numerical methods are implemented faithfully and that these methods produce physically relevant results.
To maintain this confidence as the code is inevitably modified, a focus is placed on automated tests that can be rerun after every change.
A set of automated unit tests verifies the expected order of accuracy when computing individual terms of the discretized equations for prescribed velocity fields and when applying the time-integration algorithms to ordinary differential equations.
A set of automated integration tests verifies that the solution to canonical transient two-dimensional laminar flows can be simulated with the expected order of accuracy.
Finally, fully turbulent flow solutions are computed and compared to published results produced with different codes.
These tests are not automated since they require significant computational resources and have no analytical solution to compare against so there is some degree of judgment required to evaluate the quality of the solution.

## individual terms: exact solutions {.outline}

The automated tests of individual terms make use of the fact that the implemented numerical methods are exact for certain velocity fields.
The diffusion term and the pressure solver are exact for a function that is the product of truncated Fourier series along horizontal dimensions and a quadratic polynomial in vertical direction.
The advection term is only exact for a linear function in vertical direction due to linear interpolations, although the term is still second-order accurate.
By constructing such a function with randomized parameters, each term can be computed numerically as well as analytically and matching values give a high degree of confidence in the correctness of the implementation.
Furthermore, we can verify the order of convergence when computing the terms at different grid resolutions for a velocity field that cannot be handled exactly by the implemented methods.
The time integration methods are verified in a similar way by solving ordinary differential equations that have analytical solutions with different time steps.
We also verify that the tridiagonal solver is exact for random inputs.


## laminar flow solutions {.outline}

laminar flow validation
––
To test the full solver including time integration, the automated tests include a number of laminar flow problems, currently the transient Poiseuille and Couette flows as well as decaying Taylor–Green vortices.
Numerical solutions computed at different resolutions are then compared to the analytical solution to ensure that the order of convergence corresponds to formal order of the numerical methods, as shown in Fig. \ref{fig:validation-laminar}.
Dimensional parameters such as domain sizes and velocity scales are again chosen randomly since parameters that are zero or unity can mask errors in the solution.
For Poiseuille and Couette flows, this includes the horizontal direction of the flow.
The two-dimensional Taylor–Green vortices are oriented both in horizontal and vertical planes.
For the former, the spatial discretization is exact so the test case verifies the order of convergence of the time integration method.

## general details about automated tests {.outline}

The above test cases have been verified in single-process (serial) mode as well as in multi-process (parallel) mode.
Multi-process tests are run both with a vertical resolution greater than and equal to the number of processes since those configurations sometimes rely on different code paths.
Since the tests are automated and run within minutes on consumer hardware, they can be rerun whenever changes are made to the code to ensure that any future version of the code still satisfies all the tested properties.

## turbulent flow validation {.outline}

turbulent validation flows: real flows with moderate complexity
- turbulence-resolving simulations too expensive for running interactively during development
- good to validate actual types of flows that code is made for
- validation flows chosen to be small enough that they can be re-run easily
- take about 2h using 32 cores on a system of the Intel Skylake generation
- points that could also be made: not only validates numerics but also physics
––
Turbulence-resolving flow simulations require substantially more computation and are therefore not included in the automated tests that are meant to be run routinely during code development.
The validation cases presented below are chosen such that they represent scientifically relevant flow systems while keeping the computational cost moderate.
Each case can be simulated in about two hours using 32 MPI processes on a single compute node of the Intel Skylake generation.

## validation dns {.outline}

rationale of dns validation
- idea: validate against published turbulent flow data that is relatively low cost so it can be re-run easily
––
The results of direct numerical simulations are not supposed to depend on the exact method used for modeling the flow, at least for lower-order flow statistics.
For relatively low Reynolds numbers, simulations have been run with many different codes and with a wide range of parameters such as domain sizes, aspect ratios, and grid resolutions, so the expected simulation results are well-established and have been validated against wind-tunnel measurements [@KimMoinMoser1987; @DelAlamoJimenez2003; @LeeMoser2015].

\begin{figure*}[!t]%
  \centering
  \includegraphics[width=\textwidth]{fig02.pdf}
  \caption{
    Direct numerical simulation (DNS) of a turbulent channel flow at $\re_\tau\approx180$, validated against data published by \citet{LeeMoser2015}.
    Mean profiles are shown for the streamwise velocity $u_1^+$, the advective transport $u_1^+ u_3^+$, the diffusive transport $\pdvinline{u_1^+}{x_3^+}$, as well as the production $\mathcal{P}^+$ and (pseudo)dissipation $\varepsilon^+$ of turbulent kinetic energy.
    The last panel shows contours of the premultiplied turbulent kinetic energy spectra $E_{ii}^+$ along the streamwise ($k_1^+ = 2 \pi \kappa_1/L_1^+$) and cross-stream ($k_2^+ = 2 \pi \kappa_2/L_2^+$) direction.
    The superscript ${}^+$ marks values in inner units, i.e. non-dimensionalized with the friction velocity $u_\tau$ and the kinematic viscosity $\nu$.
  }\label{fig:validation-dns}
\end{figure*}

details of dns validation
- validation against lee moser 2015, ret=180
- simulation setup: closed channel (ui=0 at vertical boundaries) flow with Reb=1/3.5e-4, ub held constant with dynamically adjusted pressure forcing, domain size 4πδ×2πδ×2δ, 255×191 horizontal modes, 96 vertical grid points with grid stretching η=0.97, starting from perturbed laminar solution, averaged over 17.5 turnover times after 3.5 spin-up
––
In Fig. \ref{fig:validation-dns}, we show a comparison of a closed-channel flow at $\re_\tau \approx 180$ with data published by @LeeMoser2015.
The friction Reynolds number $\re_\tau = u_\tau \delta/\nu$ is based on the half-channel height $\delta$ and the friction velocity $u_\tau^2 = \bigl. \nu \pdv{u_1}{x_3} \bigr|_{x_3=0}$ here.
The simulation is run with a bulk Reynolds number of $\re_\mathrm{b} = U_\mathrm{b} \delta / \nu = 20000/7$, where the vertically averaged bulk velocity $U_\mathrm{b}$ is held constant by a dynamically adjusted pressure forcing.
The solution is computed in a domain of size $4 \pi \delta$ in streamwise and $2 \pi \delta$ in cross-stream direction.
The velocity field is discretized with $255 \times 191$ Fourier modes at $96$ vertical grid points that are spaced according to a sinusoidal grid transform
$x_3(\zeta) = \delta + \delta \sin\left((2\zeta-1)\,\eta\,\pi/2\right) / \sin(\eta\,\pi/2)$
with $\eta = 0.97$.
The mean statistics computed over ${\sim}17.5$ large-eddy turnover times $T_\tau = \delta / u_\tau$ after a spin-up time of ${\sim}3.5\,T_\tau$ closely match the results from @LeeMoser2015.

## validation les {.outline}

rationale for les validation
- more difficult to validate since model differences can give different results, hard to distinguish legitimate model differences from implementation errors
- idea: validate against independent implementation of same model
––
For large-eddy simulation, validation is not as straightforward since results remain relatively sensitive to differences in the modeling approach and in the grid resolution.
To validate the new implementation, we limit ourselves to a comparison with a pre-existing Fortran implementation of the same physical and numerical models [@Giometto+2017a] and refer to previous publications for validation studies and discussions of limitations of the modeling approach [@PorteAgelMeneveauParlange2000; @Yue+2007a; @Yue+2008; @Giometto+2016].

\begin{figure*}[!t]%
  \centering
  \includegraphics[width=\textwidth]{fig03.pdf}
  \caption{
    Large-eddy simulation (LES) of a turbulent channel flow at $\re_\tau=10^8$ with an aerodynamically rough wall and a channel height of $h/z_0 = 10^4$, validated against a tried and tested Fortran code with the same numerical approach \citep{Giometto+2017a}.
    All values are non-dimensionalized with the friction velocity $u_\tau$ and the roughness length $z_0$.
    Mean profiles are shown for the streamwise velocity $u_1$, the resolved transport $u_1 u_3$, the subgrid-scale transport $\tau_{13}^\mathrm{sgs}$, as well as the production $\mathcal{P}$ and (pseudo-)dissipation $\varepsilon$ of resolved turbulent kinetic energy.
    The last panel shows contours of the resolved turbulent kinetic energy spectra $E_{ii}$ along the streamwise direction, premultiplied with the wavenumber $k_1 = 2\pi \kappa_1/L_1$.
  }\label{fig:validation-les}
\end{figure*}

details of les validation
- simulation setup: open channel flow (free slip above), ret=1e8, z0=1e-4δ, constant pressure forcing
––
In Fig. \ref{fig:validation-les} we show the results of this comparison for an open channel flow at $\re_\tau = 10^8$ driven by a constant body force $f_1$.
The friction Reynolds number $\re_\tau = u_\tau h/\nu$ is based on the channel height $h$ and the friction velocity $u_\tau^2 = h f_1$ here.
The lower surface is characterized by a roughness length $z_0$ that results in a non-dimensional channel height of $h/z_0 = 10^4$.
The solution is computed in a domain of size $2 \pi h$ in streamwise and $(4/3) \pi h$ in cross-stream direction.
The velocity field is discretized with $63 \times 63$ Fourier modes at $64$ equidistant vertical grid points.
The mean statistics computed over $200$ large-eddy turnover times $T_\tau = h/u_\tau$ after a spin-up time of $50\,T_\tau$ closely match for the two separate implementations.

## conclusion of validation {.outline}

conclusion
- validation presented here quite comprehensive
- confident in correct implementation of numerics, and that numerics can successfully represent physics
- goal not just validating the current state, but also make it easy to validate future development
- this should help not only adding new features, but also making changes to existing features and not get locked into design decisions that might prove suboptimal once more features are added
––
Combined, these validation efforts provide ample evidence that the implementation matches the mathematical formulation of the methods and that those methods are capable of accurately simulating flow physics, within the limitations of the physical models.
A comprehensive set of easily repeatable validation tests serves both to verify the current implementation and to ensure that future developments do not jeopardize correctness.
This should not only facilitate adding new functionality but also help making changes to existing functionality and avoid getting locked into design decisions that might prove suboptimal for future developments.


# Performance and scaling {#sec:performance}

## introduction: what is and isn’t discussed? {.outline}

introduction
- complete question: quality of solution, computational resources, time to solution
- broad question outside scope of this paper
- actual questions answered: quality = just available settings of solver, resources = just number of nodes of a fairly typical hpc cluster, time to solution
- in this case: 2× 40 core Intel Xeon Platinum 8380 ("Ice Lake") per node, intel omnipath network
––
Defined in a broad way, performance can be understood as the time required to obtain a solution at the required quality given the available computational resources.
We can examine how the time changes as a function of the required quality and the available resources (relative performance) or how fast different methods arrive at a solution for fixed quality and resources (absolute performance).
However, it is difficult to measure the overall quality of a turbulent flow simulation in a quantitative way since the system is chaotic and analytic solutions are not available.
Furthermore, there is a great diversity of computational resources that vary along important dimensions such as floating point operations per second (FLOPS), memory bandwidth and latency, network bandwidth and latency, and many more.
We therefore narrow the scope of the performance analysis to the question of the time required to obtain a solution given the specified simulation parameters and the number of compute nodes, as measured on a fairly typical high-performance computing (HPC) system.

which simulation parameters are relevant?
- time discretization: only number of RHS evaluation matters
- space discretization: resolution matters, with variation in vertical direction having strongest impact since domain is distributed between processes along that direction
- other parameters: only dns/les has strong impact, other settings such as type of boundary conditions, grid stretching, or forcing have virtually no impact on the amount of work of the solver
––
For the implemented explicit time integration schemes, the computational cost of a single evaluation of the right-hand side of $\dd\vb{u}/\dd t=f(\vb{u})$ fully characterizes the overall cost of a simulation, which is a simple function of the number of steps and the evaluations per step (i.e. stages of a Runge–Kutta method).
The computational cost of a single evaluation of $f(\vb{u})$ depends primarily on the number of Fourier modes and vertical grid points and whether a subgrid-scale term is modeled (LES) or not (DNS).
The impact of other parameters such as the type of pressure forcing (constant-flux vs. constant-force), boundary conditions, and grid transformations is imperceptible.

\begin{figure*}[!t]%
  \centering
  \includegraphics[width=\textwidth]{fig04.pdf}
  \caption{Performance and scaling of individual terms and full time step, as measured on Intel Xeon Ice Lake nodes of the Stampede2 system at the Texas Advanced Computing Center. Individual terms are computed at a resolution of $256 \times 256 \times \lambda N_\mathrm{p}$, where $\lambda \in \{1,2,4,8,16\}$ is the number of vertical grid points per MPI process and $N_\mathrm{p}$ is the total number of MPI processes. Dotted lines indicate weak scaling, dashed lines indicate strong scaling, and the grid lines correspond to perfect scaling. The last panel shows the overall performance for the resolution $256 \times 256 \times 1280$ (highlighted in orange on other panels) in comparison to a pre-existing Fortran implementation of the same numerical approach \citep{Giometto+2017a}. DNS performance is shown with the symbol $+$, LES performance with the symbol $\times$.}\label{fig:performance}
\end{figure*}


## relative performance {.outline}

overview of relative performance
- what is discussed?
- main observations: advection is dominant, pressure has worst scaling properties
––
To assess relative performance, Fig. \ref{fig:performance} shows the time required to compute the advection, diffusion, and pressure terms for different numbers of compute nodes and vertical grid points.
The figure displays both strong scaling, where the number of processes is varied for a problem of fixed size, and weak scaling, where the problem size is varied in proportion to the computational resources.
These results shows that the advection term contributes most to the overall cost while the pressure term exhibits the most problematic scaling behavior.

scaling of advection term
- almost perfect scaling, since cost is dominated by FFTs that are not distributed
- communication small part (otherwise strong scaling would look worse)
––
Computing the advection term is a global operation in horizontal direction but only involves neighboring nodes in vertical direction.
The bulk of the computational work consists of computing discrete Fourier transforms, which are local to each MPI process and scale as $\order{N \log N}$ where $N$ is the number of modes.
It appears that this cost dominates over the cost of communication, resulting in near-perfect strong and weak scaling.
When computing subgrid-scale stresses with a static Smagorinsky model, additional transforms are required and the cost increases to almost twice as much without affecting the scaling behavior.

scaling of pressure term
- no horizontal coupling (each mode separate), vertical coupling is global (thomas algorithm for tri-diagonal system) → strongest impact of vertical scaling as expected
- imperfect weak scaling: local work per layer the same but larger global size & more processes → impact of information having to travel up & down before completion?
- imperfect strong scaling, especially as 1 layer per process is approached: computation per communication goes down
––
Computing the pressure term has no data dependency between horizontal modes but is a sequential, global process in vertical direction (Thomas algorithm).
Horizontal modes can be processed in batches to stagger the sequential passes up and down the domain, where the size those batches is a tuning parameter that represents a trade-off between maximizing parallelism and minimizing per-batch overhead.
The resulting performance shows imperfect weak and strong scaling.
Scaling appears to improve when there are more vertical grid points per process, increasing the work-to-communication ratio.
While the overall cost appears to remain at most about a quarter of the cost of the advection term, it is possible that the two costs are even closer for some combinations of hardware configurations and simulation parameters, in which case it could be worth optimizing the batch size parameter of the pressure solver.

scaling of diffusion term
- no horizontal coupling (each mode separate), vertical coupling with neighbors only (local support of finite differences)
- weak scaling perfect (constant work to communication ratio since communication only with neighbors)
- strong scaling better at larger vertical sizes since ratio of work to communication is better, but still decent
––
Computing the diffusion term only involves neighboring vertical grid points and has no global data dependencies.
This results in near-perfect weak scaling.
Strong scaling is not quite perfect, which is explained by the fact that there is very little work to do for each grid point so the work-to-communication ratio is low.
This has no discernible effect on the overall scaling behavior however, as computing the diffusion term is always at least an order of magnitude less work than computing the advection term.

## absolute performance {.outline}

overview of absolute performance
- comparison with Fortran code
- answer practical question of substituting new code for existing one
- main observation: similar to Fortran code
- secondary observation: slight advantage of Julia code
––
To assess absolute performance, Fig. \ref{fig:performance} includes a comparison with a Fortran code that implements the same numerical methods [@Giometto+2017a].
While such a comparison does not answer the question of whether either code is making optimal use of the computational resources, it does respond to the practical question of whether there are any performance trade-offs when substituting the new code for a codebase that has been actively used for turbulence research for over two decades.
The comparison shows that the overall performance of both implementations is of a similar order of magnitude, with the new Julia implementation showing somewhat better scaling and significantly faster DNS performance.

details of performance comparison
- introducing excessive overhead was avoided without much effort and micro-optimization
- why slightly faster & better scaling? formulation in Fourier space avoids some FFTs, some communication–work trade-offs chosen differently, likely will inform some performance improvements in the Fortran code
––
It appears that the new Julia code has avoided introducing excessive overhead without much effort devoted to performance optimization.
That it even surpasses the performance of the Fortran implementation is likely explained by two factors.
First, the new code is formulated with the Fourier domain representation at its center, which makes it easier to avoid unnecessary Fourier transforms than in the physical-space formulation of the Fortran implementation.
Second, the new code makes different trade-offs between work and communication which appear to be more suitable for modern hardware.
Some of these insights will flow back to the Fortran implementation, reducing the performance discrepancy between the two codes.

## conclusion {.outline}

summary of findings
- main point: same order as Fortran
- main points: performance dominated by ffts, but pressure not too far off if many processes are involved & work per process is low, main goal to keep ffts to minimum & avoid introducing too much overhead
––
Overall, the performance characteristics of the new code are as expected.
The computational cost is dominated by the Fourier transforms necessary to compute the non-linear term while the pressure solver shows the least favorable scaling properties, and the overall performance is comparable to a Fortran implementation of an equivalent numerical scheme.
The analysis shows that for the implemented numerical scheme and current HPC hardware, performance is optimized by reducing the number and size of Fourier transforms and choosing an efficient implementation of the fast Fourier transform algorithm.
Other details matter less as long as the computational cost can be kept significantly below the cost of the non-linear term.

## future work {.outline}

future improvements should be weighted against benefits
- scale beyond 1 process per layer, either shared-memory parallelism or mpi-fftw
- move computation to GPU
- room for optimizations of communication, might become more important when parallelism is further increased
- but need to be careful with optimizations that make code harder to follow & change: computational cost at least quartic function of resolution, so large differences of performance required before there are meaningful differences in science that can be done
––
Future improvements are likely to focus on parallelism along horizontal coordinate directions, either through multi-threaded CPU code or through GPGPU computation, allowing the code to scale to larger systems.
There is also room for optimization in how communication is handled, which might become important if the work-to-communication is decreased through additional parallelism.
However, performance optimizations always have to be weighed against their impact on code simplicity and ease of adaptation.
Since the computational cost of flow simulations is a strongly non-linear function of the grid resolution, large performance differences are required for a practical difference in the scientific problems that can be tackled.

## On the suitability of Julia for high-performance computation

conceptual differences in Julia & Fortran
- two main differences to language like Fortran or C: types and memory
- types: usually not specified → when functions are called, the types of all arguments are known and version of the function is compiled for the specific combination of argument types; as long as all other types can be derived from the argument types (type stability), the compiler should be able to generate the same machine code as a language with static types
- memory: Julia has automatic memory management, i.e. it automatically allocates memory whenever needed and periodically interrupts the program to determine what memory is still in use and free unused memory (garbage collection); Julia provides good support for writing code without allocating intermediate values, but requires some attention to avoid allocating and deallocating lots of memory, which kills performance
––
Conceptually, there are two main differences in the performance characteristics of Julia compared to languages like Fortran, C, and \Cpp.
The first difference concerns the handling of types.
For a statically typed language like Fortran, the type of every variable is specified and therefore available to the compiler, which can use the information to generate efficient machine code.
In Julia, machine code is generated when a function is called for the first time, at which point the types of the function arguments are known.
For subsequent function calls with the same argument types, the compiled code is reused while a new copy of the function is compiled for calls with different argument types.
If the types of all variables inside a function can be derived from the type of function arguments, the compiler has access to the same information as for a statically typed language and can in principle generate the same or equivalent machine code.
The second difference concerns the handling of memory allocations.
For a language with (semi-)manual memory management like Fortran, memory is allocated and freed either manually with explicit commands or following deterministic rules.
In Julia, memory is automatically allocated whenever required for the operation that is performed and the program is periodically interrupted to examine which of the memory is still in use and which memory can be freed (garbage collection).
Vectorized Julia code written in a naive way often allocates large amounts of memory for intermediate results, but writing Julia code that operates in-place and avoids such allocations is relatively easy with experience.
Therefore, there is no *a priori* reason that code written in Julia should be significantly slower (or faster) than code written in Fortran, C, or \Cpp, although Julia does make it much easier to accidentally include expensive operations that result in poor performance.

philosophy comment, suitability of Julia
- strategy: collect measurements & make best trade-offs, avoid excessive overhead
- overhead that is intrinsic to Julia (GC & some function calls) should be minimal for this type of code, so much of the remaining differences would be due to how well performance optimizations in LLVM subsystem work compared to Fortran compiler
- Julia makes it easier to accidentally add overhead, requires regular checks and occasional fixes
- but Julia might help for measurements and clean organization that helps experimenting with methods, libraries, hardware (GPUs)
- Julia might be net advantage
––
Julia was chosen as implementation language for this code with the goal of improving the ease of development, hoping that the negative effect on performance could be kept minimal.
However, the experience so far has shown that the net impact on performance might even be positive.
Performance optimizations can be seen as a continuum from low-level (e.g. vectorized CPU instructions, optimal use of CPU cache) to high-level (e.g. choice of algorithms, speed–accuracy trade-offs).
At the lower end of this range, Julia relies on the LLVM compiler framework and a number of pre-existing libraries, making use of countless hours of optimization work, but Julia also makes it rather easy to write code it cannot optimize very well.
Maintaining close-to-optimal performance therefore requires regular measurements and occasional fixes.
High-level optimizations on the other hand require understanding the performance characteristics of different approaches and choosing the right one, often by measuring the performance of different implementations.
This type of optimization work benefits from the Julia language features and the ease of integrating packages from a growing ecosystem.
The impact of language choice on performance depends not so much on what optimizations are theoretically possible but which ones are simple enough that they are done in practice.
It is therefore possible that the choice of Julia will be a benefit rather than a drawback for the performance of this code over its lifetime.
