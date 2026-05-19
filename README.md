# Temperature Distribution in a 2D Composite Wall Using the Finite Element Method

This project presents a finite element analysis of steady-state heat conduction in a two-dimensional composite wall. The temperature distribution is solved numerically using the Finite Element Method (FEM) with four-node quadrilateral elements.

The project was developed as part of the ME622 Finite Element Method course at the New Jersey Institute of Technology.

## Objectives

The main objectives of this project are:

- To model steady-state heat conduction in a 2D composite wall.
- To implement the finite element method using Q4 quadrilateral elements.
- To apply Dirichlet and Neumann boundary conditions.
- To compute the temperature distribution across the wall.
- To post-process and visualize heat flux and temperature contours.



## Repository Structure

```text
temperature-distribution-2D-wall/
│
├── main_composite_wall_2D.m
├── mesh_structured_Q4.m
├── element_Q4_K.m
├── apply_dirichlet.m
├── apply_neumann_edges.m
├── postprocess_flux.m
├── plot_results.m
├── report/
│   └── Temperature_Distribution_2D_Composite_Wall.pdf
└── README.md
