# Code for "Bridging developmental and statistical approaches to variation and evolution" 
by Lisandro Milocco and Tobias Uller

The folders in this repository are named after the figures of the paper.

---

## Figure 1
Contains the script **`Figure1.m`**, which generates the time-series data for the bistable switch and calculates both the regression-based average effect and the sensitivity.

---

## Figure 2
Contains the script **`Figure2.m`**, which generates the histogram.

The data for the histogram are generated with **`generateData_Figure2.m`**, which allows selecting which parameter varies (population size, developmental noise, or measurement noise) while keeping the others fixed.

The adaptive Kalman filter implementation is in **`adaptiveKF.m`**.

---

## Figure 3
The script **`Simulate_G_P_Matrix_Development.m`** regenerates populations for different values of the parameter *p* in Figure 2, which is the minor allele frequency underlying parameter $\lambda_2$.

For each value of *p*, the script generates population genotypes and then calculates phenotypes using the ODE model of the toggle switch gene network. The P and G matrices calculated here are then read by **`Figure3_Ellipses.m`**, which generates the ellipses and calculates the angles between $G_{max}$ and $P_{max}$.

To test evolutionary predictions using G or P matrices, the script **`LandePrediction_Generate.m`** generates a population with genotypic frequencies given allele frequency *p* (see panel C of Figure 3), applies a single round of selection toward an optimum of $(4,4)$, and predicts the change using the G matrix and P as a proxy.

**`LandePrediction_Analysis.m`** calculates prediction errors and generates the plots in panel C.

---

## Figure 4
For this figure, the **Tooth Model** (see original publication for details) was used.

Small perturbations in developmental parameters were introduced, and sensitivities were calculated by regression (see main text). The script **`Figure4.m`** plots the sensitivities, calculates the angles between them, and shows that these angles are smaller than those between random vectors.
