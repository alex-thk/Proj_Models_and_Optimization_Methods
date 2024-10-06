# Optimal Control of an F/A-18 Hornet Landing Procedure

## Project Overview

This project focuses on the optimal control of an F/A-18 Hornet fighter jet during its landing and braking on an aircraft carrier. The system is nonlinear and was linearized to be solvable using an **automatic control approach**. The solution ensures safe landing and braking under constraints like pitch angle, touchdown velocity, and wind gusts.

The process is split into two phases:
1. **Landing Phase**: Control of the descent, managing aerodynamic forces and wind gusts.
2. **Braking Phase**: Control of deceleration after touchdown using a tailhook mechanism.

## MATLAB Scripts

The following MATLAB scripts are used for the simulations:
- `controllo_jet_land.m`: Simulates the landing phase of the jet.
- `controllo_jet_break.m`: Simulates the braking phase of the jet.
- `build_qp.m`: Constructs the quadratic program to solve the optimal control problem.

## Full Explanation

For a complete description of the problem and the detailed solution, please refer to:
- **Papi_proj-5.pdf**: Detailed solution and explanation.
- **Progetto-4.pdf**: Initial problem setup.

