# PD-Control-SIPR-Lab
## PD Joint Space Control of a Soft Robot in PyElastica Python (Case where 3 segments are attached in series and kept in a plane and fixed at one end)

### Link to Progress Report/ Slides 

https://drive.google.com/file/d/1u_5f9a9J3cb-0_Ousl0klyvUvEig7lPT/view?usp=share_link

### Dynamic Model Used: Cosserat Rod Theory
The material of the rod was Silicone Rubber with Young's Modulus of $E = 3.79\ MPa$, density of $\rho = 1180\ kg/m^3$ and viscosity of $\eta = 12,500\ cps$.  
Using Experimental data, i.e., performing series of simulations to record bending angle and rate of bending angle data at different constant torque inputs, the stiffness and damping of the rod was estimated to be $k = 0.942$ N-m/rad and $\beta = 0.07$ N-m-s/rad. 


### Installation
See [2] for instructions on Installation and Documentation of PyElastica. You may need to downgrade (install earlier version) numba version apart from
```
$ pip install pyelastica
```
You can follow the error you get after running the pip install and running any PyElastica sample code to know which numba version to install and how.

### Results: Joint Space PD Control Kp = 3.3275 and Kd = 0.0075
### $q_1$, $q_2$, $q_3$ are estimated bending angles of rods 1, 2, and 3 respectively.
### $q_{des1}$, $q_{des2}$, $q_{des3}$ are desired bending angles of rods 1, 2, and 3 respectively.
### q1:
![fixed_pd_q_des_sin_cos_Kp_3_3275_q_1](https://user-images.githubusercontent.com/34472717/209483234-7a60c474-4f14-47fc-b29d-6aa20cd9653b.png)

### q2:
![fixed_pd_q_des_sin_cos_Kp_3_3275_q_2](https://user-images.githubusercontent.com/34472717/209483239-877afcc0-a385-4281-9e77-ccd802a79914.png)

### q3:
![fixed_pd_q_des_sin_cos_Kp_3_3275_q_3](https://user-images.githubusercontent.com/34472717/209483247-9593443b-56ce-46fc-9db8-a9b881708e5a.png)

### Positions
![fixed_pd_q_des_sin_cos_Kp_3_3275_Positions](https://user-images.githubusercontent.com/34472717/209483250-45b5cf22-2e0e-4cda-b323-aaf25900abaf.png)


### References
[1] Tekinalp, Arman and Kim, Seung Hyun and Parthasarathy, Tejaswin and Bhosale, Yashraj, $\textbf{PyElastica: A computational framework for Cosserat rod assemblies}$, https://github.com/GazzolaLab/PyElastica, 2022

[2] https://docs.cosseratrods.org/en/latest/


