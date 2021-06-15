# Basic ecosystem box model with N isotope

This ecosystem box model simulates the basic NPZ (nitrate, phytoplankton, zooplankton) food web dynamics with N isotope. The N isotope model simulates the evolution of the two stable N isotopes, <sup>14</sup>N and <sup>15</sup>N in all reservoirs. The fluxes of <sup>14</sup>N and <sup>15</sup>N between N reservoirs are simulated based on Rayleigh isotope fractionation kinetics for irreversible reactions in a closed system. 

For non-fractionating physical processes, N fluxes are calculated as a concentration-weighted volume balance. For fractionating processes, such as nitrate assimilation and zooplankton grazing, the flux of <sup>14</sup>N is simulated following the Michaelis Menten dynamics, while the corresponding flux of <sup>15</sup>N depends on the flux of <sup>14</sup>N and the isotope effect of the process. The fraction of the grazed phytoplankton that is digested to be zooplankton biomass has an isotope effect, while the remainder that is exported as fecal pellets has no isotope fractionation.

<p align="middle">
  <img src="https://user-images.githubusercontent.com/47376014/122094595-57879500-cdda-11eb-9b0f-be88ccf26d6b.PNG" width="400">
</p>

### This model is sensitve to the prescribed zooplankton maximum grazing rate, gmax
1. With maximum grazing rate, gmax = 0.6 /day, the model reaches a steady state
<p align="middle">
  <img src="https://user-images.githubusercontent.com/47376014/122096498-a3d3d480-cddc-11eb-908e-c9ef121eb255.png" width="400"> 
  <img src="https://user-images.githubusercontent.com/47376014/122096337-738c3600-cddc-11eb-9d44-080029fb8c76.png" width="400">
</p>


2. With maximum grazing rate, gmax = 0.7 /day, the model is periodic
<p align="middle">
  <img src="https://user-images.githubusercontent.com/47376014/122096517-a7675b80-cddc-11eb-9e6c-49f3fbdf75a1.png" width="400"> 
  <img src="https://user-images.githubusercontent.com/47376014/122096533-aafae280-cddc-11eb-8071-d57bd2f36866.png" width="400">
</p>
