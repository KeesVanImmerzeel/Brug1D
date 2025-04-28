# Simulate time-dependent head and flow in a semi-infinite aquifer adjacent to open water.

C.H. van Immerzeel

24/4/2025

Simulate time-dependent head and flow in a semi-infinite aquifer adjacent to 
open water (for example a river), where the boundary condition changes from t = 0 according to a 
specified course.

The equation used describes one-dimensional flow for a situation with 
spatially constant kD and S, without recharge (precipitation, evaporation, leakage).

The response to a sequence of stress values (a) may be studied (simulation option).

![](https://github.com/user-attachments/assets/48ca5d4c-1f16-4f37-b254-1eb2972da429)

![](https://github.com/user-attachments/assets/a2abc864-4ecc-48bd-8d36-50992e3bf623)

![](https://github.com/user-attachments/assets/c7b7ab7f-a3e2-4bc7-aff5-5cec2c32def5)

![](https://github.com/user-attachments/assets/f8b6d4f8-76a5-49a1-8064-9a89808d5961)

![](https://github.com/user-attachments/assets/c9568185-5a05-4aa8-8f04-f14444804f52)


## Link to the app
<https://sweco.shinyapps.io/Brug1D/>

## Source code
R-source code of the app:

<https://github.com/KeesVanImmerzeel/Brug1D>

## Simulation
A sequence of stress values (a) may be uploaded. For this purpuse, prepare a spreadsheet with two columns ('Time' and 'a').  


## Numerical stability
Be aware that exotic input values may lead to numerical instability. In that case, no results or plots are presented. Instead, an error message appears.


## References
- Olsthoorn, T.N. (2006). Van Edelman naar Bruggeman. Stromingen 12 (2006) nummer 1, blz. 5-12. 
  <https://edepot.wur.nl/13730>
- Veling, E.J.M. (2006). Over de erfc-functie in ‘Van Edelman naar Bruggeman'. Stromingen 12 (2006) nummer 2, blz. 56-57. <https://www.nhv.nu/wp-content/uploads/2020/06/2006-2_brieven.pdf>
- <http://www.grondwaterformules.nl/index.php/formules/waterloop/peilverandering>  
  
## Van Edelman naar Bruggeman

![](https://github.com/user-attachments/assets/80a9dcff-c016-45ea-8948-4a08b783c1eb)
![](https://github.com/user-attachments/assets/8cc312f2-43a4-4647-bfcb-fd514d459d4c)
![](https://github.com/user-attachments/assets/5eb4ca1d-3de4-4190-81dd-03795d0a0085)
![](https://github.com/user-attachments/assets/2ce390ed-9c19-44ac-8eee-bfce2e978bfa)
![](https://github.com/user-attachments/assets/ee18954a-5d0b-4bfd-bec1-cd117e1332fa)
![](https://github.com/user-attachments/assets/91050b06-3978-4e05-a933-11ebee677b4b)
![](https://github.com/user-attachments/assets/c02509db-33fd-497e-a338-4d80ea360cd8)
![](https://github.com/user-attachments/assets/37877edf-6de4-4d39-9706-6a271ab6605c)



