# ReweightedPartialSparseRecovery
### Target equation: the Aizawa equations
                     udot_{1} = (u_{3}-0.7)*u_{1} - 3.5*u_{2};
                     udot_{2} = 3.5*u_{1} + (u_{3}-0.7)*u_{2};
                     udot_{3} = 0.6+0.95*u_{3} - u_{3}^3/3 -...
                                (u_{1}^2+u_{2}^2)*(1+0.25*u_{3}) + 0.1*u_{3}*u_{1}^3

### Auxiliary functions:
    * library.m (The Candidate Functions)
    * derivative.m (Numerical derivation) 
    * RWPSTRidge.m (Reweighted Partial STRidge Algorithm)
    * aizawa.m (ODE to Test)
#### Authors: Xiaofan Lu, Huimei Ma, and Linan Zhang
#### Data: September 23, 2022
