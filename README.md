# Covid-19 simulation usign network analysis

This repo contains a simplified simulation of the spread of Covid-19. It is based on the epidemicalogical SEIR model making use of random geomatric graph instroduced by Waxman (1988) to mimick human social interaction. The simulation is able to capture the effect of social disctacing, lockdowns, herd immunity in addition to large and small (social) events.

1) `main.m` : Runs the simulation and contains options to test the various scenarios mentiond above. Also produces 
2) `SEIR.m` : Contains the SEIR algorithm.


Example results: 
### No intervention at all.
![N_1000_rep_10_T_60_r0_2 2](https://user-images.githubusercontent.com/23720435/148948664-5877dfe6-3513-4446-ae2a-6ade6c7864eb.png)

### Mandatory social distancing
![N_1000_rep_10_T_60_r0_2 2_social_distancing_at_0_1](https://user-images.githubusercontent.com/23720435/148948760-e9be06e3-3611-428b-ba8f-7db681bd67e0.png)


### Early governent intervention, but resctictions lifted too early.

![N_1000_rep_10_T_60_r0_2 2_lift_events_at_0_1](https://user-images.githubusercontent.com/23720435/148948810-b4661db1-8271-4dd4-bf69-e6d42ff0f1ff.png)
