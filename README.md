# Unit_Ops_Distillation_Col
Simulation of plate distillation column using custom thermodynamical model

Things to note:

1) Code takes a long time to run, up to 5min to converge to exact result

2) NRTL.m is the Themodynamical Model that predicts composition given a temperature and pressure (It uses Raoult's Law to find an initial approximate solution to converge)

3) Could be optimized furthur by changing the solver to non-symbolic!

How it works:

Key in parameters for distillation column ---> Get distillation column temperature and concentration profile!
