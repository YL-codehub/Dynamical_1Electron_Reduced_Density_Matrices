# Dynamical 1electron RDM: model fitting under Matlab
by Yoann Launay*, initiated by Benjamin De Bruyne*.
Former CentraleSupélec students at Paris-Saclay University.

Feel free to contact for help and please cite the following articles when using the software for your research:  

"N-representable 1-electron reduced density matrices at non-zero temperatures", Y.Launay
and J-M. Gillet, Acta cristallographica B, 2021.

"Inferring the one-electron reduced density matrix of molecular crystals
from experimental data sets through semidefinite programming", B. De Bruyne and J-M Gillet, Acta cristallographica A, 2020

Note that all calculations are performed in atomic units (a.u).

## Installation

### Packages
- You need to either re-install mosek and yalmip completely in this repertory or to add them to your path
- You need to get a Mosek license and put it in your home dir as explained in https://docs.mosek.com/9.0/install/installation.html (choose free 1 year academic license)
- For any dilation, you need the Optimization Toolbox which can be installed at the installation of Matlab.

### Before you start
You might need to add to path every file of the required packages (especially
Yalmip interface and Mosek packages). If necessary, uncomment and modify
this command line to go over it (change the source path) :

```
addpath(genpath('D:\Cristallographie quantique\cleanRDM\RDM_SDP'));
```

If you keep the packages already installed by your predecessor, you might have an "invalid MEX-file" issue see then section "troubleshooting" at https://docs.mosek.com/9.2/toolbox/install-interface.html to
update without install again. For example you may need the following command line
to be uncommented in main.m:

```
setenv('PATH', [getenv('PATH') ';your_path\mosek\8\tools\platform\win64x86\bin']);
```

### How to start
You need to uncomment or edit main.m with your preferences. After that, all is left is running it. 

See below for some possible runs.

## Functionalities

Note that functions are documented (inputs, outputs, goal and methods).
In Matlab you can call : help('functionName') in the command window.

### EXAMPLE 1 : Import expectation values and add noise

Different files are needed to read data on (see readCrystal.m's doc).
An example at 0K (see readCrystal.m's doc for more details): SF data in the Ewald sphere of radius 8 and DCPs under 6 u.a momentum with a resolution of 10 points per a.u (max 100 but impossible in real data).(3 directions in our default Crystal outputs).

```
[SF,CP] = readCrystal(primVectors,'InputFiles/prop_ter.d3','InputFiles/prop_output_ter.out','InputFiles/PROF_ter.DAT',0,8, 6, 10); 
```

Or same example at 100K 
```
[SF,CP] = readCrystal(primVectors,'InputFiles/prop_ter.d3','InputFiles/prop_output_ter_100K.out','InputFiles/PROF_ter.DAT', 0,8, 6,10);
```

You can also add some noise to these two with a given amplitude, here 3by default:
```
[SF,CP] = addNoise(SF,CP);
```
Note that:
- You'll soon face many problem to have a correct consideration of both
SFs and DCPs in your Chi^2. It is important to keep enough DCPs compared
to SFs. As in the example, 7of DCPs point is correct. This even more
true when going in 100K surveys. This is why you usually need to keep a
small Ewald sphere for your SFs (<= radius 9, or a spherical chelle for higher values 9<= ..<=12 for instance)
- 100K files are hand-made : 'prop_ter.d3' and
'PROF_ter.DAT' do not change from 0 to 100K but
'prop_output_ter_100K.out' is 'prop_output_ter.out' computed at 0K where
we just replaced the static SFs by their dynamic values computed with
Crystal (in order to keep DCPs at 0K).
- The level of noise is usually set at 1to get better results at 100K.
- Any variable can be saved and load later thanks to 'save' and 'load'
functions (see Matlab's doc).

### EXAMPLE 2 (S0): Refine a population matrix and plot the associated 0K maps

From previous variables:
```
[solP,C] = optimizeP_quadratic(atomicOrbitals,R,T,SF,CP);
```

Some info appears in the Command Window during the refinement :
- the status of the refinement(s) 
- the dialog box of Mosek solver
- Performances (proportion of SF and CP in data or in Chi^2, Chi^2,
CHi^2/ndf)
Note that : 
- the bigger the data, the more time it takes, especially for the creation of predictors.
- for a reasonable choice of data, the solver only needs a few seconds to conclude.

You can now plot the results:
```
plotALL(solP,atomicOrbitals,R,T,SF,CP,volume,'InputFiles/DENSPatoNew.dat'); 
```
Note that the desorthogonalization of P is included here.

The order of the figure is the following :
- R-factors displayed in the Command Winfdow are equal to sum(abs(abs(Data)-abs(Predicted)))/sum(abs(Data));
- Directional Compton Profiles (as many figures as there are directions, 3
by default). Dashed blue curve = refined. Continuous red curve =
pseudo-experimental.
- Difference between previous dashed blue and continuous red ones for
each directions.
- Anisotropies with the last direction of the DCPs. Dashed blue curve = refined. Continuous red curve =
pseudo-experimental. (Warning, legend has to be adapted if directions
change, see plotCP.m file).
- Deformation density or the difference between the refined density (diagonal terms of the 1RDM) and the pseudo-atom density map at 0K. 
(Warning, there is no such map, like 'InputFiles/DENSPatoNew.dat', for a given temperature)
- Refined Fourier map. Contours at +1 (blue) or -1 (red, due to truncation) x0.01x2^n (n=0 to 20) a.u^-3. 
- Pseudo-experimental Fourier map. Contours at +1 (blue) or -1 (red, due to truncation) x0.01x2^n (n=0 to 20) a.u^-3. 
- Difference of the previous maps : pseudo-data minus refined Fourier. Contours at 0 (green) or +1 (blue) or -1 (red) x0.0001x2^n (n=1 to 20) a.u^-3. 
maps.
- 1-RDM map twice along CO2 axis. Contours at +1 (blue) or -1 (red) x0.01x2^n (n=0 to 20) a.u^-3. 
See articles for more precision on figures.


### EXAMPLE 3 (S0 x N): Launch a Monte-Carlo campaign for the 0K SDP method and visualize results
N refinements are made based on random and different noise but with same amplitude:

```
N = 10;
[SFmc,CPmc] = readCrystal(primVectors,'InputFiles/prop_ter.d3','InputFiles/prop_output_ter.out','InputFiles/PROF_ter.DAT',0,8, 6,10 ); 
SFmc_n = [];
CPmc_n = []; %noise storage (0 index is non-noisy)
SOLS_P = [];
for i=1:N
    disp('iteration');
    disp(i);
    [SFmc_n2,CPmc_n2] = addNoise(SFmc,CPmc);
    SFmc_n = [SFmc_n SFmc_n2];%concatenate
    CPmc_n = [CPmc_n CPmc_n2];%concatenate
    [solP_mc2,C_mc2] = optimizeP_quadratic(atomicOrbitals,R,T,SFmc_n2,CPmc_n2);
    SOLS_P = [SOLS_P solP_mc2];
end

plotMonteCarloRDM(SOLS_P,atomicOrbitals); plot mean and std maps
```

### EXAMPLE 4 : Save or import a 1-RDM and plot

Don't forget the desorthogonalization of the population matrix !
```
write_output(atomicOrbitals, C.'*solP*C, T, R, 'ExitMatlab.txt')
[atomicOrbitals, solP, R, T] = read_exit('ExitMatlab.txt');
```

Note that :
- the format of the file can be changed by editing write_output
- to read the 0K Gamess calculus with 3-21(d) basis set, you have to
%uncomment the permutation of orbitals in read_exit code.
- Such a format is also used by my software written in Python language.

### EXAMPLE 5 (S100): Refine a 100K population matrix 
Combine what we just learned
```
[SF,CP] = readCrystal(primVectors,'InputFiles/prop_ter.d3','InputFiles/prop_output_ter_100K.out','InputFiles/PROF_ter.DAT', 0,8, 6,10);
[SF,CP] = addNoise(SF,CP,0.01); 1%-noisy
[solP,C] = optimizeP_quadratic(atomicOrbitals,R,T,SF,CP);
plotALL(solP,atomicOrbitals,R,T,SF,CP,volume,'InputFiles/DENSPatoNew.dat');
```

### EXAMPLE 6 (S100CD100S100): Refine both the basis set and population matrix at 100K 
Introduce intermediary refinements
```
[SF,CP] = readCrystal(primVectors,'InputFiles/prop_ter.d3','InputFiles/prop_output_ter_100K.out','InputFiles/PROF_ter.DAT', 0,8, 6,10);
[SF,CP] = addNoise(SF,CP,0.01); 1%-noisy
[solP,C] = optimizeP_quadratic(atomicOrbitals,R,T,SF,CP);
atomicOrbitals_final = atomic_least_squares(atomicOrbitals,solP,R,T,SF,CP); or use valence_least_squaresà améliorer et fusionner
[solP2,C2] = optimizeP_quadratic(atomicOrbitals_final,R,T,SF,CP);
```

### EXAMPLE 7 (S100DW100): Have a try on Debye-Waller terms refined on 100K-refined population matrices  

```
[SF6,CP6] = readCrystal(primVectors,'InputFiles/prop_ter.d3','InputFiles/prop_output_ter_100K.out','InputFiles/PROF_ter.DAT',0,8, 6,10 ); 
 ```

Get a spherical shell with no DCPs for DW refinement:
```
[SF6bis,CP6bis] = readCrystal(primVectors,'InputFiles/prop_ter.d3','InputFiles/prop_output_ter_100K.out','InputFiles/PROF_ter.DAT', 8.001,12, -0.1,10);

[solP, C] = optimizeP_quadratic(atomicOrbitals,R,T,SF6,CP6);
B = optimizeB(atomicOrbitals,solP,R,T,SF6bis)
```

You'll see that unefficient DW ===> Strong S100 !

### EXAMPLE 8 (G0DW100S^S100): Refine a static map from Debye-Waller and SDP refinement

Import Gamess results with 3-21G(d) basis set
```
[atomicOrbitals_G, solP_G, R, T] = read_exit('CO2gamess.txt',true);
Visualize RDM and observable, comparison with Crystal ones
[SF,CP] = readCrystal(primVectors,'InputFiles/prop_ter.d3','InputFiles/prop_output_ter.out','InputFiles/PROF_ter.DAT',0,15, 6, 10); 
plotALL(solP_G,atomicOrbitals_G,R,T,SF,CP,volume,'InputFiles/DENSPatoNew.dat');
```
and the 100K data
```
[SF6,CP6] = readCrystal(primVectors,'InputFiles/prop_ter.d3','InputFiles/prop_output_ter_100K.out','InputFiles/PROF_ter.DAT',0,8, 6,10 ); 
[SF6bis,CP6bis] = readCrystal(primVectors,'InputFiles/prop_ter.d3','InputFiles/prop_output_ter_100K.out','InputFiles/PROF_ter.DAT', 8.001,12, -0.1,10);%no CP?
Add noise or not:
[SF6,CP6] = addNoise(SF6,CP6);
[SF6bis,CP6bis] = addNoise(SF6bis,CP6); %no CP data
```

Refine DW and Pop. Matrix
```
B = optimizeB(atomicOrbitals_G,solP_G,R,T,SF6bis)
[solP,C] = optimizeP_quadratic(atomicOrbitals_G,R,T,SF6,CP6,B);

plotALL(solP,atomicOrbitals_G,R,T,SF6,CP6,volume,'InputFiles/DENSPatoNew.dat',B);
```

Now is your turn to combine ... Other possibilities with softwares in 'Explorations' or see article.

## Hierarchy for future developers
<p align="center">
  <img src="https://github.com/YL-codehub/Dynamical_1Electron_Reduced_Density_Matrices/blob/main/SoftwareArchitecture.png" alt="Header" style="width:75%;"/>
</p>
