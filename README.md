Calculation of optical properties for MoS2 (later TMDs (transition metal dichalcogenides)) using matlab.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The calculations are based on the optical bloch equations for semiconductors. 

At first the electronic properties are computed with a tight binding simulation (Liu & Shan). Later the GW results from AG Wehling (Malte, Gunnar) will be used. These results are better in the entire Brillouin zone,
because the TB model is only fitted to the important symmetrie points (Gamma, K, K*). 
The electron spin is realized with an easy approach, a k-dependent Russel-Saunders coupling. (In the TB model the SOC is k-independent.)

The dipole transition are calculated using the Peierls approximation.

The coulomb interaction is implemented by using a fit formula (GW results) for the bare Coulomb matrix elements U. The screened Coulomb matrix elements are obtained with the two dimensional dielectric function,
which is also recieved through fitting to GW results.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The program is subdivided in a main program (MoS2_main.m) for control options. Calculations and plotting ist performed by subroutines in subfolders.