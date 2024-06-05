# KePASSA 2024 - On the validity of the Gaussian assumption in the Jovian system: evaluating linear and nonlinear filters for measurement-sparse estimation
The KePASSA2024 repository, along with the [GBEES repository](https://github.com/bhanson10/GBEES) contains the codebase accompanying the KePASSA 2024 presentation "On the validity of the Gaussian assumption in the Jovian system: evaluating linear and nonlinear filters for measurement-sparse estimation". Below is an in-depth summary explaining the proper usage of all the software included in said repository. In all, the codebase provides the computational framework necessary to perform and in-depth analysis of the results of comparing the results UKF, EnKF, and GBEES to a high-resolution Particle Filter (PF) used as the truth model.
These filters are applied to the Restricted Two-Body Problem (R2BP) and the Planar Circular Restricted Three-Body Problem __without measurement updates__. <br>

The code that performs the actual application of the UKF, EnKF, PF, and GBEES can be found in the [GBEES repository](https://github.com/bhanson10/GBEES), but this repository serves as the location for the MATLAB code used to create the analysis figures used in the presentation given at KePASSA 2024. From time to time I will reference the PDF of the slides to explain which script created which figure. Some of these MATLAB scripts require global functions I have defined in my [Functions repository](https://github.com/bhanson10/Functions). 

In alphabetical order...

## create_pdf_movie_pcr3bp.m
This script takes the probability distributions at the given epochs and creates an animation that maps how the uncertainty changes in the PCR3BP. Included in the animation are the UKF and PF. This animation was not used in the presentation, for time-restriction reasons. 

## create_pdf_movie_r2bp.m
See 'create_pdf_movie_pcr3bp.m', but for the R2BP. This animation was not used in the presentation, for time-restriction reasons.

## create_title_figure_pcr3bp.m
This script creates the figure on __Slides 1__ and __21__, by taking the probability distributions at given epochs and mapping them into a 3D space, comparing the UKF, PF, and GBEES. 

## create_title_figure_r2bp.m
See 'create_title_figure_pcr3bp.m', but for the R2BP. The image created by this script is used in __Slide 12__.

## jaccard_coef.m
This script takes two 2D/3D isosurfaces, with each isosurface either being a mean and covariance or an alpha shape, and approximates the Jaccard coefficient by using a discretized Monte Carlo grid. This function is used in many of the other scripts in the repository. 

## mc_vs_pf.m
This script creates a simple schematic that depicts the difference between a MC distribution and a PF distribution, used on __Slide 5__. 

## pcr3bp_filter_analysis.m
This script performs all of the comparison between filters for the PCR3BP. The results are used in __Slides 19__ and __20__. 

## pf_vs_gbees_resolution.m
This script plots the resolution of the PF vs. GBEES for the PCR3BP example. The resulting figure is utilized in __Slide 22__.

## r2bp_filter_analysis.m
See 'pcr3bp_filter_analysis.m' but for the R2BP system. These figures are used in __Slides 10__ and __11__.

## volume_overlap_slide.m
This slide creates the schematics used for explaining the process of taking a 3D point mass representation to a 3D isosurface. Figures are used in __Slides 8__ and __9__. 

<br><br>
For further information about code usage, please contact blhanson@ucsd.edu
