# sphericity
Step-down method for assessing sphericity of a VCV matrix

Initially sparked by Kirkpatrick's effective dimensions (Kirkpatrick, Genetica 2009), which is a function of both the evenness of variances 
and the strength of covariances. This procedure breaks the sphericity into either component to assess where deviations from spherical emerge 
from given a matrix of interest.

Code here uses matrices estimated using MCMCglmm, and uses the posterior distribution to evaluate divergence from spherical. 

Plots are used to assess overlap visually.
