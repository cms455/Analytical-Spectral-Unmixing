<h1 align = 'center'> Optimal and Fast Algorithm for Wavelength Selection for Spectral Unmixing </h1>
<p align="center">
 <img src="/README_images/collage_v2.png" width="600" height="400">
</p>

## Summary: 
We present a novel method for determining the optimal wavelengths to perform spectral unmixing in photoacoustic experiments.
Our method utilizes an implementation of the Bourgain-Tzafriri algorithm from [ref Tropp] to perform subset-column selection
on the absorption coefficient matrix to construct submatrices of optimal wavelengths. The algorithm performs significantly better
than current methods and is able to find better wavelength selections more quickly in many cases. In fact for matrices that optimal
solutions can be reasonably calculated, the algorithm is able to find the same the true optimal results as found by brute force searches
with a high probability. The algorithm generalizes better for more complicated cases, can be adjusted to look for other qualities
such as condition number, matrix norm, and can be tuned to find good approximations for specific time constraints.
## Notable Features:
- Implementations of Brute-force Algorithm, Luke Algorithm (greedy maximazation of the least sigaulr value), Bourgain Tzafriri Algorithm, modified-BT algorithm, and more!
- Code for test cases
- Examples on Hb, HbO, ICG, ...
- Database of absorption coefficient spectrum
## Table of Contents:
- [Installation](#installation)
- [Examples](#examples)
## Installations:
### Requirements
- Matlab 2024b
## DOxygen Documentation
For more thorough documentation there is a generated wiki genreate through DOxygen. Open Docs/html/index.html and search functions. 
<p align="center">
 <img src="/README_images/Doxygen_picture.png" width="1000" height="400">
</p>


## Getting Started:


## Acknowledgements
Special thanks to David Veysset for mentoring me through this process. Thank you to Bhaskara Rao Chintada as well for help with the delay and sum algorithm.
Thank you to the Bouma Lab for supporting me for the past few months!



