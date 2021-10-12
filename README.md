# Analysis of Microsoft SEAL in context of the SecuSim projekt

## Description

The SecuSim project deals with secure cloud-based simulations. 
The benefit of cloud computing may evade the acquisition of expensive hardware that often stands in the beginning of a simulation project. In many cases 3D models shall
be simulated that are confidential. Such cases require a comprehensive consideration of the security of cloud computing. Data in a cloud are publicly accessible, thus a third party can try to get access.
Most cloud providers ensure a certain degree of security, which might not be sufficient for some confidential simulation data. Therefore, it was investigated how simulations can be performed on encrypted data, 
such that the confidential data never exist unencrypted on the cloud. One possible way to achieve this is the technology of Homomorphic Encryption (HE). This technology makes it possible
to perform mathematical operations on encrypted data and obtain the result (after decryption that is) one would get if the operation would have been performed on the unencrypted data.
The different HE schemas support different types of mathematical operations and different number of sequential operations. The encryption comes with a memory overhead that is reduced if vectors of data
are encrypted. For that reason, two libraries were investigated, the SEAL and HEAAN library (https://github.com/SanJomyi/SecuSim_HEAANDemo). Both libraries support fixed point arithmetic, which requires special handling if it shall be used for a simulation.
With this code the memory overhead of the encrypted data and the performance of a matrix-vector multiplication were analysed. It was planned to implement the encryption into the matrix based FIT TD solver (https://github.com/SanJomyi/EMW_FIT_Solver), which 
is being developed for the OpenSource simulation platform OpenTwin (link is following).


## Conclusion
The result of this investigation is that SEAL is currently unsuitable for simulations. Detailed information can be found in the publication (in progress). Nevertheless, HE seems very promising and SEAL is well documented.
Therefore, this software project can be used to analyse the suitability of upcoming library versions for simulations.

## Dependencies

* SEAL VERSION 3.6.6

## References

* SEAL project: https://github.com/microsoft/SEAL
* Implementation of matrix vector multiplication: https://github.com/MarwanNour/SEAL-FYP-Logistic-Regression
