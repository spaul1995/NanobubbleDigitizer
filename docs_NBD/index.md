# PyNano 

NanobubbleDigitizer is a MATLAB based nanopore data analysis interface, which can be used to analysis the current-time data of bubble generation in solid-state nanopore.	


## Installation 

Download the repository files as a .zip and extract in a folder in your PC. MATLAB (Mathworks.Inc) should be pre-installed to run the code.


## NanobubbleDigitizer user manual

![ui1](Slide2.jpg)

Fig.1 (a) Current time trace of nanopore bubble blockage. By running the algorithm, we have automated the detection of the nucleation point and collapse point. From these two points, 4 features are extracted. The current values at the nucleation point and collapse point gives nucleation current and collapse current respectively. The time separation between the nucleation point of a bubble and the collapse point of the previous bubble gives the waiting time for bubble nucleation. The time separation between the nucleation point and collapse point of the same bubble gives the blockage duration. (b) The heat map shows the correlation coefficient between the 4 features. (c) The histogram shows two peaks related to blockage duration, indicating two types of nanopore bubbles. The second peak which happens at a higher value of blockage duration indicates homogeneous bubbles nucleating at the pore center while the first peak at lower value of blockage duration indicates heterogenous bubbles on the cylindrical pore surface.


## Website Under Construction