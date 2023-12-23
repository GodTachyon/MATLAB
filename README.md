# MATLAB - Data Analysis and Uncertainity Assignment (CFD)

This is **Aswath Ashok**'s repository for Data Analysis and Uncertainity module 2023-24. 
The Assignment files are in the assignment folder.

Changes made

15/11/23 - 
Created Repository and required folder for assignment titled "Assignment".
Uploaded files from Dr. Jozsa's matlab scripts.

12/5/23 - 
Uploaded the Assignment Folder onto the Repository and sent and invite to Dr. Tamas Jozsa to act as an observer.

21/12/23 -
Covered Literature Review, visualising coherent structures on Tecplot360 and established Q-criterion for the 3 velocity components.
This serves as a reference for Task 2.

For Task 3, Updates to (Assignment/assignment.mlx) has been modified to read data samples, calculate the mean and standard deviation of the velocity components and ...
plot them on figures. The gradient has also been calculated and plotted. The accuracy of the gradient values is in doubt and further investigation is needed.

The (Assignment/symbolic.mlx) file has been newly added today. This is an attempt at Task 5 in which I attempt to use Symbolic functions to obtain the wall - normal 
mean velocity componenet and plot it.

23/12/23 - 
[https://github.com/GodTachyon/MATLAB/blob/main/Assignment/FFT.m](FFT.m) has been uploaded. This file reads the DNS data and analyses the wall normal velocity component (u). ...
At the x = 0.06 location, FFT methods are used to understand the flow behaviour. The effects of filtering are also employed to understand the flucations and turbulence charactersitics ...
of the flow. Moving - average filtering is used here. Hamming filtering was tried but does not work completely yet. Regular and log-log plots are used.
Furthermore, statistical analysis is carried out for the wall normal velocity component using mean, standard deviation, skewness and kurtosis measures. The Probability Density Function (pdf) ...
obtained and plotted. Pearson correlation is implemented to see if there is any linear relationship between the **x** and **y** direction velocity components.
