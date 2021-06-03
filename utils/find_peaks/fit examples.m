load chrom

>> [FitResults,MeanFitError]=peakfit(chrom,12.1,0.98,1,5,5.2,1,[12.1 0.2])
FitResults =
            1       11.872       539.12      0.24051       137.55
MeanFitError =0.4553
       

>> [FitResults,MeanFitError]=peakfit(chrom,5.8,2.5,3,5,5.2,10,[5.16 0.17 5.8 0.17 6.44 0.17])
FitResults =
            1       5.2464        21.46      0.24923       5.6937
            2       5.7185       35.901      0.29594        11.31
            3       6.2041       32.948      0.23499       8.2409
MeanFitError =5.0218
       

load mike 

Peak Shape = Gaussian
Autozero ON
Number of peaks = 3
Fitted x range = 4.695 - 5.675 (dx=0.98)  (5.185)  
Percent Error = 1.707
         Peak#   Position       Height       Width         Area
            1       5.0608       116.74      0.37559       46.164
            2       5.0565       164.05      0.12052       21.047
            3       5.3923       245.24      0.16777       43.796


[FitResults,MeanFitError]=peakfit(mike,5.185,0.98,4,1,46.0192,1,[4.9024 0.21058 5.0592 0.14844 5.3169 0.26528 5.4033 0.14041]) 
Peak Shape = Gaussian
Autozero ON
Number of peaks = 4
Fitted x range = 4.695 - 5.675 (dx=0.98)  (5.185)  
Percent Error = 1.1325
         Peak#   Position       Height       Width         Area
            1       4.9024       60.076      0.21058       13.329
            2       5.0592       257.25      0.14844       40.649
            3       5.3169       107.06      0.26528       30.213
            4       5.4033       178.31      0.14041       26.652

            
[FitResults,MeanFitError]=peakfit(datamatrix,5.135,0.98,4,1,46.0192,1,[4.8946 0.24188 5.0627 0.1566 5.247 0.13053 5.3935 0.15815]) 
Peak Shape = Gaussian
Autozero ON
Number of peaks = 4
Fitted x range = 4.645 - 5.625 (dx=0.98)  (5.135)  
Percent Error = 1.0516
         Peak#   Position       Height       Width         Area
            1       4.8946       63.378      0.24188       16.196
            2       5.0627       263.72       0.1566       43.965
            3        5.247       71.331      0.13053       9.9119
            4       5.3935       258.34      0.15815       43.481
            

[FitResults,MeanFitError]=peakfit(mike,5.17,1.03,5,1,1.6639,1,[4.8118 0.12767 4.9157 0.1 5.0578 0.16 5.2508 0.14 5.3957 0.16])
FitResults =
            1        4.802       41.302      0.11668       5.1223
            2       4.9181       59.565      0.12024       7.6246
            3       5.0596       281.85      0.15515       46.553
            4       5.2412       70.516      0.13501       10.135
            5       5.3933       260.47      0.16225       44.987
MeanFitError =0.682883


Peaks near 4
Fitting Error 0.92384
          Peak#     Position     Height      Width         Area  
            1       3.8137       83.295      0.17501       15.508
            2       4.0559       375.66      0.20587        82.33
            3         4.29       113.52      0.18763       22.665
[FitResults,MeanFitError]=peakfit(datamatrix,4.065,0.98,3,1,46.0192,1,[3.8137     0.17501      4.0559     0.20587        4.29     0.18763]) 

Peak Shape = Gaussian
Autozero ON
Number of peaks = 4
Fitted x range = 3.575 - 4.555 (dx=0.98)  (4.065)  
Percent Error = 0.89381
         Peak#   Position       Height       Width         Area
            1       3.8161       80.107      0.19006       16.184
            2       4.2484       56.705      0.09318       5.6248
            3       4.0604       377.15      0.21092       84.681
            4       4.3269       88.417      0.15011       14.126
[FitResults,MeanFitError]=peakfit(mike,4.065,0.98,4,1,3.8169,1,[3.8161     0.19006      4.2484     0.09318      4.0604     0.21092      4.3269     0.15011]) 
[FitResults,MeanFitError]=peakfit(mike,4.065,0.98,4,1,3.8169,1,[3.5772    0.021602      3.8137     0.17498      4.0559     0.20588        4.29     0.18762]) 

[FitResults,MeanFitError]=peakfit(mike,4.065,0.98,4,1,3.8169,1,[3.7768     0.10019      4.0625     0.16176      4.0327      0.3531      4.2977     0.18318]) 
FitResults =
            1       3.7649       44.125      0.10711       5.0313
            2       4.0578       373.76      0.19355       77.009
            3       3.8644       71.654      0.17769       13.553
            4       4.2825       117.72      0.19748       24.734
MeanFitError =
       0.5827

Peaks near 7.4
Peak Shape = Gaussian
Autozero ON
Number of peaks = 3
Fitted x range = 6.905 - 7.855 (dx=0.95)  (7.38)  
Percent Error = 1.2684
         Peak#   Position       Height       Width         Area
            1       7.2244       303.51      0.20017       64.667
            2       7.4947       90.178      0.16864       16.189
            3       7.6689       100.18      0.12877       13.728
            
            Peak Shape = Gaussian
Autozero ON
Number of peaks = 4
Fitted x range = 6.905 - 7.855 (dx=0.95)  (7.38)  
Percent Error = 0.62434
         Peak#   Position       Height       Width         Area
            1       7.0923       32.912      0.10626       3.7229
            2       7.4909       90.244      0.19504       18.737
            3       7.2279       307.64      0.18154       59.453
            4        7.672       95.881      0.12549       12.80
[FitResults,MeanFitError]=peakfit(datamatrix,7.38,0.95,4,1,3.8169,1,[7.0923     0.10626      7.4909     0.19504      7.2279     0.18154       7.672     0.12549]) 

[FitResults,MeanFitError]=peakfit(mike,7.38,0.95,4,1,3.8169,1,[7.0923     0.10626      7.4909     0.19504      7.2279     0.18154       7.672     0.12549]) 
FitResults =
            1        7.099       40.628      0.11559       4.9991
            2       7.4923       89.661      0.21323       20.352
            3       7.2289       306.32      0.17626       57.477
            4       7.6746       92.155      0.12074       11.842
MeanFitError =
      0.59875