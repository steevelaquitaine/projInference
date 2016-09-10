## Tables

**Table 1**

| Subject | 01 | 02 | 03 |04 |05 |06 |07 |08 |09 |10 |11 |12 |
|---------|----|----|----|----|----|----|----|----|----|----|----|----|
| Bayesian obs| 77182.68 | 74175.23|84145.32|46451.09|61913.16|75357.9|54584.24|62731.75|84927.68|62633.63|63233.81|66099.47|
| Likelihood obs|100803.6|92739.7|110810.02|56504.8|68159.3|88925.5|68253.5| 68253.5|101627.7|71161.2|76223.3|76223.3|


**Table 2**

| Subject | ll1 | ll2 |ll3|p1|p2|p3|p4|c|r|m|t|
| ------- | --- | --- | --- | --- | --- | --- | --- | ---- | ---- | ---- | ----- |
| 01 | 80  | 40  |  20 | 1.74| 4.77|10.74|34.25| NaN |0.001| 15| NaN|     
| 02 | 101 |70| 22| 1| 3| 9.74| 34| NaN| 0.0007| 6.9| NaN| 
|03 |80| 40| 16| 1.7| 4.8| 11| 33| NaN| 0.00039| 17| NaN| 
|04| 79.3| 38.9| 9.19| 1.7| 4.8| 11| 33| NaN| 0.00039| 17| NaN| 
|05| 80| 5| 1| 0| 0| 0| 0| NaN| 0.001| 13.5| NaN|
|06| 30| 3| 1| 0.01| 1.001| 1.7| 2| NaN| 0.001| 15| NaN|
|07| 80| 35| 15| 1.74| 4.77| 10.74| 34.25| NaN| 0.001| 13.5| NaN| 
|08| 60| 10| 5| 1.535| 4.3425| 10.62| 40| NaN| 0.00062| 13.975| NaN| 
|09| 25| 3| 1.5| 0.01| 0.2| 2| 3.2| NaN| 0.001| 15| NaN|
|10| 200| 1| 1| 0.18| 1| 21.63| 300.63| NaN| 0.00104| 5| NaN|
|11| 25| 3| 1.5| 0.06| 0.2| 2.1| 3.4| NaN| 0.001| 15| NaN|
|12| 5| 2| 1| 0.01| 1.00015| 1.035| 3| NaN| 0.001| 15| NaN|

ll1/2/3 are circular von Mises likelihood k parameters for 24%,12%,6% coherence; p1/2/3/4 prior von Mises k parameters for 80, 40,20, 10 deg priors;c is the k parameter of each von Mises in the mixture of von Mises cardinal prior; r proportion of random estimation; m is the k parameter of the convolved von Mises motor noise; t is a tail parameter (mixture coefficient of von Mises prior and uniform) that controls prior kurtosis (tailedness).
Different sets of initial parameters were tested (set 1: above; set 2 : set 1 with 8*p; set 3 : set 1 with p/8; set 4 : set 1 with ll*8; set 5: set 1 with ll and p*8; set 6 : set 1 with l*8 and p/8; set 7 : set 1 with l/8; set 8 : set 1 with l/8 and p*8; set 9 : set 1 with l and p/8; set 10: set 1 with same averaged l and p).

**Table 3 : Learning Bayesian initial parameters**

Initial parameters were best fit parameters obtained from teh basic Bayesian observer

| Subject | ll1 | ll2 |ll3|p1|p2|p3|p4|c|r|m|t|tau|
| ------- | --- | --- | --- | --- | --- | --- | --- | ---- | ---- | ---- | ----- |----- |
| 01 |23.78| 7.56| 1.88| 0.66| 1.71| 2.77| 5.30| NaN| 0.001| 83.63| NaN| 20|
| 02 |16.2| 0.91| 0.75 |0.65 |1.28 |12 |32.1| NaN |0.0007| 11.7| NaN| 20 |
| 03 |34.74| 5.5| 1.9 |0.26 |1.56| 3.5 |8.9 |NaN |0.0003| 65.4 |NaN |20|
| 04 |12.6| 2.2| 1 |0.4 |0.9| 3| 10.4| NaN| 0.0006| 19.2| NaN| 20|
| 05 |98.8| 2.35| 1 |0.0004| 0.0005| 0.00098| 0.0005| NaN| 0.0008| 9.9 |NaN| 20|
| 06 |36.8| 2.8| 1| 3.9| 0.8| 2.16| 6.07| NaN| 0.007| 11.3| NaN| 20|
| 07 |20.28| 4.28| 1.55| 0.038| 0.78| 2.7| 13.4| NaN| 0.002| 36.75| NaN| 20|
| 08 |9.1| 1.7| 0.9| 0.6| 4.4| 13.6| 39.4| NaN| 0.0005| 2.2| NaN| 20|
| 09 |103| 3| 1.4| 0.76| 0.61| 2.6| 3.3| NaN| 0.004| 14.6| NaN| 20|
| 10 |5.7| 0.6| 0.6| 3.2^-5 |0.4 |2.2 |8.3 |NaN| 0.01| 9.2| NaN | 20|
| 11 |29.4| 3.1| 1.4| 0.02 |0.25 |1.3 |6| NaN| 0.001| 27| NaN| 20|
| 12 |3.4| 1.6| 1.1| 0.007| 0.7| 1| 4.6| NaN| 0.001| 22| NaN| 20|

**Table 4 : Bayesian learning rate tau from fit**

| Subject | 01 | 02 | 03 |04 |05 |06 |07 |08 |09 |10 |11 |12 |
|---------|----|----|----|----|----|----|----|----|----|----|----|----|
|Best tau|0.29| 0.08| 4.00| 18.06| 20.00| 20.35| 1.68| 20.00| 21.20| 55.11| 16.69| 1.57|
