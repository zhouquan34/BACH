_Copyright 12/04/2015 Quan Zhou, advised by Yongtao Guan, Baylor College of Medicine_

This program calculates the p-value of a weighted sum of chi-squared varibles with d.f. = 1.
The algorithm was first described by Bausch. 

According to your operating system, you can run bach-linux or bach-mac. 
If they don't work or you want to compile yourself, please refer to INSTALL for instructions.

In the following, we describe the basic usage of bach. 

##  Single Computation   
Suppose Y = X1 + 0.8X2 + 0.6X3 + 0.4X4 where each X~chi-sq-1, 
and you want to compute P(Y>10) and P(Y>20).

Then type
       ./bach 1 0.8 0.6 0.4 10,20

Two p-values will be printed to the stdout.

Rules of input:
	Coefficients must be separated by space. 
	The values of the test statisitc must be separated by comma. 
	Coefficients must be given before the values of the test statistic.
	At least one value of the test statistic must be given.

##   Batch Mode   
To calculate multiple density functions, you can specify the input (and output) files to invoke the batch mode. 

For example,
	./bach -i test_in.txt 

In the input file, each line specifies the coefficients and the values of the test statistic, according to the rules above.

For example, if the input file contains:
	1 0.9 10,20
	1 0.8 5

You will see three p-values printed to stdout:
	P(Q1 + 0.9Q2 > 10), P(Q1 + 0.9Q2 > 20)
	P(Q1 + 0.8Q2 > 5)

if '-o' argument is NOT provided. 
If the output file is specified, p-values and other information will be written to the output.
	
You can also find useful examples in 'test_data'. If you don't know how to run tests, please refer to file Run_Tests.

## Output format  
In batch mode, if the output file name is provided by '-o' argument, you will see an output file composed of 4 columns:

	p-value    error-bound    coefficients    statistic


## Argument list 
The type and default value of each argument are given in the square brackets.  

-i [string=]: the input filename (the only required argument for batch mode).
-o [string=]: the output filename.
-c [int=6]: precision (the number of significant digits) for output.
-e [int]: desired relative error bound of P. '-e 2' means error/P < 1e-2.

## References 
1. Johannes Bausch, On the efficient calculation of a linear combination of chi-squared random variables with an application in counting string vacua. _Journal of Physics A: Mathematical and Theoretical_, 46(50):505202, 2013. 
2. Quan Zhou and Yongtao Guan. On the null distribution of Bayes factors in linear regression. _Journal of the American Statistical Association_, 113(523):1362-1371, 2017.
