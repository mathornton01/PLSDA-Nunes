# PLSDA-Nunes

In order to install R-libraries directly from Github, the user must first install the devtools library with: 
```r
install.packages('devtools'); 
```
Once the devtools package has been installed, the user can load the library using the standard 'library' or 'require' functions
then from the devtools package the user can install the PLSDA-R-Package by using the 'install_github' function. 
```r
library(devtools);
install_github('mathornton01/PLSDA-Nunes'); 
```
Now R will display several messages that give various updates of the status of the download and installation.  Once the installation 
has been completed, the user can load the library with the typical 'library' and 'require' functions. 
```r
library(YinGenomicDFTDistances); 
```

## Checking the Documentation 

Users proficient in R already know that the man command or '?' prior to a function will allow the user to pull up a reference on 
a particular function in a package if that documentation exists.  To check the documentation for the functions in the Nunes PLSDA package, the same ? notation may be used.  For instance, to check the documentation on how to Use PLSDA on some data
one may use the following command: 

```r
?PLS
``` 
