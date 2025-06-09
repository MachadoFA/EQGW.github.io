
|        |        |        |    |
|--------|---------------------------------------------|--------------------|------------------------------------------|
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [2025 Workshop](/index.html) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [SCHEDULE](/2025/schedule.html) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [RESOURCES](/2025/resources.html) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [PREVIOUS YEARS](2025/previous.html) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; |


<div align="left">
<img src="/media/SWVirginiaMtns.jpg" alt="[Southwest Virginia Mountains]">
</div>


<table><tr><td>&larr; <a href="/2025/lecture1-3/lecture1-3.html">Previous</a></td><td width="772">&nbsp;</td><td> <a href="/2025/lecture2-1/lecture2-1.html">Next &rarr;</a></td></tr></table>
[//]: # (This is a comment. Edit the Next and Previous links above to go the right links)  

## Exercise 1-1: Selection & Inheritance Lab Exercise ##

### Instructors: Pat Carter, Steve Arnold, & Josef Uyeda ###
  
In this computer exercise we will improve our intuition about heritability in a couple of
ways. First, we will use actual data from a natural population to obtain a point estimate
of heritability using the statistical procedure of regression. We will examine data on
garter snakes and use regression to estimate heritability for 6 traits in our dataset. Second,
we will use bootstrapping to see how our point estimate of heritability (2 times the
regression slope) is affected by the particular sample of mothers and daughters that we
use.

#### Background readings:  ####
You can read more about these data and analyses in these papers from which the data
were originally presented.
* [Arnold & Phillips, 1999](/2025/exercise1-1/Arnold-Phillips-1999.pdf)
* [Phillips & Arnold, 1999](/2025/exercise1-1/Phillips-Arnold-1999.pdf)

#### Data ####

* [Dataset](/2025/exercise1-1/R_inland_snake_data.txt)
Save these data onto your computer. The data are the values of 6 traits in a 
sample of mothers and their daughters from one population of a garter snake, 
Thamnophis elegans. The first line is the header which lists the names of the 
variables (8 columns):
* local = the name of the population
* family = the name of the family (the first row with each family name is the
mother, the rows that follow with the same name are her daughters)
* body = the number of vertebrae in the body
* tail = the number of vertebrae in the tail
* mid = the number of scale rows around the snake, assessed at the middle point of
the body
* ilab = the number of scales on the lower lip (infralabials, sum of left and right
sides)
* slab = the number of scales on the upper lip (supralabials, sum of left and right
sides)
* post = the nmber of scales behind the eye (postoculars, sum of left and right sides)

#### Code ####
We will use an R markdown file. 
* [Heritability R notebook](/2025/exercise1-1/Exercise-1_1-Heritability-of-snake-vertebral-numbers_vers11.nb.Rmd)
* [Notebook w/answers R notebook](/2025/exercise1-1/Exercise-1_1-Heritability-of-snake-vertebral-numbers_vers11_answers.nb.Rmd)
* [Notebook w/answers web page](/2025/exercise1-1/Exercise-1_1-Heritability-of-snake-vertebral-numbers_vers11_answers.nb.html)

