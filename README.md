# Reproducible research: version control and R

\# INSERT ANSWERS HERE #
#Question 4

```{r}
#This code below is copied from the question
#install.packages("ggplot2")
#install.packages("gridExtra")

library(ggplot2)
library(gridExtra)

random_walk  <- function (n_steps) {
  
  df <- data.frame(x = rep(NA, n_steps), y = rep(NA, n_steps), time = 1:n_steps)
  
  df[1,] <- c(0,0,1)
  
  for (i in 2:n_steps) {
    
    h <- 0.25
    
    angle <- runif(1, min = 0, max = 2*pi)
    
    df[i,1] <- df[i-1,1] + cos(angle)*h
    
    df[i,2] <- df[i-1,2] + sin(angle)*h
    
    df[i,3] <- i
    
  }
  
  return(df)
  
}

data1 <- random_walk(500)

plot1 <- ggplot(aes(x = x, y = y), data = data1) +
  
  geom_path(aes(colour = time)) +
  
  theme_bw() +
  
  xlab("x-coordinate") +
  
  ylab("y-coordinate")

data2 <- random_walk(500)

plot2 <- ggplot(aes(x = x, y = y), data = data2) +
  
  geom_path(aes(colour = time)) +
  
  theme_bw() +
  
  xlab("x-coordinate") +
  
  ylab("y-coordinate")

grid.arrange(plot1, plot2, ncol=2)
```
- This code has produced two figures of a simulation of a random walk. 
- The X axis of the figures is the X-coordinate, and the Y-axis is the Y-coordinates.
- For both of the simulations, they start at coordinates (0,0). At this point in time the line which shows the walk is a dark blue colour, as it transitions through time, it becomes lighter.
- The walk is simulated for 500 steps. For each of these steps, the simulation moves a set distance of 0.25 at an angle which is randomly generated. It works out the change in the X and Y coordinates by taking the sine and cosine of the angle, and plots this 500 times.
- The two simulations have different coordinates on the axes due to them moving in different directions randomly.
- The one on the left moves further down and to the left, whereas the one on the right stays more central.

- A random seed is a starting point for generating a reproducible sequence of random numbers. So, each time the code is run, you will get the same numbers which have been randomly generated.
- The random number generator function requires a seed to the run the algorithm, so if you know the seed and the generator used, you can predict and reproduce the output of the random number generator. 
- If the seed is not specified, R instead uses the clock of the system to establish the seed.
- It is particularly useful for reproducing the same output in simulations.

```{r}
#Copied the code below from the question but made it be able to reproduce the same outputs every time, by setting a seed for the random numbers.
random_walk  <- function (n_steps) {
  
  df <- data.frame(x = rep(NA, n_steps), y = rep(NA, n_steps), time = 1:n_steps)
  
  df[1,] <- c(0,0,1)
#This here is where I set the seed to 5 (Doesn't have to be this specific number.), and this ensures the same output is reproducible.
  set.seed(5)
  
  for (i in 2:n_steps) {
    
    h <- 0.25
    
    angle <- runif(1, min = 0, max = 2*pi)
    
    df[i,1] <- df[i-1,1] + cos(angle)*h
    
    df[i,2] <- df[i-1,2] + sin(angle)*h
    
    df[i,3] <- i
    
  }
  
  return(df)
  
}

data1 <- random_walk(500)

plot1 <- ggplot(aes(x = x, y = y), data = data1) +
  
  geom_path(aes(colour = time)) +
  
  theme_bw() +
  
  xlab("x-coordinate") +
  
  ylab("y-coordinate")

data2 <- random_walk(500)

plot2 <- ggplot(aes(x = x, y = y), data = data2) +
  
  geom_path(aes(colour = time)) +
  
  theme_bw() +
  
  xlab("x-coordinate") +
  
  ylab("y-coordinate")

grid.arrange(plot1, plot2, ncol=2)
```
#Question 5

```{r}
#naming the data frame from cui_etal2014 virusData
virusData <- read.csv("Cui_etal2014.csv")
virusData
#Counting the rows and columns
rows <- nrow(virusData)
columns <- ncol(virusData)
#Printing the number of rows
rows
#Printing the number of columns
columns
```
- There are 33 rows and 13 columns in the Cui_etal2014.csv file.
- Taking the natural logarithm of both the virion volume and the genome length would allow a linear model to be fitted to the transformed data.

```{r}
#Making a new column in virusData called logVolume with the natural logarithm of the virion volumes.
virusData$logVolume <- log(virusData$Virion.volume..nm.nm.nm.)
#Making a new column in virusData called logGenome, with the natural logarithm of the genome lengths.
virusData$logGenome <- log(virusData$Genome.length..kb.)
#Plotting a scatter graph with logGenome on the x-axis and logVolume on the y-axis.
ggplot(virusData, aes(x = logGenome, y = logVolume))+
  labs(x = "Log[Genome Length (kb)]", y = "Log[Virion Volume mm^3]", title = "Transformation to Fit a Linear Model")+
  geom_point()
```
```{r}
#LogVolume = alpha * logGenomeLength + logBeta
#Doing a linear regression to find alpha and beta
linearModel <- lm(logVolume ~ logGenome, data = virusData)
summary(linearModel)
```
- The summary of the linear model shows that alpha is 1.5152, as it is the coefficient of logGenome. The P value for alpha is 6.44e-10, so is statistically significant.
- The intercept is logBeta.
- logBeta = 7.0748
- beta = e^7.0748 
- beta = 1181.81
- The P value for beta is 2.28e-10 so is also statistically significant. 
- In table 2 of the paper it shows that for dsDNA viruses, alpha = 1.52, the alpha I worked out when rounded to 3 significant    figures is also 1.52. The value it showed for the scaling factor beta was 1182, and I also found it to be 1182 when rounded. 
- So, it can be concluded that the linear model here produced the same values as the paper and is statistically significant. 

```{r}
ggplot(virusData, aes(x = logGenome, y = logVolume))+
  labs(x = "Log[Genome Length (kb)]", y = "Log[Virion Volume mm^3]")+
#This geom_smooth function adds a linear regression line in blue, and the se = TRUE adds the shading around the lines for the confidence intervals.
  geom_smooth(method = "lm", se = TRUE, color = "blue")+
  geom_point()
```
- Estimated volume of a 300 kb dsDNA virus
- V = Beta*Genomelength^alpha
- V = 1181.81 * 300^1.5152
- V = 6697079


## Instructions

The homework for this Computer skills practical is divided into 5 questions for a total of 100 points (plus an optional bonus question worth 10 extra points). First, fork this repo and make sure your fork is made **Public** for marking. Answers should be added to the # INSERT ANSWERS HERE # section above in the **README.md** file of your forked repository.

Questions 1, 2 and 3 should be answered in the **README.md** file of the `logistic_growth` repo that you forked during the practical. To answer those questions here, simply include a link to your logistic_growth repo.

**Submission**: Please submit a single **PDF** file with your candidate number (and no other identifying information), and a link to your fork of the `reproducible-research_homework` repo with the completed answers. All answers should be on the `main` branch.

## Assignment questions 

1) (**10 points**) Annotate the **README.md** file in your `logistic_growth` repo with more detailed information about the analysis. Add a section on the results and include the estimates for $N_0$, $r$ and $K$ (mention which *.csv file you used).
   
2) (**10 points**) Use your estimates of $N_0$ and $r$ to calculate the population size at $t$ = 4980 min, assuming that the population grows exponentially. How does it compare to the population size predicted under logistic growth? 

3) (**20 points**) Add an R script to your repository that makes a graph comparing the exponential and logistic growth curves (using the same parameter estimates you found). Upload this graph to your repo and include it in the **README.md** file so it can be viewed in the repo homepage.
   
4) (**30 points**) Sometimes we are interested in modelling a process that involves randomness. A good example is Brownian motion. We will explore how to simulate a random process in a way that it is reproducible:

   - A script for simulating a random_walk is provided in the `question-4-code` folder of this repo. Execute the code to produce the paths of two random walks. What do you observe? (10 points)
   - Investigate the term **random seeds**. What is a random seed and how does it work? (5 points)
   - Edit the script to make a reproducible simulation of Brownian motion. Commit the file and push it to your forked `reproducible-research_homework` repo. (10 points)
   - Go to your commit history and click on the latest commit. Show the edit you made to the code in the comparison view (add this image to the **README.md** of the fork). (5 points)

5) (**30 points**) In 2014, Cui, Schlub and Holmes published an article in the *Journal of Virology* (doi: https://doi.org/10.1128/jvi.00362-14) showing that the size of viral particles, more specifically their volume, could be predicted from their genome size (length). They found that this relationship can be modelled using an allometric equation of the form **$`V = \beta L^{\alpha}`$**, where $`V`$ is the virion volume in nm<sup>3</sup> and $`L`$ is the genome length in nucleotides.

   - Import the data for double-stranded DNA (dsDNA) viruses taken from the Supplementary Materials of the original paper into Posit Cloud (the csv file is in the `question-5-data` folder). How many rows and columns does the table have? (3 points)
   - What transformation can you use to fit a linear model to the data? Apply the transformation. (3 points)
   - Find the exponent ($\alpha$) and scaling factor ($\beta$) of the allometric law for dsDNA viruses and write the p-values from the model you obtained, are they statistically significant? Compare the values you found to those shown in **Table 2** of the paper, did you find the same values? (10 points)
   - Write the code to reproduce the figure shown below. (10 points)

  <p align="center">
     <img src="https://github.com/josegabrielnb/reproducible-research_homework/blob/main/question-5-data/allometric_scaling.png" width="600" height="500">
  </p>

  - What is the estimated volume of a 300 kb dsDNA virus? (4 points)

**Bonus** (**10 points**) Explain the difference between reproducibility and replicability in scientific research. How can git and GitHub be used to enhance the reproducibility and replicability of your work? what limitations do they have? (e.g. check the platform [protocols.io](https://www.protocols.io/)).
