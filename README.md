# SVAR-MSH-ID
### R Codes for Bayesian Inference for Structural Vector Autoregressions Identified with Markov-Switching Heteroskedasticity

by Tomasz Woźniak

> **Summary.** The project provides source codes and reproduction files for the paper mentioned below. Utility functions for the Gibbs sampler for Bayesian Structural Vector Autoregressions with Markov-Switching Heteroskedasticity, Savage-Dickey density ratio for uniqueness conditions and homoskedasticity hypothesis, and a marginal data density estimator are provided.
>
> **Keywords.** SVAR-MSH, identification through heteroskedasticity, Savage-Dickey density ratio for uniquenss conditions, Gibbs sampler

## Citation

To refer to the codes in publications, please, cite the following paper:

> Lütkepohl, H., Woźniak, T. (2020) Bayesian Inference for Structural Vector Autoregressions Identified with Markov-Switching Heteroskedasticity, *Journal of Economic Dynamics and Control*, 113, DOI: [10.1016/j.jedc.2020.103862](https://doi.org/10.1016/j.jedc.2020.103862).

## Corrigendum

A correct version of the first equation on the top of page 17 from the Appendix of the paper is given by:
![](corrigendum.png)

## Project contents

The project includes:

- folder `codes` with the whole source code
- R files for the reproduction of most of the results from the paper
- file `dataBI2015.RData` with data used in the paper

## Downloading the codes

To download the codes follow the steps (requires a [GitHub account](https://github.com) and [GitHub Desktop](https://desktop.github.com/)):

1. **On this page:** fork the project by clicking the icon on the top of this page ([more info](https://guides.github.com/activities/forking/))

   ![](https://github-images.s3.amazonaws.com/help/bootcamp/Bootcamp-Fork.png)

2. **In GitHub Desktop:** go to the application menu to select: `File -> Clone Repository...`, choose the forked repository from the list of available repositories from GitHub.com, select your local folder to store the repository, and click `Clone`. ([even more info](https://docs.github.com/en/get-started/quickstart/fork-a-repo))

3. **Optional:** if you find a bug or if you improve the code, please feel free to submit a pull request. (more info [here](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests) and [here](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request))
