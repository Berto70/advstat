## I. for Bern.
- Binom: Conj prior $\rightarrow Beta(\alpha, \beta)$
	$\alpha' = \alpha + \sum x_{i}$    ,    $\beta' = \beta + n - \sum x_{i}$
	$mean = \frac{\alpha + r}{\alpha + \beta + n}$    , $var = \frac{\alpha \beta (\alpha + \beta + n)}{(\alpha + \beta)^2 (\alpha + \beta + 1)}$   (need to substitute $\alpha'$ and $\beta'$)
	
	If $m$ and $s$ of the Beta are given
	$\alpha = \bigl(\frac{1-m}{s^2} - \frac{1}{m}\bigr)m^2$    ,     $\beta = \alpha \bigl(\frac{1}{m} - 1\bigr)$
- Geom: Beta prior $\Rightarrow$ Beta post w/    $\alpha'=\alpha+n \quad,\quad \beta'=\beta+\sum(k_{i}-1)$

## I. for Pois.
The likelihood is $\propto$ to a $Gamma(y, \alpha, \lambda)$ with $\alpha = \sum y_{i} + 1$    ,    $\lambda = n$
#### Priors
- Uniform $\Rightarrow$ post is Gamma w/    $\alpha = \sum y_{i} + 1 \quad,\quad\lambda = n$
- Jeffrey $\Rightarrow$ post is Gamma w/    $\alpha = \sum y_{i} + \frac{1}{2} \quad,\quad \lambda = n$
- **Gamma** $\rightarrow$ conj prior $\rightarrow$ $\alpha = r + 1 \quad,\quad \lambda=n$
	- single obs: $\Rightarrow$ post is Gamma w/    $\alpha' = \alpha + y \quad,\quad \lambda'=\lambda+1$
	- n obs {$y_{i}$}: $\Rightarrow$ post is Gamma w/    $\alpha' = \alpha + \sum y_{i} \quad,\quad \lambda'=\lambda+n$
	$$\begin{equation}
	E[n|y] = \frac{\alpha'}{\lambda'} \quad,\quad var=\frac{\alpha'}{\lambda'^2}
\end{equation}$$
	If $m$ and $s$ of the Gamma are given
	$\alpha = \left( \frac{m}{s} \right)^2 \quad,\quad \lambda=\frac{m}{s^2}$
If likelihood is $Exp$ and we use $Gamma$ prior $\Rightarrow$ post is Gamma w/    $\alpha'=\alpha+n \quad,\quad \beta'=\beta+n\bar{x}$
## I. for Norm.
$\sigma^2$ is known, while $\mu$ is a parameter. Likelihood is $Norm(\mu, \sigma^2)$
#### Priors
- Uniform $\Rightarrow$ post is Norm w/    $\mu_{0} = \frac{1}{N} \sum y_{i}$ (true empirical mean)    ,    $s^2 = \frac{\sigma^2}{N}$
	  The inference will be $\mu = \mu_{0} \pm \frac{\sigma}{\sqrt{N}}$  
	- If the data has errors $\rightarrow$ $\mu_{0} = \frac{\sum y_{i}/\sigma^2_{i}}{\sum 1/\sigma^2_{i}} \quad , \quad s^2 = \left( \sum 1/\sigma^2_{i} \right)^{-1}$
- $Norm(m, s^2)$ $\Rightarrow$ post is Norm w/    
	- single obs: $m'=\frac{\sigma^2 m + s^2 y}{\sigma^2 + s^2} \quad, \quad s'^2 = \frac{\sigma^2 s^2}{\sigma^2 + s^2}$
	- n obs {$y_{i}$}: $m' = \frac{1/s^2}{n/\sigma^2 + 1/s^2}m + \frac{n/\sigma^2}{n/\sigma^2 + 1/s^2}\bar{y} \quad , \quad s'^2 = \frac{\sigma^2 s^2}{\sigma^2 + n s^2}$

$\mu$ is known, while $\sigma^2$ is a parameter.
$\rightarrow$ prior is the Inverse Gamma $\rightarrow$ $E[x]=\frac{\beta}{\alpha - 1} \quad,\quad var(x)=\frac{\beta^2}{(\alpha-1)^2(\alpha-2)}$
$\Rightarrow$ post is inverse Gamma w/    $\alpha' = \alpha + \frac{n}{2} \quad,\quad \beta'=\beta+\frac{\sum(x_{i} - \mu)^2}{2}$

## Distributions
- Binom: x successes in n trials. $E[x]=np \quad,\quad var=np(1-p)$
- Geom: x failures to get the 1st success. $E[x]=\frac{1}{p} \quad,\quad var=\frac{1-p}{p^2}$
- Multinomial: generalize Binom $\rightarrow$ outcome $A_{i}$ with prob $p_{i}$. $E[x_{i}] = np_{i} \quad,\quad var(x)=np_{i}(1-p_{i})$
- Pois: $E[x]=\lambda \quad,\quad var=\lambda$
- NegBinom: prob of obtaining r-th success in n trials. $E[x] = \frac{r}{p} \quad,\quad var(x)=\frac{r(1-p)}{p^2}$
- Exp: prob of the distance between events in a Pois process. $E[x]=\frac{1}{\lambda} \quad,\quad var=\frac{1}{\lambda^2}$. $Exp(\lambda)\sim Gamma(1, \lambda)$
- Erlang: prob of waiting time x for the n-th event to occur. $E[x]=\frac{n}{\lambda} \quad,\quad var(x)=\frac{n}{\lambda^2}$. $Eral(n, \alpha)\sim Gamma(n, \lambda)$ 
- Gamma: $E[x]=\frac{\alpha}{\beta} \quad,\quad var(x)=\frac{\alpha}{\beta^2}$
- Beta: $E[x]=\frac{\alpha}{\alpha+\beta} \quad,\quad var(x)=\frac{\alpha \beta}{(\alpha+\beta)^2(\alpha + \beta +1)}$. 
## Combinatory
- Unique pairs (no order): $\frac{n(n-1)}{2}$
- Unique ordering: $n!$
- Permutations (order <u>matter</u>):
	- yes rep: $n^r$ (seq of r obj from n)
	- no rep: $\frac{n!}{(n-r)!}$ (unique selection of r obj among n)
- Combinations (order doesn't matter):
	- yes rep: $\frac{(n+r+1)!}{r!(n-1)!}$ (n° of ways of choosing r obj from n, w/ replacement)
	- no rep: $\frac{n!}{r!(n-r)!}$ (n° of ways of choosing r obj from n, w/o regard of order)
## Ineq.
Very useful when we don't have enough infos about the distr of random variables, but for which we can calculate $E[x]$ and/or $var(x)$
- Markov: $X\ge0$ w/ $E[x]=\mu$ $\Rightarrow$ $P(X\ge k)\le \frac{\mu}{k}$
- Jensen: $E[x^2]\ge (E[x])^2$ since $var(c)\ge0$ $\rightarrow$ $X$ w/ \\
$E[x]=\mu$ and $g(x)$ is a convex func $\Rightarrow$ $g(E[x])\le E[g(x)]$
- Cheby:  $X\ge0$ w/ $E[x]=\mu$ and $var(x)=\sigma^2$ $\Rightarrow$ $P(|X-\mu|\ge k)\le \frac{\sigma^2}{k^2}$ \\
  if $k=r\sigma\quad\rightarrow\quad P(|X-\mu|\ge r\sigma)\le \frac{1}{r^2}$ 