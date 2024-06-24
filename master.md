## I. for Bern.
- Conj prior $\rightarrow Beta(\alpha, \beta)$
	$\alpha' = \alpha + \sum x_{i}$    ,    $\beta' = \beta + n - \sum x_{i}$
	$mean = \frac{\alpha + r}{\alpha + \beta + n}$    , $var = \frac{\alpha \beta (\alpha + \beta + n)}{(\alpha + \beta)^2 (\alpha + \beta + 1)}$   (need to substitute $\alpha'$ and $\beta'$)
	
	If $m$ and $s$ of the Beta are given
	$\alpha = \bigl(\frac{1-m}{s^2} - \frac{1}{m}\bigr)m^2$    ,     $\beta = \alpha \bigl(\frac{1}{m} - 1\bigr)$  (for a Binom. distr)

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

