# Bessel-Hankel

**Bessel function** and **Hankel transform**

[中文版](README.cn.md)

## License

[The MIT License](http://tchel.mit-license.org/)

## Author

Tche LIU, seistche@gmail.com, USTC

## Methodology

The repository includes some examples for Hankel transform, and involves calculations of Bessel functions of the 1st kind and 2nd kind, their derivatives of the 1st order and 2nd order, and their zeros. Next, I will show you some mathematical formulas about calculations of Bessel function's derivatives and zeros, and Hankel transform.

### Bessel function's derivatives

Represent $J_\nu(x)$ or $Y_\nu(x)$ by $Z_\nu(x)$ in the follows, and there are recurrence relations of Bessel function:

$$ x Z_{\nu - 1}(x) + x Z_{\nu + 1}(x) = 2\nu Z_\nu(x) $$

$$ Z_\nu'(x) = \frac{\nu}{x} Z_\nu(x) - Z_{\nu + 1}(x) $$

so that

$$ Z_{\nu + 2}(x) = \frac{2(\nu + 1)}{x} Z_{\nu + 1}(x) - Z_\nu(x) $$

$$ Z_{\nu + 1}'(x) = \frac{\nu + 1}{x} Z_{\nu + 1}(x) - Z_{\nu + 2}(x) $$

And it's not difficult to reduce the 2nd derivative:

$$ \begin{align*} 
     Z_\nu''(x) & = \frac{d}{dx} \left[ \frac{\nu}{x} Z_\nu(x) - Z_{\nu + 1}(x) \right] \newline
                & = - \frac{\nu}{x^2} Z_\nu(x) + \frac{\nu}{x} Z_\nu'(x) - Z_{\nu + 1}'(x) \newline
                & = - \frac{\nu}{x^2} Z_\nu(x) + \frac{\nu}{x} \left[ \frac{\nu}{x} Z_\nu(x) - Z_{\nu + 1}(x) \right] - \left[ \frac{\nu + 1}{x} Z_{\nu + 1}(x) - \frac{2(\nu + 1)}{x} Z_{\nu + 1}(x) + Z_\nu(x) \right] \newline
                & = \left( - \frac{\nu}{x^2} + \frac{\nu^2}{x^2} - 1 \right) Z_\nu(x) + \left( - \frac{\nu}{x} - \frac{\nu + 1}{x} + \frac{2\nu + 2}{x} \right) Z_{\nu + 1}(x) \newline
                & = \frac{1}{x^2} \left[ (\nu^2 - \nu - x^2) Z_\nu(x) + x Z_{\nu + 1}(x) \right]
   \end{align*} ​$$

and

$$ \begin{align*} 
     Z_\nu''(x) & = \left( \frac{\nu^2}{x^2} - 1 \right) Z_\nu(x) - \frac{1}{x} \left[ \frac{\nu}{x} Z_\nu(x) - Z_{\nu + 1}(x) \right] \newline
                & = \left( \frac{\nu^2}{x^2} - 1 \right) Z_\nu(x) - \frac{1}{x} Z_\nu'(x)
   \end{align*} ​$$

### Bessel function's zeros

#### Fixed point iteration

According to [Wikipedia](https://en.wikipedia.org/wiki/Fixed-point_iteration#Applications), the fixed point iteration to solve zeros of Bessel function is:

$$ \begin{align*}
     x_{n + 1} & = x_n - \frac{Z_\nu(x_n)}{Z_\nu'(x_n)} \newline
               & = x_n - \frac{x_n Z_\nu(x_n)}{\nu Z_\nu(x_n) - x_n Z_{\nu + 1}(x_n)}
   \end{align*} ​$$

#### Halley's method

Based on the above formulas, these results of product are:

$$ Z_\nu(x) Z_\nu'(x) = Z_\nu(x) \left[ \frac{\nu}{x} Z_\nu(x) - Z_{\nu + 1}(x) \right] = \frac{1}{x^2} x Z_\nu(x) \left[ \nu Z_\nu(x) - x Z_{\nu + 1}(x) \right] $$

$$ [Z_\nu'(x)]^2 = \frac{\nu^2}{x^2} Z_\nu(x) - \frac{2\nu}{x} Z_\nu(x) Z_{\nu + 1}(x) + [Z_{\nu + 1}(x)]^2 = \frac{1}{x^2} \left[ \nu^2 Z_\nu(x) - 2\nu x Z_\nu(x) Z_{\nu + 1}(x) + x^2 [Z_{\nu + 1}(x)]^2 \right] $$

$$ Z_\nu(x) Z_\nu''(x) = \frac{1}{x^2} \left\lbrace (\nu^2 - \nu - x^2) [Z_\nu(x)]^2 + x Z_\nu(x) Z_{\nu + 1}(x) \right\rbrace $$

According to [Wikipedia](https://en.wikipedia.org/wiki/Halley's_method#Method), the iteration of Halley's method to solve zeros of Bessel function is:

$$ \begin{align*}
     x_{n + 1} & = x_n - \frac{2 Z_\nu(x_n) Z_\nu'(x_n) }{ 2 [Z_\nu'(x_n)]^2 - Z_\nu(x_n) Z_\nu''(x_n) } \newline
               & = x_n - \frac{ 2 x_n Z_\nu(x_n) [ \nu Z_\nu(x_n) - x_n Z_{\nu + 1}(x_n) ] }{ 2 x_n^2 [Z_{\nu + 1}(x_n)]^2 - (4\nu + 1) x_n Z_\nu(x_n) Z_{\nu + 1}(x_n) + (\nu^2 + \nu + x_n^2) [Z_\nu(x_n)]^2}
   \end{align*} ​$$

### Hankel transform

#### Definition

According to [Wikipedia](https://en.wikipedia.org/wiki/Hankel_transform#Definition), the Hankel transform of order $\nu$ of a function $f(r)$ is given by:

$$ F_\nu(k) = \int_0^\infty f(r) J_\nu(k r) r dr $$

#### (Guptasarma and Singh, 1997)

The type of Hankel transform of this paper is (equations 2 and 3):

$$ f(r) = \int_0^\infty K(\lambda) J_i(r \lambda) d\lambda \approx \frac{1}{r} \sum_{i = 1}^n K(\lambda_i) W_i $$

The formula is **only** applicable for Hankel transform based on Bessel function of order 1 and 2.

There are some examples of standard Hankel transform to verify our programs in this paper:

|                         $K(\lambda)$                         | $J_\nu$ |                            $f(r)$                            |
| :----------------------------------------------------------: | :-----: | :----------------------------------------------------------: |
|                      $e^{ - c \lambda}$                      |  $J_0$  |                 $\frac{1}{\sqrt{c^2 + r^2}}$                 |
|                 $\lambda e^{ - c \lambda^2}$                 |  $J_0$  |                $\frac{1}{2c} e^{ - r^2/(4c)}$                |
|                  $\lambda e^{ - c \lambda}$                  |  $J_0$  |               $\frac{c}{\sqrt{(c^2 + r^2)^3}}$               |
| $\lambda e^{ - c \lambda} + \alpha \lambda^2 e^{ - c \lambda^2}$ |  $J_1$  | $\frac{r}{\sqrt{(c^2 + r^2)^3}} + \alpha \frac{r e^{ - r^2/(4c)}}{4 c^2}$ |
|                  $\lambda e^{ - c \lambda}$                  |  $J_1$  |               $\frac{r}{\sqrt{(c^2 + r^2)^3}}$               |
|                $\lambda^2 e^{ - c \lambda^2}$                |  $J_1$  |              $\frac{r e^{ - r^2/(4c)}}{4 c^2}$               |
|                      $e^{ - c \lambda}$                      |  $J_1$  |      $\frac{\sqrt{c^2 + r^2} - c}{r \sqrt{c^2 + r^2}}$       |

#### (Ogata, 2005)

The approximation of integrals of the Hankel transformation type of this paper is (equation 5.2):

$$ \begin{align*}
  \int_0^\infty f(x) J_\nu(x) dx & \approx \pi \sum_{k = 1}^{\infty} \omega_{\nu k} f(x) J_\nu(x) \psi'(t) \newline
                                 & = \pi \sum_{k = 1}^{\infty} \omega_{\nu k} f \left( \frac{\pi}{h} \psi(h \xi_{\nu k}) \right) J_\nu \left( \frac{\pi}{h} \psi(h \xi_{\nu k}) \right) \psi'(h \xi_{\nu k}) \newline
\end{align*}$$

where $ x = \pi / h \cdot \psi(t) $, $ \psi(t) = t \tanh(\pi / 2 \cdot \sinh t) $ and $ t = h \xi_{\nu k} $.

Set $ x = r \lambda ​$, so that $ \lambda = x/r ​$ and $ d\lambda = 1/r \cdot dx ​$, and $ x = \pi/h \cdot \psi(h \xi_{\nu k}) ​$, we can reduce:

$$ \begin{align*}
     \int_0^\infty f(\lambda) J_\nu(r \lambda) \lambda d\lambda & = \int_0^\infty f \left( \frac{x}{r} \right) J_\nu(x) \frac{x}{r} \cdot \frac{1}{r} dx \newline
                                                                & = \frac{1}{r^2} \int_0^\infty f \left( \frac{x}{r} \right) J_\nu(x) x dx \newline
                                                                & \approx \frac{1}{r^2} \left[ \pi \sum_{k = 1}^{\infty} \omega_{\nu k} f \left( \frac{x}{r} \right) J_\nu(x) \psi'(h \xi_{\nu k}) \right] x
   \end{align*} ​$$

The formula is applicable for Hankel transform based on Bessel function of an **arbitrary** order.

## Implementation

### Bessel_Function

[Bessel_Function](Bessel_Function.F90) is a `fortran` module for calculations of Bessel functions of the 1st and 2nd kind, their derivative of the 1st and 2nd order, and zeros of original functions and derivations of the 1st order. The module is mainly modified from [mjyzo.f90](http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/mjyzo_f90.txt).

There are two macros for preprocessing in this source file: 
- **MOREAU** for using Moreau's method to calculate Bessel functions of the 1st and 2nd kind and their derivatives of the 1st order; and if not defined, using `fortran` intrinsic functions to calculate.
- **HALLEY** for using Halley's method to solve zeros of Bessel function of the 1st and 2nd kind; and if not defined, using the fixed point iteration to solve.

Because the source file has no main program, as an example, you can compile it by `$ gfortran -DHALLEY -c Bessel_Function.F90`, but CAN'T run it.

Besides, [Bessel Zero Solver](https://ww2.mathworks.cn/matlabcentral/fileexchange/48403-bessel-zero-solver) is an efficient and reliable `matlab` function to solve zeros of Bessel function of the 1st and 2nd kind based on Halley's method.

### Guptasarma_1997

[Guptasarma_1997](Guptasarma_1997.F90) is a independent program to implement the methodology of (Guptasarma and Singh, 1997), and includes all examples in this paper.

There is a macro for preprocessing in this source file: **FAST** for using less sample points to faster finish the Hankel transform; and if not defined, using more sample points to more exactly finish.

For example, you can compile the source file by `$ gfortran -DFAST Guptasarma_1997.F90 -o Guptasarma`, and run the program by `$ ./Guptasarma`.

### Ogata_2005

[Ogata_2005](Ogata_2005.F90) is a program based on [the Bessel_Function module](Bessel_Function.F90) to apply the methodology of (Ogata, 2005) to Hankel transform.

For example, you can compile the source file by `$ gfortran -DHALLEY Bessel_Function.F90 Ogata_2005.F90 -o Ogata`, and run the program by `$ ./Ogata`.

## References

- Guptasarma and Singh, 1997. **New digital linear filters for Hankel J_0 and J_1 transforms**. Geophysical Prospecting, 45, 745-762. [https://doi.org/10.1046/j.1365-2478.1997.500292.x](https://onlinelibrary.wiley.com/doi/abs/10.1046/j.1365-2478.1997.500292.x)
- Ogata, 2005. **A numerical integration formula based on the Bessel functions**. Publ. RIMS, 41, 949-970. [https://doi.org/10.2977/prims/1145474602](https://www.ems-ph.org/journals/show_abstract.php?issn=0034-5318&vol=41&iss=4&rank=8)

