# Bessel-Hankel

**Bessel 函数**和 **Hankel 变换**

[Switch to English](README.md)

## 开源协议

[MIT 协议](http://tchel.mit-license.org/)

## 联系作者

Tche LIU, seistche@gmail.com, USTC

## 写在前面

本文中包含有些许 MathJax 公式，由于 GitHub 官方不支持此类公式渲染，公式显示略显凌乱。幸于 GitHub 网友用爱发电，针对 Chrome 和 Firefox 浏览器，分别安装 [MathJax Plugin for GitHub](https://chrome.google.com/webstore/detail/mathjax-plugin-for-github/ioemnmodlmafdkllaclgeombjnmnbima)（[GitHub 仓库](https://github.com/orsharir/github-mathjax)）和 [github-mathjax-firefox](https://github.com/traversaro/github-mathjax-firefox/releases/download/v0.2.3/github_with_mathjax-0.2.3.xpi)（[GitHub 仓库](https://github.com/traversaro/github-mathjax-firefox)）插件即可查看完美显示的公式了。

## 理论方法

本仓库包含了一些 Hankel 变换的例子，牵涉到了第一类和第二类 Bessel 函数、其一阶和二阶导数及其零点的计算。下文中将给出一些与 Bessel 函数导数和零点、以及 Hankel 变换的计算有关的数学公式。

### Bessel 函数的导数

下文中以 $Z_\nu(x)$ 表示 $J_\nu(x)$ 或 $Y_\nu(x)$，根据 Bessel 函数的递推关系：

$$ x Z_{\nu - 1}(x) + x Z_{\nu + 1}(x) = 2\nu Z_\nu(x) $$

$$ Z_\nu'(x) = \frac{\nu}{x} Z_\nu(x) - Z_{\nu + 1}(x) $$

可得

$$ Z_{\nu + 2}(x) = \frac{2(\nu + 1)}{x} Z_{\nu + 1}(x) - Z_\nu(x) $$

$$ Z_{\nu + 1}'(x) = \frac{\nu + 1}{x} Z_{\nu + 1}(x) - Z_{\nu + 2}(x) $$

不难得出：

$$ \begin{align*} 
     Z_\nu''(x) & = \frac{d}{dx} \left[ \frac{\nu}{x} Z_\nu(x) - Z_{\nu + 1}(x) \right] \newline
                & = - \frac{\nu}{x^2} Z_\nu(x) + \frac{\nu}{x} Z_\nu'(x) - Z_{\nu + 1}'(x) \newline
                & = - \frac{\nu}{x^2} Z_\nu(x) + \frac{\nu}{x} \left[ \frac{\nu}{x} Z_\nu(x) - Z_{\nu + 1}(x) \right] - \left[ \frac{\nu + 1}{x} Z_{\nu + 1}(x) - \frac{2(\nu + 1)}{x} Z_{\nu + 1}(x) + Z_\nu(x) \right] \newline
                & = \left( - \frac{\nu}{x^2} + \frac{\nu^2}{x^2} - 1 \right) Z_\nu(x) + \left( - \frac{\nu}{x} - \frac{\nu + 1}{x} + \frac{2\nu + 2}{x} \right) Z_{\nu + 1}(x) \newline
                & = \frac{1}{x^2} \left[ (\nu^2 - \nu - x^2) Z_\nu(x) + x Z_{\nu + 1}(x) \right]
   \end{align*} $$

和

$$ \begin{align*} 
     Z_\nu''(x) & = \left( \frac{\nu^2}{x^2} - 1 \right) Z_\nu(x) - \frac{1}{x} \left[ \frac{\nu}{x} Z_\nu(x) - Z_{\nu + 1}(x) \right] \newline
                & = \left( \frac{\nu^2}{x^2} - 1 \right) Z_\nu(x) - \frac{1}{x} Z_\nu'(x)
   \end{align*} $$

### Bessel 函数的零点

#### 固定点迭代

根据 [Wikipedia](https://en.wikipedia.org/wiki/Fixed-point_iteration#Applications)，可用于求解 Bessel 函数零点的固定点迭代式为：

$$ \begin{align*}
     x_{n + 1} & = x_n - \frac{Z_\nu(x_n)}{Z_\nu'(x_n)} \newline
               & = x_n - \frac{x_n Z_\nu(x_n)}{\nu Z_\nu(x_n) - x_n Z_{\nu + 1}(x_n)}
   \end{align*} $$

#### Halley 方法

基于以上公式，以下乘积的结果为：

$$ Z_\nu(x) Z_\nu'(x) = Z_\nu(x) \left[ \frac{\nu}{x} Z_\nu(x) - Z_{\nu + 1}(x) \right] = \frac{1}{x^2} x Z_\nu(x) \left[ \nu Z_\nu(x) - x Z_{\nu + 1}(x) \right] $$

$$ [Z_\nu'(x)]^2 = \frac{\nu^2}{x^2} Z_\nu(x) - \frac{2\nu}{x} Z_\nu(x) Z_{\nu + 1}(x) + [Z_{\nu + 1}(x)]^2 = \frac{1}{x^2} \left[ \nu^2 Z_\nu(x) - 2\nu x Z_\nu(x) Z_{\nu + 1}(x) + x^2 [Z_{\nu + 1}(x)]^2 \right] $$

$$ Z_\nu(x) Z_\nu''(x) = \frac{1}{x^2} \left\lbrace (\nu^2 - \nu - x^2) [Z_\nu(x)]^2 + x Z_\nu(x) Z_{\nu + 1}(x) \right\rbrace $$

根据 [Wikipedia](https://en.wikipedia.org/wiki/Halley's_method#Method)，可用于求解 Bessel 函数零点的 Halley 法迭代式为：

$$ \begin{align*}
     x_{n + 1} & = x_n - \frac{2 Z_\nu(x_n) Z_\nu'(x_n) }{ 2 [Z_\nu'(x_n)]^2 - Z_\nu(x_n) Z_\nu''(x_n) } \newline
               & = x_n - \frac{ 2 x_n Z_\nu(x_n) [ \nu Z_\nu(x_n) - x_n Z_{\nu + 1}(x_n) ] }{ 2 x_n^2 [Z_{\nu + 1}(x_n)]^2 - (4\nu + 1) x_n Z_\nu(x_n) Z_{\nu + 1}(x_n) + (\nu^2 + \nu + x_n^2) [Z_\nu(x_n)]^2}
   \end{align*} $$

### Hankel 变换

#### 定义

根据 [Wikipedia](https://en.wikipedia.org/wiki/Hankel_transform#Definition)，函数 $f(r)$ 的 $\nu$ 阶 Hankel 变换为：

$$ F_\nu(k) = \int_0^\infty f(r) J_\nu(k r) r dr $$

#### (Guptasarma and Singh, 1997)

文中的 Hankel 变换类型为（式 2 和 3）：

$$ f(r) = \int_0^\infty K(\lambda) J_\nu(r \lambda) d\lambda \approx \frac{1}{r} \sum_{i = 1}^n K(\lambda_i) W_i^\nu $$

此近似公式**仅**适用于基于第一类零阶或一阶 Bessel 函数的 Hankel 变换。

文中给出了一些可用于检验程序的标准 Hankel 变换示例：

| 公式编号 |                         $K(\lambda)$                         | $J_\nu$ |                            $f(r)$                            |
| :------: | :----------------------------------------------------------: | :-----: | :----------------------------------------------------------: |
|    4     |                      $e^{ - c \lambda}$                      |  $J_0$  |                 $\frac{1}{\sqrt{c^2 + r^2}}$                 |
|    5     |                 $\lambda e^{ - c \lambda^2}$                 |  $J_0$  |                $\frac{1}{2c} e^{ - r^2/(4c)}$                |
|    6     |                  $\lambda e^{ - c \lambda}$                  |  $J_0$  |               $\frac{c}{\sqrt{(c^2 + r^2)^3}}$               |
|    7     | $\lambda e^{ - c \lambda} + \alpha \lambda^2 e^{ - c \lambda^2}$ |  $J_1$  | $\frac{r}{\sqrt{(c^2 + r^2)^3}} + \alpha \frac{r e^{ - r^2/(4c)}}{4 c^2}$ |
|    8     |                  $\lambda e^{ - c \lambda}$                  |  $J_1$  |               $\frac{r}{\sqrt{(c^2 + r^2)^3}}$               |
|    9     |                $\lambda^2 e^{ - c \lambda^2}$                |  $J_1$  |              $\frac{r e^{ - r^2/(4c)}}{4 c^2}$               |
|    10    |                      $e^{ - c \lambda}$                      |  $J_1$  |      $\frac{\sqrt{c^2 + r^2} - c}{r \sqrt{c^2 + r^2}}$       |

#### (Ogata, 2005)

文中的 Hankel 变换型积分的近似公式为（式 5.2）：

$$ \begin{align*}
  \int_0^\infty f(x) J_\nu(x) dx & \approx \pi \sum_{k = 1}^{\infty} \omega_{\nu k} f(x) J_\nu(x) \psi'(t) \newline
                                 & = \pi \sum_{k = 1}^{\infty} \omega_{\nu k} f \left( \frac{\pi}{h} \psi(h \xi_{\nu k}) \right) J_\nu \left( \frac{\pi}{h} \psi(h \xi_{\nu k}) \right) \psi'(h \xi_{\nu k}) \newline
\end{align*}$$

其中 $ x = \pi / h \cdot \psi(t) $，$ \psi(t) = t \tanh(\pi / 2 \cdot \sinh t) $ 且 $ t = h \xi_{\nu k} $.

设 $ x = r \lambda $，则 $ \lambda = x/r $ 且 $ d\lambda = 1/r \cdot dx $，又有 $ x = \pi/h \cdot \psi(h \xi_{\nu k}) $，可得：

$$ \begin{align*}
     \int_0^\infty f(\lambda) J_\nu(r \lambda) \lambda d\lambda & = \int_0^\infty f \left( \frac{x}{r} \right) J_\nu(x) \frac{x}{r} \cdot \frac{1}{r} dx \newline
                                                                & = \frac{1}{r^2} \int_0^\infty f \left( \frac{x}{r} \right) J_\nu(x) x dx \newline
                                                                & \approx \frac{1}{r^2} \left[ \pi \sum_{k = 1}^{\infty} \omega_{\nu k} f \left( \frac{x}{r} \right) J_\nu(x) \psi'(h \xi_{\nu k}) \right] x
   \end{align*} $$

此公式适用于基于第一类**任意**阶 Bessel 函数的 Hankel 变换。

## 具体实现

### example_Functions

[example_Functions](example_Functions.F90) 是一个包含了 (Guptasarma and Singh, 1997) 文中 Hankel 变换的全部例子的 `fortran` 模块。

### Hankel_Transform

[Hankel_Transform](Hankel_Transform.F90) 是一个基于 (Guptasarma and Singh, 1997) 文中方法实现 Hankel 变换的 `fortran` 模块。

由于本源文件中不包含主程序，例如，可通过命令 `$ gfortran -c Hankel_Transform.F90` 来编译此程序，但无法运行之。

### Bessel_Function

[Bessel_Function](Bessel_Function.F90) 是一个用于计算第一类和第二类 Bessel 函数、其一阶和二阶导数及其原函数和一阶导数的零点的 `fortran` 模块。该模块主要由 [mjyzo.f90](http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/mjyzo_f90.txt) 修改而来。

在本源文件中有两个预处理宏：
- **MOREAU** 用于调用 Moreau 方法来计算第一类和第二类 Bessel 函数及其一阶导数；如未定义，则调用 `fortran` 内部子程序来计算。
- **HALLEY** 用于调用 Halley 方法来计算第一类和第二类 Bessel 函数的零点；如未定义，则调用固定点迭代法来计算。

由于本源文件不包含主程序，例如，可通过命令 `$ gfortran -DHALLEY -c Bessel_Function.F90` 来编译此程序，但无法运行之。

此外， [Bessel Zero Solver](https://ww2.mathworks.cn/matlabcentral/fileexchange/48403-bessel-zero-solver) 是一个有效且可靠的用于基于 Halley 方法求解第一类和第二类 Bessel 函数零点的 `matlab` 函数。

### Guptasarma_1997

[Guptasarma_1997](Guptasarma_1997.F90) 是一个基于 [Hankel_Transform 模块](Hankel_Transform.F90)的主程序，它实现了 (Guptasarma and Singh, 1997) 的方法论，并通过调用 [example_Functions 模块](example_Functions.F90)包含了该文中的所有例子。

在本源文件中有一个预处理宏：**FAST** 用于采用更少的采样点来更快速地完成 Hankel 变换；如未定义，则采用更多的采样点来更准确地完成。

例如，可通过命令 `$ gfortran -DFAST example_Functions.F90 Hankel_Transform.F90 Guptasarma_1997.F90 -o Guptasarma` 来编译此源文件，且以命令 `$ ./Guptasarma` 运行之。

### Ogata_2005

[Ogata_2005](Ogata_2005.F90) 是一个基于 [Bessel_Function 模块](Bessel_Function.F90)的主程序，它将 (Ogata, 2005) 的方法论应用到了 Hankel 变换上，并通过调用 [example_Functions 模块](example_Functions.F90)包含了 (Guptasarma and Singh, 1997) 中的所有例子。

例如，可通过命令 `$ gfortran -DHALLEY example_Functions.F90 Bessel_Function.F90 Ogata_2005.F90 -o Ogata` 来完成此程序的编译，并以命令 `$ ./Ogata` 运行之。

## 参考文献

- Guptasarma and Singh, 1997. **New digital linear filters for Hankel J_0 and J_1 transforms**. Geophysical Prospecting, 45, 745-762. [https://doi.org/10.1046/j.1365-2478.1997.500292.x](https://onlinelibrary.wiley.com/doi/abs/10.1046/j.1365-2478.1997.500292.x)
- Ogata, 2005. **A numerical integration formula based on the Bessel functions**. Publ. RIMS, 41, 949-970. [https://doi.org/10.2977/prims/1145474602](https://www.ems-ph.org/journals/show_abstract.php?issn=0034-5318&vol=41&iss=4&rank=8)

