# MathieuFunction

## Definition

## Eigenvalues

![eigen summary](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/eigen_plot.svg)  

### Continued Fraction and Convergence

[DLMF 28.6](http://dlmf.nist.gov/28.6)  

![a even](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/eigen_a_even.svg)  
![b even](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/eigen_b_even.svg)  
![a odd](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/eigen_a_odd.svg)  
![b odd](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/eigen_b_odd.svg)  

### Normalized

![eigen normalized](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/eigen_plot_normalized.svg)  

### Average and Difference

![eigen mean sub](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/eigen_mean_sub.svg)  
![eigen mean raw](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/eigen_plot_mean_raw.svg)  
![eigen sub raw](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/eigen_plot_sub_raw.svg)  

### Rational Coefficient Series

![eigen rcoef md](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/eigen_rcoef_md.svg)  
![eigen rcoef md grad0](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/eigen_rcoef_md_grad0.svg)  
![eigen rcoef m](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/eigen_plot_m.svg)  
![eigen rcoef d](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/eigen_plot_d.svg)  

[Numeric Table (64 digits)](https://github.com/tk-yoshimura/MathieuFunction/blob/main/results)  

### Continued Fraction Terms

128bit  
![eigen cfrac 128bit](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/needs_frac_log2_mp4.svg)  
256bit  
![eigen cfrac 256bit](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/needs_frac_log2_mp8.svg)  
512bit  
![eigen cfrac 512bit](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/needs_frac_log2_mp16.svg)  
1024bit  
![eigen cfrac 1024bit](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/needs_frac_log2_mp32.svg)  

## Fourier Series
![coef fourier](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/coef_fourier.svg)  

### Coefficient Matrix
The relationship between the coefficients is expressed using a triple diagonal matrix T.  
![coef matrix](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/coef_matrix.svg)  
![coef param](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/coef_param.svg)  

### Backward Recurse
Backward recursion is used as an approximate method to solve the coefficient sequence.  
![coef backrecur](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/coef_backrecur.svg)  

### Solve Matrix
In backward recursion, terms may cancel each other and produce a value close to zero, although this is a rare case.  
In this case, it must be solved in pieces as a linear problem because the digits drop out and subsequent values become inaccurate.  
![coef det](https://github.com/tk-yoshimura/MathieuFunction/blob/main/figures/coef_det.svg)  
