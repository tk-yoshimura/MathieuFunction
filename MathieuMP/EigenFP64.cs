// Copyright (c) T.Yoshimura 2022
// https://github.com/tk-yoshimura

using System;
using System.Diagnostics;

namespace MathieuMP {
    public static class EigenFP64 {

        /// <summary>
        /// Continue fraction eigen A.
        /// </summary>
        /// <remarks>DLMF 28.6</remarks>
        private static double FractionA(int n, double q, double a, int terms) {
            if (n > 16384 || !(q >= 0)) {
                throw new ArgumentException($"{nameof(n)}, {nameof(q)}");
            }

            int n_sq = n * n;

            double h = q * q;
            double u = ((n & 1) == 0)
                ? (n == 0) ? 0d : 2d / (a + n_sq)
                : (n == 1) ? 1d / q : 1d / (a + (n_sq - 1) - q);
            double v = 0d;

            for (int k = 2 + (n & 1); k < n; k += 2) {
                u = 1d / (a + (n_sq - k * k) - u * h);
            }
            for (int k = terms * 2 + n; k > n; k -= 2) {
                v = 1d / (a + (n_sq - k * k) - v * h);
            }

            double y = (((n > 0) ? u : v) + v) * h - a;

            return y;
        }

        /// <summary>
        /// Continue fraction eigen B.
        /// </summary>
        /// <remarks>DLMF 28.6</remarks>
        private static double FractionB(int n, double q, double a, int terms) {
            if (n < 1 || n > 16384 || !(q >= 0)) {
                throw new ArgumentException($"{nameof(n)}, {nameof(q)}");
            }

            int n_sq = n * n;

            double h = q * q;
            double u = ((n & 1) == 0)
                ? 0d
                : (n == 1) ? -1d / q : 1d / (a + (n_sq - 1) + q);
            double v = 0d;

            for (int k = 2 + (n & 1); k < n; k += 2) {
                u = 1d / (a + (n_sq - k * k) - u * h);
            }
            for (int k = terms * 2 + n; k > n; k -= 2) {
                v = 1d / (a + (n_sq - k * k) - v * h);
            }

            double y = (u + v) * h - a;

            return y;
        }

        /// <summary>
        /// If a matches the true value, return 0.
        /// </summary>
        public static double Fraction(EigenFunc func, int n, double q, double a, int terms = -1) {
            terms = (terms < 0) ? FracTerms(func, n, q) : terms;

            double y = func switch {
                EigenFunc.A => FractionA(n, q, a, terms),
                EigenFunc.B => FractionB(n, q, a, terms),
                _ => throw new ArgumentException(nameof(func)),
            };

            return y;
        }

        /// <summary>
        /// Obtain true value by the binary search and secant method.
        /// NOTE: a is within the radii of convergence.
        /// </summary>
        public static (double value, double score, bool is_convergence) SearchFit(EigenFunc func, int n, double q, double a, int frac_terms = -1) {
            if (q == 0) {
                return (0, 1, is_convergence: true);
            }

            frac_terms = (frac_terms < 0) ? FracTerms(func, n, q) : frac_terms;

            double h = Math.Max(1, n * n) / 32d;
            double truncation_thr = 2 + Math.Max(1, n * n) * 0.1;
            double heuristics_err = Math.Max(a * 1e-2, func == EigenFunc.A ? 4.02731e-3 * n * n + 2.0 : 3.80915e-3 * n * n + 2.5);

            (double ar, bool ar_convergence, double ar_score) = RootFinder.Search((a) => Fraction(func, n, q, a, frac_terms), a, h, truncation_thr);
            (double ap, bool ap_convergence, _) = RootFinder.Search((a) => 1 / Fraction(func, n, q, a, frac_terms), a, h, truncation_thr);

            (double a_likelihood, double score_likelihood) = (ar_convergence && Math.Abs(a - ar) < heuristics_err) ? (ar, ar_score) : (double.NaN, 0);

            if (ap_convergence) {
                (double apm, bool apm_convergence, double apm_score) = RootFinder.Search((a) => Fraction(func, n, q, a, frac_terms), Math.BitDecrement(ap), h, truncation_thr, SearchDirection.Minus);
                (double app, bool app_convergence, double app_score) = RootFinder.Search((a) => Fraction(func, n, q, a, frac_terms), Math.BitIncrement(ap), h, truncation_thr, SearchDirection.Plus);

                if (apm_convergence && apm_score > ar_score * 0.5 && Math.Abs(a - apm) < heuristics_err) {
                    (a_likelihood, score_likelihood) = Math.Abs(a - a_likelihood) < Math.Abs(a - apm) ? (a_likelihood, score_likelihood) : (apm, apm_score);
                }
                if (app_convergence && app_score > ar_score * 0.5 && Math.Abs(a - app) < heuristics_err) {
                    (a_likelihood, score_likelihood) = Math.Abs(a - a_likelihood) < Math.Abs(a - app) ? (a_likelihood, score_likelihood) : (app, app_score);
                }
            }

            bool is_convergence = !double.IsNaN(a_likelihood);
            if (!is_convergence && Math.Abs(a - ap) < heuristics_err) {
                a_likelihood = ap;
            }

            return (a_likelihood, score_likelihood, is_convergence);
        }

        /// <summary>
        /// Initial value for the secant method.
        /// </summary>
        public static double InitialValue(EigenFunc func, int n, double q) {
            if (func == EigenFunc.A && n < 0) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }
            if (func == EigenFunc.B && n < 1) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }

            if (q < 1d / 32) {
                return NearZero(func, n, q);
            }

            static double bump(double x, double s, double t, double a, double b) {
                double c = (x - s) / (t - s);

                if (c < 0.001) {
                    return a;
                }
                if (c > 0.988) {
                    return b;
                }

                double w = 1d / (Math.Exp(1d / c - 1d / (1d - c)) + 1d);

                return a * (1d - w) + b * w;
            }

            double y = (func, n) switch {
                (EigenFunc.A, 0) => (q <= 2.08) ? NearPeak(func, n, q) : Asymptotic(func, n, q),
                (EigenFunc.A, 1) => (q <= 4.00) ? NearPeak(func, n, q) : Asymptotic(func, n, q),
                (EigenFunc.A, 2) => (q <= 6.72) ? NearPeak(func, n, q) : Asymptotic(func, n, q),
                (EigenFunc.A, 3) => (q <= 10.27) ? NearPeak(func, n, q) : Asymptotic(func, n, q),
                (EigenFunc.A, 4) => bump(q, 5.280, 13.92, NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.A, 5) => bump(q, 14.50, 21.50, NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.A, 6) => bump(q, 16.20, 29.16, NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.A, 7) => bump(q, 15.68, 47.53, NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.A, 8) => bump(q, 16.00, 68.48, NearPeak(func, n, q), Asymptotic(func, n, q)),
                
                (EigenFunc.B, 1) => bump(q, 0.601, 5.734, NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.B, 2) => bump(q, 1.375, 9.250, NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.B, 3) => bump(q, 2.390, 13.78, NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.B, 4) => bump(q, 4.480, 8.640, NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.B, 5) => bump(q, 3.250, 20.50, NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.B, 6) => bump(q, 8.280, 29.52, NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.B, 7) => bump(q, 13.72, 39.69, NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.B, 8) => bump(q, 20.48, 51.20, NearPeak(func, n, q), Asymptotic(func, n, q)),

                (EigenFunc.A, < 11) => bump(q / (n * n), 
                                         0.21628 + (10 - n) * (10 - n) * 1.36207e-2,
                                         1.01724 - (10 - n) * (10 - n) * 1.69540e-2, 
                                         NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.A, < 50) => bump(q / (n * n), 
                                         0.43546 - (50 - n) * (50 - n) * 1.13079e-4,
                                         0.82488 + (50 - n) * (50 - n) * 1.41561e-4, 
                                         NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.A, _) => bump(q / (n * n), 0.43546, 0.82488, NearPeak(func, n, q), Asymptotic(func, n, q)),

                (EigenFunc.B, < 30) => bump(q / (n * n), 
                                         0.51191 - (30 - n) * (30 - n) * (30 - n) * 1.89525e-5,
                                         0.65511 + (30 - n) * (30 - n) * 3.02459e-4, 
                                         NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.B, < 50) => bump(q / (n * n), 
                                         0.42715 + (50 - n) * 1.91587e-3,
                                         0.81275 - (50 - n) * (50 - n) * 3.69580e-4, 
                                         NearPeak(func, n, q), Asymptotic(func, n, q)),
                (EigenFunc.B, _) => bump(q / (n * n), 0.42715, 0.81275, NearPeak(func, n, q), Asymptotic(func, n, q)),
                _ => throw new ArgumentException(nameof(func))
            };

            return y;
        }

        public static double InitialValueTest(EigenFunc func, int n, double q, double s, double t) {
            if (func == EigenFunc.A && n < 0) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }
            if (func == EigenFunc.B && n < 1) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }

            if (q < 1d / 32) {
                return NearZero(func, n, q);
            }

            static double bump(double x, double s, double t, double a, double b) {
                double c = (x - s) / (t - s);

                if (c < 0.001) {
                    return a;
                }
                if (c > 0.988) {
                    return b;
                }

                double w = 1d / (Math.Exp(1d / c - 1d / (1d - c)) + 1d);

                return a * (1d - w) + b * w;
            }

            double r = 1d / Math.Max(1, n * n);

            double y = bump(q * r, s, t, NearPeak(func, n, q), Asymptotic(func, n, q));

            return y;
        }

        /// <summary>
        /// Root finding a + b x + c x^2 + x^3 = 0
        /// </summary>
        public static double RootCubic(double a, double b, double c) {
            double p = (b * 3 - c * c) / 9;
            double q = (a * 27 - c * (b * 9 - c * c * 2)) / 54;

            double d = p * p * p + q * q;

            if (d >= 0) {
                double r = Math.Sqrt(d), u = -q + r, v = -q - r;
                double x = Math.Cbrt(u) + Math.Cbrt(v) - c / 3;

                return x;
            }
            else {
                double t = Math.Acos(-q / Math.Sqrt(-p * p * p));
                double x = 2 * Math.Sqrt(-p) * Math.Cos((t + 4 * Math.PI) / 3) - c / 3;

                return x;
            }
        }

        /// <summary>
        /// Near zero approx eigen values. (q less than 1/32)
        /// </summary>
        private static double NearZero(EigenFunc func, int n, double q) {
            double h = q * q;

            double y = (func, n) switch {
                (EigenFunc.A, 0) => h * (-1d / 2 + h * (7d / 128)),
                (EigenFunc.A, 1) => q + h * (-1d / 8 + h * (-1d / 64)),
                (EigenFunc.A, 2) => h * (5d / 12 + h * (-763d / 13824)),
                (EigenFunc.A, 3) => h * (1d / 16 + q * (1d / 64 + q * (13d / 20480))),
                (EigenFunc.B, 1) => -q + h * (-1d / 8 + h * (1d / 64)),
                (EigenFunc.B, 2) => h * (-1d / 12 + h * (5d / 13824)),
                (EigenFunc.B, 3) => h * (1d / 16 + q * (-1d / 64 + q * (13d / 20480))),
                _ => h / (2 * (n * n - 1))
            };

            return y;
        }

        /// <summary>
        /// Near peak (q/(n^2) in less than 1-2)
        /// </summary>
        /// <remarks>Fayez A. Alhargan (2000)</remarks>
        private static double NearPeak(EigenFunc func, int n, double q) {
            double h = q * q;

            static double nz_largen(int n, double h) {
                double n_sq = n * n, n_sq_m1 = n_sq - 1, n_sq_m4 = n_sq - 4, n_sq_m9 = n_sq - 9;

                double c2 = 1 / (2 * n_sq_m1);
                double c4 = (7 + n_sq * 5) / (32 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m4);
                double c6 = (29 + n_sq * (58 + n_sq * 9)) / (64 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m4 * n_sq_m9);

                return h * (c2 + h * (c4 + h * c6));
            }

            double y = (func, n) switch {
                (EigenFunc.A, 0) => 2 - Math.Sqrt(4 + h * 2),
                (EigenFunc.A, 1) => 4 + (q - Math.Sqrt(64 + q * (-16 + q * 5))) / 2,
                (EigenFunc.A, 2) => RootCubic(h * 20, -48 - h * 3, -8),
                (EigenFunc.A, 3) => RootCubic(h * (8 + q), -128 + q * 16 - h * 2, -q - 8),
                (EigenFunc.B, 1) => 4 - (q + Math.Sqrt(64 + q * (16 + q * 5))) / 2,
                (EigenFunc.B, 2) => 6 - Math.Sqrt(36 + h),
                (EigenFunc.B, 3) => RootCubic(h * (8 - q), -128 - q * 16 - h * 2, q - 8),
                (EigenFunc.A, 4) => nz_largen(n, h) + h * h * (1 / 2304d),
                (EigenFunc.B, 4) => nz_largen(n, h) - h * h * (1 / 2304d),
                _ => nz_largen(n, h)
            };

            return y;
        }

        /// <summary>
        /// Asymptotic eigen values. (q/(n^2) greater than 1-2)
        /// </summary>
        /// <remarks>DLMF 28.6</remarks>
        private static double Asymptotic(EigenFunc func, int n, double q) {
            double u = Math.Sqrt(q), v = 1d / u;
            double s = func == EigenFunc.A ? (2 * n + 1) : (2 * n - 1), s_sq = s * s;

            double c0 = -Math.ScaleB(1 + s_sq, -3);
            double c1 = -Math.ScaleB(s * (3 + s_sq), -7);
            double c2 = -Math.ScaleB(9 + s_sq * (34 + s_sq * 5), -12);
            double c3 = -Math.ScaleB(s * (405 + s_sq * (410 + s_sq * 33)), -17);
            double c4 = -Math.ScaleB(486 + s_sq * (2943 + s_sq * (1260 + s_sq * 63)), -20);
            double c5 = -Math.ScaleB(s * (41607 + s_sq * (69001 + s_sq * (15617 + s_sq * 527))), -25);

            double c6 = (n < 4)
                ? 0d
                : func == EigenFunc.A
                    ? (-5.682576740891 * Math.Pow(n, 4.7) - 0.005055889635 * Math.Pow(n, 7.9))
                    : (+0.517862332643 * Math.Pow(n, 5.5) - 0.002906986648 * Math.Pow(n, 8.0));

            double y = (2 * (-q + s * u) - n * n) + c0 + v * (c1 + v * (c2 + v * (c3 + v * (c4 + v * (c5 + v * c6)))));

            return y;
        }

        /// <summary>
        /// Compute eigan value
        /// </summary>
        /// <param name="func">eigen function</param>
        /// <param name="n">order</param>
        /// <param name="q">q</param>
        /// <param name="zero_shift">remove bias (=n^2)</param>
        /// <returns></returns>
        public static (double value, double score, bool is_convergence) Value(EigenFunc func, int n, double q, int frac_terms = -1, bool zero_shift = false) {
            frac_terms = (frac_terms < 0) ? FracTerms(func, n, q) : frac_terms;

            double a0 = InitialValue(func, n, q);
            (double value, double score, bool is_convergence) = SearchFit(func, n, q, a0, frac_terms);

            if (!zero_shift) {
                value += n * n;
            }

            return (value, score, is_convergence);
        }

        public static (double value, int terms) ConvergenceFracTerms(EigenFunc func, int n, double q, int init_terms) {
            double a0 = Value(func, n, q, frac_terms: init_terms, zero_shift: true).value;
            double a1 = Value(func, n, q, frac_terms: init_terms + 1, zero_shift: true).value;
            double err01 = Math.Abs(a0 - a1);

            for (int frac_terms = init_terms; frac_terms <= 4096; frac_terms += 2) {
                double a2 = Value(func, n, q, frac_terms + 2, zero_shift: true).value;
                double err12 = Math.Abs(a1 - a2);

                double relative_err = err01 / (Math.Abs(a0) + double.Epsilon);

                if ((relative_err < 1e-14) && (err01 <= err12) && ((a0 >= a1 && a1 <= a2) || (a0 <= a1 && a1 >= a2))) {
                    return (a2, frac_terms);
                }

                (a0, a1, err01) = (a1, a2, err12);
            }

            return (double.NaN, -1);
        }

        public static int FracTerms(EigenFunc func, int n, double q) {
            if (!(n >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }
            if (!(q >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(q));
            }

            double intercept = 2.13317 + Math.Sqrt(n + 1) * 0.104619;
            double slope = 0.239519 * Math.Pow(n + 1, 0.0192171);

            int terms = Math.Max(n, (int)Math.Ceiling(Math.Pow(2, intercept + slope * Math.Log2(q))));

            return terms;
        }
    }
}
