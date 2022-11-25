// Copyright (c) T.Yoshimura 2022
// https://github.com/tk-yoshimura

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
        public static double Fraction(EigenFunc func, int n, double q, double a, int terms = 32) {
            return func switch {
                EigenFunc.A => FractionA(n, q, a, terms),
                EigenFunc.B => FractionB(n, q, a, terms),
                _ => throw new ArgumentException(nameof(func)),
            };
        }

        /// <summary>
        /// The true value is obtained by the secant method.
        /// NOTE: a is within the radii of convergence.
        /// </summary>
        private static (double value, double error) SearchFit(EigenFunc func, int n, double q, double a, int maxiter = 16) {
            if (q == 0) {
                return (0d, 0d);
            }

            double a0 = a, d0 = Fraction(func, n, q, a0), da = double.NaN;

            for (int iter = 0; iter < maxiter; iter++) {
                double ah = Math.CopySign(Math.Max(Math.ScaleB(1, -40), Math.ScaleB(Math.Abs(a0), -10)), a0);
                double a1 = a0 + ah, d1 = Fraction(func, n, q, a1);

                double dh = d1 - d0;

                //Console.WriteLine($"{d0}, {d1}");

                if (dh == 0) {
                    break;
                }
                da = ah / dh * d0;
                a0 -= da;

                if (Math.Abs(da / a0) <= 1e-15) {
                    break;
                }

                d0 = Fraction(func, n, q, a0);

                if (d0 == 0) {
                    return (a0, 0d);
                }
            }

            return (a0, da);
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

            double y = (func, n) switch {
                (EigenFunc.A, 0 or 1) => (q <= 4) ? NearPeak(func, n, q) : Asymptotic(func, n, q),
                (EigenFunc.A, 2) => (q <= 3) ? NearPeak(func, n, q) : Asymptotic(func, n, q),
                (EigenFunc.A, 3) => (q <= 6.25) ? NearPeak(func, n, q) : Asymptotic(func, n, q),
                (EigenFunc.A, <= 64) => (q <= 0.49 * n * n)
                                        ? NearPeak(func, n, q)
                                        : (q <= 0.7225 * n * n)
                                            ? (NearPeak(func, n, q) + Asymptotic(func, n, q)) / 2
                                            : Asymptotic(func, n, q),
                (EigenFunc.A, _) => ((n & 1) == 0)
                                        ? Value(func, n / 2, q / 2, zero_shift: false).value * 4
                                        : Value(func, n / 2, q / 2, zero_shift: false).value * (n * n / (double)(n / 2 * n / 2)),
                (EigenFunc.B, 1) => (q <= 4) ? NearPeak(func, n, q) : Asymptotic(func, n, q),
                (EigenFunc.B, 2) => (q <= 5) ? NearPeak(func, n, q) : Asymptotic(func, n, q),
                (EigenFunc.B, 3) => (q <= 6.25) ? NearPeak(func, n, q) : Asymptotic(func, n, q),
                (EigenFunc.B, <= 64) => (q <= 0.49 * n * n)
                                        ? NearPeak(func, n, q)
                                        : (q <= 0.7225 * n * n)
                                            ? (NearPeak(func, n, q) + Asymptotic(func, n, q)) / 2
                                            : Asymptotic(func, n, q),
                (EigenFunc.B, _) => ((n & 1) == 0)
                                        ? Value(func, n / 2, q / 2, zero_shift: false).value * 4
                                        : Value(func, n / 2, q / 2, zero_shift: false).value * (n * n / (double)((n / 2) * (n / 2))),
                _ => throw new ArgumentException(nameof(func))
            };

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
                (EigenFunc.A, _) => nz_largen(n, h),
                (EigenFunc.B, _) => nz_largen(n + 1, h),
                _ => throw new ArgumentException(nameof(func))
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

            double y = (2 * (-q + s * u) - n * n) + c0 + v * (c1 + v * (c2 + v * (c3 + v * (c4 + v * c5))));

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
        public static (double value, double error) Value(EigenFunc func, int n, double q, bool zero_shift = false) {
            double a0 = InitialValue(func, n, q);
            (double value, double error) = SearchFit(func, n, q, a0);

            if (!zero_shift) {
                value += n * n;
            }

            return (value, error);
        }
    }
}
