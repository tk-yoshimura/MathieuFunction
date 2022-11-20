// Copyright (c) T.Yoshimura 2022
// https://github.com/tk-yoshimura

namespace MathieuMP {
    public static class EigenFP64 {
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

        // If a matches the true value, return 0.
        public static double Fraction(EigenFunc func, int n, double q, double a, int terms = 32) {
            return func switch {
                EigenFunc.A => FractionA(n, q, a, terms),
                EigenFunc.B => FractionB(n, q, a, terms),
                _ => throw new ArgumentException(nameof(func)),
            };
        }

        // The true value is obtained by the secant method.
        public static (double value, double error) SearchFit(EigenFunc func, int n, double q, double a, int maxiter = 16) {
            if (q == 0) {
                return (0d, 0d);
            }

            double a0 = a, d0 = Fraction(func, n, q, a0), da = double.NaN;

            for (int iter = 0; iter < maxiter; iter++) {
                double ah = Math.CopySign(Math.Max(Math.ScaleB(1, -40), Math.ScaleB(Math.Abs(a0), -10)), a0);
                double a1 = a0 + ah, d1 = Fraction(func, n, q, a1);

                double dh = d1 - d0;

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
                    break;
                }
            }

            return (a0, da);
        }
    }
}
