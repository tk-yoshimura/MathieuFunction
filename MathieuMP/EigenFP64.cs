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
        private static (double value, double error) SearchFit(EigenFunc func, int n, double q, double a, int maxiter = 16) {
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

                Console.WriteLine($"{ah/dh},{d0 - d1},{a0}");

                d0 = Fraction(func, n, q, a0);

                if (d0 == 0) {
                    return (a0, 0d);
                }
            }

            return (a0, da);
        }

        // Initial value for the secant method.
        public static double InitialValue(EigenFunc func, int n, double q) {
            if (!Enum.IsDefined(func)) {
                throw new ArgumentException(nameof(func));
            }

            if (func == EigenFunc.A && n < 0) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }
            if (func == EigenFunc.B && n < 1) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }

            double y = (q < 2 * Math.Max(1, n * n)) ? NearZero(func, n, q) : Asymptotic(func, n, q);
            y = Math.Min(y, UpperBound(func, n, q));

            return y;
        }

        private static double NearZero(EigenFunc func, int n, double q) {
            double h = q * q;

            double y = (func, n) switch {
                (EigenFunc.A, 0) => h * (-1d / 2),
                (EigenFunc.A, 1) => q + h * (-1d / 8),
                (EigenFunc.B, 1) => -q + h * (-1d / 8),
                (EigenFunc.A, 2) => h * (5d / 12),
                (EigenFunc.B, 2) => h * (-1d / 12),
                _ => h / (2 * (n * n - 1))
            };

            return y;
        }

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

        private static double UpperBound(EigenFunc func, int n, double q) {
            double maxrate = (func, n) switch {
                (EigenFunc.A, 0) => 0,
                (EigenFunc.A, 1) => 1.52169,
                (EigenFunc.A, 2) => 1.03126,
                (EigenFunc.A, _) => 0.846,
                (EigenFunc.B, <= 2) => 0,
                (EigenFunc.B, 3) => 0.0291174,
                (EigenFunc.B, 4) => 0.0883491,
                (EigenFunc.B, _) => 0.4,
                (_, _) => throw new ArgumentException(nameof(func))
            };

            double y = n * n * maxrate;

            return y;
        }

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
