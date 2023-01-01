// Copyright (c) T.Yoshimura 2022
// https://github.com/tk-yoshimura

using MultiPrecision;

namespace MathieuMP {
    public static class EigenMP<N> where N : struct, IConstant {
        /// <summary>
        /// Continue fraction eigen A.
        /// </summary>
        /// <remarks>DLMF 28.6</remarks>
        private static MultiPrecision<N> FractionA(int n, MultiPrecision<N> q, MultiPrecision<N> a, int terms) {
            if (n > 16384 || !(q >= 0)) {
                throw new ArgumentException($"{nameof(n)}, {nameof(q)}");
            }

            int n_sq = n * n;

            MultiPrecision<N> h = q * q;
            MultiPrecision<N> u = ((n & 1) == 0)
                ? (n == 0) ? 0d : 2d / (a + n_sq)
                : (n == 1) ? 1d / q : 1d / (a + (n_sq - 1) - q);
            MultiPrecision<N> v = 0d;

            for (int k = 2 + (n & 1); k < n; k += 2) {
                u = 1d / (a + (n_sq - k * k) - u * h);
            }
            for (int k = terms * 2 + n; k > n; k -= 2) {
                v = 1d / (a + (n_sq - k * k) - v * h);
            }

            MultiPrecision<N> y = (((n > 0) ? u : v) + v) * h - a;

            return y;
        }

        /// <summary>
        /// Continue fraction eigen B.
        /// </summary>
        /// <remarks>DLMF 28.6</remarks>
        private static MultiPrecision<N> FractionB(int n, MultiPrecision<N> q, MultiPrecision<N> a, int terms) {
            if (n < 1 || n > 16384 || !(q >= 0)) {
                throw new ArgumentException($"{nameof(n)}, {nameof(q)}");
            }

            int n_sq = n * n;

            MultiPrecision<N> h = q * q;
            MultiPrecision<N> u = ((n & 1) == 0)
                ? 0d
                : (n == 1) ? -1d / q : 1d / (a + (n_sq - 1) + q);
            MultiPrecision<N> v = 0d;

            for (int k = 2 + (n & 1); k < n; k += 2) {
                u = 1d / (a + (n_sq - k * k) - u * h);
            }
            for (int k = terms * 2 + n; k > n; k -= 2) {
                v = 1d / (a + (n_sq - k * k) - v * h);
            }

            MultiPrecision<N> y = (u + v) * h - a;

            return y;
        }

        /// <summary>
        /// If a matches the true value, return 0.
        /// </summary>
        public static MultiPrecision<N> Fraction(EigenFunc func, int n, MultiPrecision<N> q, MultiPrecision<N> a, int terms = -1) {
            terms = (terms < 0) ? FracTerms((double)q) : terms;

            MultiPrecision<N> y = func switch {
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
        public static (MultiPrecision<N> value, MultiPrecision<N> score, bool is_convergence) SearchRoot(
            EigenFunc func, int n, MultiPrecision<N> q, MultiPrecision<N> a,
            int frac_terms = -1, bool lowscore_interpolate = true) {

            if (q == 0) {
                return (0, 1, is_convergence: true);
            }

            frac_terms = (frac_terms < 0) ? FracTerms((double)q) : frac_terms;

            MultiPrecision<N> heuristics_err = MultiPrecision<N>.Max(1, MultiPrecision<N>.Abs(a)) * 1e-6;

            MultiPrecision<N> h = MultiPrecision<N>.Max(1, n * n) / 256d;
            MultiPrecision<N> truncation_thr = 2 + MultiPrecision<N>.Max(1, n * n) * 0.1;

            (MultiPrecision<N> ar, bool ar_convergence, MultiPrecision<N> ar_score) =
                RootFinderMP<N>.Search((a) => Fraction(func, n, q, a, frac_terms), a, h, truncation_thr);
            (MultiPrecision<N> ap, bool ap_convergence, _) =
                RootFinderMP<N>.Search((a) => 1 / Fraction(func, n, q, a, frac_terms), a, h, truncation_thr);

            (MultiPrecision<N> a_likelihood, MultiPrecision<N> score_likelihood) =
                (ar_convergence && MultiPrecision<N>.Abs(a - ar) < heuristics_err) ? (ar, ar_score) : (MultiPrecision<N>.NaN, 0);

            if (ap_convergence) {
                (MultiPrecision<N> apm, bool apm_convergence, MultiPrecision<N> apm_score) =
                    RootFinderMP<N>.Search((a) => Fraction(func, n, q, a, frac_terms), MultiPrecision<N>.BitDecrement(ap), h, truncation_thr, SearchDirection.Minus);
                (MultiPrecision<N> app, bool app_convergence, MultiPrecision<N> app_score) =
                    RootFinderMP<N>.Search((a) => Fraction(func, n, q, a, frac_terms), MultiPrecision<N>.BitIncrement(ap), h, truncation_thr, SearchDirection.Plus);

                if (apm_convergence && apm_score > ar_score * 0.5 && MultiPrecision<N>.Abs(a - apm) < heuristics_err) {
                    (a_likelihood, score_likelihood) =
                        MultiPrecision<N>.Abs(a - a_likelihood) < MultiPrecision<N>.Abs(a - apm) ? (a_likelihood, score_likelihood) : (apm, apm_score);
                }
                if (app_convergence && app_score > ar_score * 0.5 && MultiPrecision<N>.Abs(a - app) < heuristics_err) {
                    (a_likelihood, score_likelihood) =
                        MultiPrecision<N>.Abs(a - a_likelihood) < MultiPrecision<N>.Abs(a - app) ? (a_likelihood, score_likelihood) : (app, app_score);
                }
            }

            if (score_likelihood < 0.5 && lowscore_interpolate) {
                MultiPrecision<N> h_interpolate = MultiPrecision<N>.Ldexp(MultiPrecision<N>.Max(1, n * n), -MultiPrecision<N>.Bits / 2);
                MultiPrecision<N> q_m = MultiPrecision<N>.Min(q, h_interpolate), q_p = h_interpolate;
                (MultiPrecision<N> value_m, MultiPrecision<N> score_m, _) = SearchRoot(func, n, q - q_m, InitialValue(func, n, q - q_m), lowscore_interpolate: false);
                (MultiPrecision<N> value_p, MultiPrecision<N> score_p, _) = SearchRoot(func, n, q + q_p, InitialValue(func, n, q + q_p), lowscore_interpolate: false);

                while (score_m < 0.75) {
                    q_m = MultiPrecision<N>.Min(q, q_m + h_interpolate);
                    (value_m, score_m, _) = SearchRoot(func, n, q - q_m, InitialValue(func, n, q - q_m), lowscore_interpolate: false);
                }
                while (score_p < 0.75) {
                    q_p += h_interpolate;
                    (value_p, score_p, _) = SearchRoot(func, n, q + q_p, InitialValue(func, n, q + q_p), lowscore_interpolate: false);
                }

                MultiPrecision<N> a_interpolate = (value_p * q_m + value_m * q_p) / (q_m + q_p);

                return (a_interpolate, 0.5, is_convergence: false);
            }

            bool is_convergence = a_likelihood.IsFinite;
            if (!is_convergence && MultiPrecision<N>.Abs(a - ap) < heuristics_err) {
                a_likelihood = ap;
            }

            return (a_likelihood, score_likelihood, is_convergence);
        }

        /// <summary>
        /// Initial value for the secant method.
        /// </summary>
        public static MultiPrecision<N> InitialValue(EigenFunc func, int n, MultiPrecision<N> q) {
            double y = EigenFP64.Value(func, n, (double)q, zero_shift: true).value;

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
        public static (MultiPrecision<N> value, MultiPrecision<N> score, bool is_convergence) Value(
            EigenFunc func, int n, MultiPrecision<N> q,
            int frac_terms = -1, bool zero_shift = false, bool lowscore_interpolate = true) {

            frac_terms = (frac_terms < 0) ? FracTerms((double)q) : frac_terms;

            MultiPrecision<N> a0 = InitialValue(func, n, q);
            (MultiPrecision<N> value, MultiPrecision<N> score, bool is_convergence) = SearchRoot(func, n, q, a0, frac_terms, lowscore_interpolate);

            if (!zero_shift) {
                value += n * n;
            }

            return (value, score, is_convergence);
        }

        public static (MultiPrecision<N> value, int terms) ConvergenceFracTerms(EigenFunc func, int n, MultiPrecision<N> q, int init_terms) {
            MultiPrecision<N> a0 = Value(func, n, q, frac_terms: init_terms, zero_shift: true, lowscore_interpolate: false).value;
            MultiPrecision<N> a1 = Value(func, n, q, frac_terms: init_terms + 1, zero_shift: true, lowscore_interpolate: false).value;
            MultiPrecision<N> err01 = MultiPrecision<N>.Abs(a0 - a1);

            for (int frac_terms = init_terms; frac_terms <= 4096; frac_terms++) {
                MultiPrecision<N> a2 = Value(func, n, q, frac_terms + 2, zero_shift: true, lowscore_interpolate: false).value;
                MultiPrecision<N> err12 = MultiPrecision<N>.Abs(a1 - a2);

                if ((err01.Exponent < a0.Exponent - MultiPrecision<N>.Bits + 6) && (err01 <= err12) && ((a0 >= a1 && a1 <= a2) || (a0 <= a1 && a1 >= a2))) {
                    return (a2, frac_terms);
                }

                (a0, a1, err01) = (a1, a2, err12);
            }

            return (MultiPrecision<N>.NaN, -1);
        }

        public static int FracTerms(double q) {
            if (!(q >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(q));
            }

            int mp_length = MultiPrecision<N>.Length;

            double logq = Math.Log2(q);
            
            if (mp_length <= 4) {
                return (q < 1e-3)
                    ? 6
                    : (int)Math.Ceiling(Math.Pow(2, 3.3479 + 1.6171e-1 * logq + 6.4659e-3 * logq * logq));
            }
            if (mp_length <= 8) {
                return (q < 1e-3)
                    ? 10
                    : (int)Math.Ceiling(Math.Pow(2, 4.0646 + 1.3941e-1 * logq + 5.9064e-3 * logq * logq));
            }
            if (mp_length <= 16) {
                return (q < 1e-3)
                    ? 18
                    : (int)Math.Ceiling(Math.Pow(2, 4.8178 + 1.2211e-1 * logq + 5.2052e-3 * logq * logq));
            }
            if (mp_length <= 32) {
                return (q < 1e-3)
                    ? 29
                    : (int)Math.Ceiling(Math.Pow(2, 5.6300 + 1.1971e-1 * logq + 4.0530e-3 * logq * logq));
            }
            if (mp_length <= 64) {
                return (q < 1e-3)
                    ? 50
                    : (int)Math.Ceiling(Math.Pow(2, 6.4471 + 1.0543e-1 * logq + 2.4224e-3 * logq * logq));
            }
            if (mp_length <= 128) {
                return (q < 1e-3)
                    ? 94
                    : (int)Math.Ceiling(Math.Pow(2, 7.3228 + 9.3942e-2 * logq + 1.5475e-3 * logq * logq));
            }

            throw new NotImplementedException();
        }
    }
}
