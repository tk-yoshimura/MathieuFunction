using MultiPrecision;

namespace MathieuMP {
    public static class RootFinderMP<N> where N : struct, IConstant {
        static readonly int max_secant_iters = MultiPrecision<N>.Bits / 3;

        public static (MultiPrecision<N> v, bool is_convergence) BothSearch(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> x0, MultiPrecision<N> h0, MultiPrecision<N> truncation_thr) {
            if (!(h0 >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h0));
            }

            MultiPrecision<N> h = h0, x = x0, y = f(x0);
            bool is_convergenced = false, is_clamp;
            int secant_iters = 0;

            while (h.Exponent >= x.Exponent - MultiPrecision<N>.Bits + 4) {
                (MultiPrecision<N> yn3, MultiPrecision<N> yn2, MultiPrecision<N> yn1, MultiPrecision<N> yp1, MultiPrecision<N> yp2, MultiPrecision<N> yp3) =
                    (f(x - h * 3), f(x - h * 2), f(x - h), f(x + h), f(x + h * 2), f(x + h * 3));

                if (SequenceUtilMP<N>.IsReciprocalCurve(yn3, yn2, yn1, y, yp1, yp2, yp3)) {
                    MultiPrecision<N> dx = (yp1 == yn1) ? 0 : 2 * h / (yp1 - yn1) * y;
                    MultiPrecision<N> dx_clamp = MultiPrecision<N>.Max(-h * 32, MultiPrecision<N>.Min(h * 32, dx));

                    (dx, is_clamp) = (dx_clamp, dx != dx_clamp);
                    MultiPrecision<N> x_next = x - dx, y_next = f(x_next);

                    if (yn3 * yp3 <= 0 || SequenceUtilMP<N>.IsMonotone(y_next, f(x - dx * 3 / 4), f(x - dx / 2), f(x - dx / 4), y)) {
                        (x, y) = (x_next, y_next);

                        if (is_clamp) {
                            h = MultiPrecision<N>.Min(h * 2.5, h0);
                        }
                        else if (yn1 * yp1 <= 0) {
                            secant_iters++;
                            h /= 4;
                        }

                        if (dx.Exponent <= x.Exponent - MultiPrecision<N>.Bits + 4 || secant_iters >= max_secant_iters) {
                            if (y == 0 || is_convergenced) {
                                is_convergenced = true;
                                break;
                            }

                            is_convergenced = true;
                        }

                        if (MultiPrecision<N>.Abs(x0 - x) > truncation_thr) {
                            return (MultiPrecision<N>.NaN, false);
                        }

                        continue;
                    }
                }

                h /= 4;
            }

            return (x, is_convergenced);
        }

        public static (MultiPrecision<N> v, bool is_convergence) MinusSearch(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> x0, MultiPrecision<N> h0, MultiPrecision<N> truncation_thr) {
            if (!(h0 >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h0));
            }

            MultiPrecision<N> s = MultiPrecision<N>.Max(MultiPrecision<N>.Epsilon, x0 - MultiPrecision<N>.BitDecrement(x0));
            MultiPrecision<N> h = s, x, y;

            do {
                x = x0 - s;
                y = f(x);

                if (MultiPrecision<N>.IsFinite(y) && y < 0) {
                    break;
                }
                s *= 2;
            } while (s < h0);

            if (s >= h0) {
                return (MultiPrecision<N>.NaN, false);
            }

            while (h < MultiPrecision<N>.PositiveInfinity) {
                if (y * f(x - h * 12) <= 0) {
                    x -= h * 3;
                    y = f(x);
                    break;
                }

                h = MultiPrecision<N>.Max(MultiPrecision<N>.BitIncrement(h), h * 1.25);
            }

            if (!MultiPrecision<N>.IsFinite(x) || !MultiPrecision<N>.IsFinite(h)) {
                return (MultiPrecision<N>.NaN, false);
            }
            if (h.Exponent < x.Exponent - MultiPrecision<N>.Bits + 4) {
                return (x, true);
            }

            h = MultiPrecision<N>.Min(h, h0);

            bool is_convergenced = false, is_clamp;
            int secant_iters = 0;

            while (h.Exponent >= x.Exponent - MultiPrecision<N>.Bits + 4) {
                (MultiPrecision<N> yn3, MultiPrecision<N> yn2, MultiPrecision<N> yn1, MultiPrecision<N> yp1, MultiPrecision<N> yp2, MultiPrecision<N> yp3) =
                    (f(x - h * 3), f(x - h * 2), f(x - h), f(x + h), f(x + h * 2), f(x + h * 3));

                if (SequenceUtilMP<N>.IsReciprocalCurve(yn3, yn2, yn1, y, yp1, yp2, yp3)) {
                    MultiPrecision<N> dx = (yp1 == yn1) ? 0 : 2 * h / (yp1 - yn1) * y;
                    MultiPrecision<N> dx_clamp = MultiPrecision<N>.Max(-h * 32, MultiPrecision<N>.Min(h * 32, dx));

                    (dx, is_clamp) = (dx_clamp, dx != dx_clamp);
                    MultiPrecision<N> x_next = MultiPrecision<N>.Min(x0, x - dx), y_next = f(x_next);

                    if (yn3 * yp3 <= 0 || SequenceUtilMP<N>.IsMonotone(y_next, f(x - dx * 3 / 4), f(x - dx / 2), f(x - dx / 4), y)) {
                        (x, y) = (x_next, y_next);

                        if (is_clamp) {
                            h = MultiPrecision<N>.Min(h * 2.5, h0);
                        }
                        else if (yn1 * yp1 <= 0) {
                            secant_iters++;
                            h /= 4;
                        }

                        if (dx.Exponent <= x.Exponent - MultiPrecision<N>.Bits + 4 || secant_iters >= max_secant_iters) {
                            if (y == 0 || is_convergenced) {
                                is_convergenced = true;
                                break;
                            }

                            is_convergenced = true;
                        }

                        if (MultiPrecision<N>.Abs(x0 - x) > truncation_thr) {
                            return (MultiPrecision<N>.NaN, false);
                        }

                        continue;
                    }
                }

                h /= 4;
            }

            return (x, is_convergenced);
        }

        public static (MultiPrecision<N> v, bool is_convergence) PlusSearch(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> x0, MultiPrecision<N> h0, MultiPrecision<N> truncation_thr) {
            if (!(h0 >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h0));
            }

            MultiPrecision<N> s = MultiPrecision<N>.Max(MultiPrecision<N>.Epsilon, MultiPrecision<N>.BitIncrement(x0) - x0);
            MultiPrecision<N> h = s, x, y;

            do {
                x = x0 + s;
                y = f(x);

                if (MultiPrecision<N>.IsFinite(y) && y > 0) {
                    break;
                }
                s *= 2;
            } while (s < h0);

            if (s >= h0) {
                return (MultiPrecision<N>.NaN, false);
            }

            while (h < MultiPrecision<N>.PositiveInfinity) {
                if (y * f(x + h * 12) <= 0) {
                    x += h * 3;
                    y = f(x);
                    break;
                }

                h = MultiPrecision<N>.Max(MultiPrecision<N>.BitIncrement(h), h * 1.25);
            }

            if (!MultiPrecision<N>.IsFinite(x) || !MultiPrecision<N>.IsFinite(h)) {
                return (MultiPrecision<N>.NaN, false);
            }
            if (h.Exponent < x.Exponent - MultiPrecision<N>.Bits + 4) {
                return (x, true);
            }

            h = MultiPrecision<N>.Min(h, h0);

            bool is_convergenced = false, is_clamp;
            int secant_iters = 0;

            while (h.Exponent >= x.Exponent - MultiPrecision<N>.Bits + 4) {
                (MultiPrecision<N> yn3, MultiPrecision<N> yn2, MultiPrecision<N> yn1, MultiPrecision<N> yp1, MultiPrecision<N> yp2, MultiPrecision<N> yp3) =
                    (f(x - h * 3), f(x - h * 2), f(x - h), f(x + h), f(x + h * 2), f(x + h * 3));

                if (SequenceUtilMP<N>.IsReciprocalCurve(yn3, yn2, yn1, y, yp1, yp2, yp3)) {
                    MultiPrecision<N> dx = (yp1 == yn1) ? 0 : 2 * h / (yp1 - yn1) * y;
                    MultiPrecision<N> dx_clamp = MultiPrecision<N>.Max(-h * 32, MultiPrecision<N>.Min(h * 32, dx));

                    (dx, is_clamp) = (dx_clamp, dx != dx_clamp);
                    MultiPrecision<N> x_next = MultiPrecision<N>.Max(x0, x - dx), y_next = f(x_next);

                    if (yn3 * yp3 <= 0 || SequenceUtilMP<N>.IsMonotone(y_next, f(x - dx * 3 / 4), f(x - dx / 2), f(x - dx / 4), y)) {
                        (x, y) = (x_next, y_next);

                        if (is_clamp) {
                            h = MultiPrecision<N>.Min(h * 2.5, h0);
                        }
                        else if (yn1 * yp1 <= 0) {
                            secant_iters++;
                            h /= 4;
                        }

                        if (dx.Exponent <= x.Exponent - MultiPrecision<N>.Bits + 4 || secant_iters >= max_secant_iters) {
                            if (y == 0 || is_convergenced) {
                                is_convergenced = true;
                                break;
                            }

                            is_convergenced = true;
                        }

                        if (MultiPrecision<N>.Abs(x0 - x) > truncation_thr) {
                            return (MultiPrecision<N>.NaN, false);
                        }

                        continue;
                    }
                }

                h /= 4;
            }

            return (x, is_convergenced);
        }

        public static MultiPrecision<N> LinearityScore(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> x) {
            MultiPrecision<N> h0 = MultiPrecision<N>.Max(MultiPrecision<N>.Ldexp(1, -256), MultiPrecision<N>.BitIncrement(x) - x);

            if (!MultiPrecision<N>.IsFinite(h0) || h0 == 0) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }

            MultiPrecision<N> score = 0;
            MultiPrecision<N> ny = f(x), iy = 1d / ny;

            for (int n = 16; n >= 2; n--) {
                MultiPrecision<N> h = h0;
                MultiPrecision<N> sw = 0, snw = 0, siw = 0, snxx = 0, snxy = 0, snyy = 0, sixx = 0, sixy = 0, siyy = 0;

                for (int i = 1, w = 2 << n; i <= n; i++, w /= 2) {
                    MultiPrecision<N> my = f(x - h), py = f(x + h);

                    MultiPrecision<N> mny = my - ny, pny = py - ny;
                    MultiPrecision<N> miy = 1d / my - iy, piy = 1d / py - iy;

                    sw += w;

                    if (MultiPrecision<N>.IsFinite(mny) && MultiPrecision<N>.IsFinite(pny)) {
                        snw += w;
                        snxx += 2 * w * h * h;
                        snxy += w * h * (pny - mny);
                        snyy += w * (mny * mny + pny * pny);
                    }

                    if (MultiPrecision<N>.IsFinite(miy) && MultiPrecision<N>.IsFinite(piy)) {
                        siw += w;
                        sixx += 2 * w * h * h;
                        sixy += w * h * (piy - miy);
                        siyy += w * (miy * miy + piy * piy);
                    }

                    h *= 2;
                }

                MultiPrecision<N> snxypw = snxy / snw, sixypw = sixy / siw;

                MultiPrecision<N> r = MultiPrecision<N>.IsFinite(snxypw) && MultiPrecision<N>.IsFinite(sixypw)
                    ? MultiPrecision<N>.Min(MultiPrecision<N>.Abs(snxypw), MultiPrecision<N>.Abs(sixypw))
                    : MultiPrecision<N>.IsFinite(snxypw) ? MultiPrecision<N>.Abs(snxypw)
                    : MultiPrecision<N>.IsFinite(sixypw) ? MultiPrecision<N>.Abs(sixypw)
                    : 0d;

                MultiPrecision<N> rn = MultiPrecision<N>.Abs(snxy / MultiPrecision<N>.Max(MultiPrecision<N>.Epsilon, MultiPrecision<N>.Sqrt(snxx * snyy))) * MultiPrecision<N>.Exp((r - MultiPrecision<N>.Abs(snxypw)) * sw);
                MultiPrecision<N> ri = MultiPrecision<N>.Abs(sixy / MultiPrecision<N>.Max(MultiPrecision<N>.Epsilon, MultiPrecision<N>.Sqrt(sixx * siyy))) * MultiPrecision<N>.Exp((r - MultiPrecision<N>.Abs(sixypw)) * sw);

                score = (MultiPrecision<N>.IsFinite(rn) ? rn : 0) - (MultiPrecision<N>.IsFinite(ri) ? ri : 0);

                if (MultiPrecision<N>.Abs(score) > 0.5) {
                    return score;
                }
            }

            return score;
        }

        public static (MultiPrecision<N> v, bool is_convergence, MultiPrecision<N> score) Search(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> x, MultiPrecision<N> h, MultiPrecision<N> truncation_thr, SearchDirection direction = SearchDirection.Both) {
            (MultiPrecision<N> v, bool is_convergence) = direction switch {
                SearchDirection.Minus => MinusSearch(f, x, h, truncation_thr),
                SearchDirection.Plus => PlusSearch(f, x, h, truncation_thr),
                _ => BothSearch(f, x, h, truncation_thr),
            };

            MultiPrecision<N> score = is_convergence ? LinearityScore(f, v) : -1;

            if (score < 0d) {
                is_convergence = false;
            }

            return (v, is_convergence, score);
        }
    }
}
