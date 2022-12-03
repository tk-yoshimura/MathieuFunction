using System.Xml.Linq;

namespace MathieuMP {
    public static class RootFinder {
        const double eps = 4e-15;
        const double truncation_thr = 256; 

        public static (double v, bool is_convergence) BothSearch(Func<double, double> f, double x0, double h0) {
            if (!(h0 >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h0));
            }

            double h = h0, x = x0, y = f(x0);
            bool is_convergenced = false, is_clamp;

            while (Math.Abs(h / (Math.Abs(x) + double.Epsilon)) >= eps) {
                (double yn3, double yn2, double yn1, double yp1, double yp2, double yp3) =
                    (f(x - h * 3), f(x - h * 2), f(x - h), f(x + h), f(x + h * 2), f(x + h * 3));

                if (SequenceUtil.IsReciprocalCurve(yn3, yn2, yn1, y, yp1, yp2, yp3)) {
                    double dx = (yp1 == yn1) ? 0 : 2 * h / (yp1 - yn1) * y;
                    double dx_clamp = Math.Max(-h * 32, Math.Min(h * 32, dx));

                    (dx, is_clamp) = (dx_clamp, dx != dx_clamp);
                    double x_next = x - dx, y_next = f(x_next);

                    if (yn3 * yp3 <= 0 || SequenceUtil.IsMonotone(y_next, f(x - dx * 3 / 4), f(x - dx / 2), f(x - dx / 4), y)) {
                        (x, y) = (x_next, y_next);

                        if (is_clamp) {
                            h = Math.Min(h * 2, h0);
                        }
                        else if (yn1 * yp1 <= 0) {
                            h /= 4;
                        }

                        if (Math.Abs(dx / (Math.Abs(x) + double.Epsilon)) <= eps) {
                            if (y == 0 || is_convergenced) {
                                is_convergenced = true;
                                break;
                            }

                            is_convergenced = true;
                        }

                        if (Math.Abs(x0 - x) > truncation_thr) {
                            return (double.NaN, false);
                        }

                        continue;
                    }
                }

                h /= 4;
            }

            return (x, is_convergenced);
        }

        public static (double v, bool is_convergence) MinusSearch(Func<double, double> f, double x0, double h0) {
            if (!(h0 >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h0));
            }

            double s = Math.Max(double.Epsilon, x0 - Math.BitDecrement(x0));
            double h = s, x, y;
            int sign = NeighborSign(f, x0, SearchDirection.Minus);

            do {
                x = x0 - s;
                y = f(x);

                if (double.IsFinite(y) && Math.Sign(y) == sign) {
                    break;
                }
                s *= 2;
            } while (s < h0);

            if (s >= h0) { 
                return (double.NaN, false);
            }

            while (h < double.PositiveInfinity) {
                if (y * f(x - h * 12) <= 0) {
                    x -= h * 3;
                    y = f(x);
                    break;
                }

                h = Math.Max(Math.BitIncrement(h), h * 1.25);
            }

            if (!double.IsFinite(x) || !double.IsFinite(h)) {
                return (double.NaN, false);
            }
            if (Math.Abs(h / (Math.Abs(x) + double.Epsilon)) < eps) {
                return (x, true);
            }

            h = Math.Min(h, h0);

            bool is_convergenced = false, is_clamp;

            while (Math.Abs(h / (Math.Abs(x) + double.Epsilon)) >= eps) {
                (double yn3, double yn2, double yn1, double yp1, double yp2, double yp3) =
                    (f(x - h * 3), f(x - h * 2), f(x - h), f(x + h), f(x + h * 2), f(x + h * 3));

                if (SequenceUtil.IsReciprocalCurve(yn3, yn2, yn1, y, yp1, yp2, yp3)) {
                    double dx = (yp1 == yn1) ? 0 : 2 * h / (yp1 - yn1) * y;
                    double dx_clamp = Math.Max(-h * 32, Math.Min(h * 32, dx));

                    (dx, is_clamp) = (dx_clamp, dx != dx_clamp);
                    double x_next = Math.Min(x0, x - dx), y_next = f(x_next);

                    if (yn3 * yp3 <= 0 || SequenceUtil.IsMonotone(y_next, f(x - dx * 3 / 4), f(x - dx / 2), f(x - dx / 4), y)) {
                        (x, y) = (x_next, y_next);

                        if (is_clamp) {
                            h = Math.Min(h * 2, h0);
                        }
                        else if (yn1 * yp1 <= 0) {
                            h /= 4;
                        }

                        if (Math.Abs(dx / (Math.Abs(x) + double.Epsilon)) <= eps) {
                            if (y == 0 || is_convergenced) {
                                is_convergenced = true;
                                break;
                            }

                            is_convergenced = true;
                        }

                        if (Math.Abs(x0 - x) > truncation_thr) {
                            return (double.NaN, false);
                        }

                        continue;
                    }
                }

                h /= 4;
            }

            return (x, is_convergenced);
        }

        public static (double v, bool is_convergence) PlusSearch(Func<double, double> f, double x0, double h0) {
            if (!(h0 >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h0));
            }

            double s = Math.Max(double.Epsilon, Math.BitIncrement(x0) - x0);
            double h = s, x, y;
            int sign = NeighborSign(f, x0, SearchDirection.Plus);

            do {
                x = x0 + s;
                y = f(x);

                if (double.IsFinite(y) && Math.Sign(y) == sign) {
                    break;
                }
                s *= 2;
            } while (s < h0);

            if (s >= h0) { 
                return (double.NaN, false);
            }

            while (h < double.PositiveInfinity) {
                if (y * f(x + h * 12) <= 0) {
                    x += h * 3;
                    y = f(x);
                    break;
                }

                h = Math.Max(Math.BitIncrement(h), h * 1.25);
            }

            if (!double.IsFinite(x) || !double.IsFinite(h)) {
                return (double.NaN, false);
            }
            if (Math.Abs(h / (Math.Abs(x) + double.Epsilon)) < eps) {
                return (x, true);
            }

            h = Math.Min(h, h0);

            bool is_convergenced = false, is_clamp;

            while (Math.Abs(h / (Math.Abs(x) + double.Epsilon)) >= eps) {
                (double yn3, double yn2, double yn1, double yp1, double yp2, double yp3) =
                    (f(x - h * 3), f(x - h * 2), f(x - h), f(x + h), f(x + h * 2), f(x + h * 3));

                if (SequenceUtil.IsReciprocalCurve(yn3, yn2, yn1, y, yp1, yp2, yp3)) {
                    double dx = (yp1 == yn1) ? 0 : 2 * h / (yp1 - yn1) * y;
                    double dx_clamp = Math.Max(-h * 32, Math.Min(h * 32, dx));

                    (dx, is_clamp) = (dx_clamp, dx != dx_clamp);
                    double x_next = Math.Max(x0, x - dx), y_next = f(x_next);

                    if (yn3 * yp3 <= 0 || SequenceUtil.IsMonotone(y_next, f(x - dx * 3 / 4), f(x - dx / 2), f(x - dx / 4), y)) {
                        (x, y) = (x_next, y_next);

                        if (is_clamp) {
                            h = Math.Min(h * 2, h0);
                        }
                        else if (yn1 * yp1 <= 0) {
                            h /= 4;
                        }

                        if (Math.Abs(dx / (Math.Abs(x) + double.Epsilon)) <= eps) {
                            if (y == 0 || is_convergenced) {
                                is_convergenced = true;
                                break;
                            }

                            is_convergenced = true;
                        }

                        if (Math.Abs(x0 - x) > truncation_thr) {
                            return (double.NaN, false);
                        }

                        continue;
                    }
                }

                h /= 4;
            }

            return (x, is_convergenced);
        }

        public static int NeighborSign(Func<double, double> f, double x, SearchDirection direction) {
            if (direction != SearchDirection.Minus && direction != SearchDirection.Plus) {
                throw new ArgumentException(nameof(direction));
            }

            double h = ((direction == SearchDirection.Minus) ? Math.BitDecrement(x) : Math.BitIncrement(x)) - x;

            if (!double.IsFinite(h) || h == 0) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }

            int signw = 0;

            for (int i = 1, w = 8; i <= 8; i++, w--) {
                double y = f(x + h);
                h *= 2;

                if (!double.IsFinite(y)) {
                    continue;
                }

                signw += w * Math.Sign(y);
            }

            int sign = (signw == 0) ? 0 : (signw > 0) ? +1 : -1;

            return sign;
        }

        public static (double v, bool is_convergence) Search(Func<double, double> f, double x, double h, SearchDirection direction = SearchDirection.Both) {
            return direction switch {
                SearchDirection.Minus => MinusSearch(f, x, h),
                SearchDirection.Plus => PlusSearch(f, x, h),
                _ => BothSearch(f, x, h),
            };
        }
    }
}
