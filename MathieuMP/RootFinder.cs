namespace MathieuMP {
    public static class RootFinder {
        const double eps = 4e-15;
        const int max_secant_iter = 32;

        private static (double v, bool is_convergence) BothSearch(Func<double, double> f, double x, double h) {
            if (!(h >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h));
            }

            double x0 = x, y0 = f(x0);
            bool is_neighbor_monotone = true;
 
            for (int secant_iter = 0; secant_iter < max_secant_iter; secant_iter++) {
                (double yn2, double yn1, double yp1, double yp2) = (f(x0 - h * 2), f(x0 - h), f(x0 + h), f(x0 + h * 2));

                if (!SequenceUtil.IsMonotone(yn2, yn1, y0, yp1, yp2)) {
                    is_neighbor_monotone = false;
                    h /= 4;

                    if (x0 + h == x0 || x0 - h == x0) {
                        break;
                    }

                    continue;
                }

                double dx = (yp1 == yn1) ? 0 : 2 * h / (yp1 - yn1) * y0;
                if (!is_neighbor_monotone) {
                    dx = Math.Max(-h, Math.Min(h, dx));
                }

                x0 -= dx;
                y0 = f(x0);

                if (y0 == 0 || Math.Abs(dx / (Math.Abs(x0) + double.Epsilon)) <= eps || f(Math.BitDecrement(x0)) * f(Math.BitIncrement(x0)) <= 0) {
                    return (x0, true);
                }
            }

            return (x0, false);
        }

        private static (double v, bool is_convergence) MinusSearch(Func<double, double> f, double x, double h) {
            if (!(h >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h));
            }

            double x0 = x, y0 = f(x0);
            bool is_neighbor_monotone = true;

            for (int secant_iter = 0; secant_iter < max_secant_iter; secant_iter++) {
                (double yn4, double yn3, double yn2, double yn1) = (f(x0 - h * 4), f(x0 - h * 3), f(x0 - h * 2), f(x0 - h));

                if (!SequenceUtil.IsMonotone(yn4, yn3, yn2, yn1, y0)) {
                    is_neighbor_monotone = false;
                    h /= 4;

                    if (x0 - h == x0) {
                        break;
                    }

                    continue;
                }

                double dx = (y0 == yn2) ? 0 : 2 * h / (y0 - yn2) * y0;
                if (!is_neighbor_monotone) {
                    dx = Math.Max(-h, Math.Min(h, dx));
                }

                x0 = Math.Min(x, x0 - dx);
                y0 = f(x0);

                if (y0 == 0 || Math.Abs(dx / (Math.Abs(x0) + double.Epsilon)) <= eps) {
                    return (x0, true);
                }
            }

            return (x0, false);
        }

        private static (double v, bool is_convergence) PlusSearch(Func<double, double> f, double x, double h) {
            if (!(h >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h));
            }

            double x0 = x, y0 = f(x0);
            bool is_neighbor_monotone = true;

            for (int secant_iter = 0; secant_iter < max_secant_iter; secant_iter++) {
                (double yp1, double yp2, double yp3, double yp4) = (f(x0 + h), f(x0 + h * 2), f(x0 + h * 3), f(x0 + h * 4));

                if (!SequenceUtil.IsMonotone(y0, yp1, yp2, yp3, yp4)) {
                    is_neighbor_monotone = false;
                    h /= 4;

                    if (x0 - h == x0) {
                        break;
                    }

                    continue;
                }

                double dx = (yp2 == y0) ? 0 : 2 * h / (yp2 - y0) * y0;
                if (!is_neighbor_monotone) {
                    dx = Math.Max(-h, Math.Min(h, dx));
                }

                x0 = Math.Max(x, x0 - dx);
                y0 = f(x0);

                if (y0 == 0 || Math.Abs(dx / (Math.Abs(x0) + double.Epsilon)) <= eps) {
                    return (x0, true);
                }
            }

            return (x0, false);
        }

        public static (double v, bool is_convergence) Search(Func<double, double> f, double x, double h, SearchDirection direction = SearchDirection.Both) { 
            return direction switch {
                SearchDirection.Minus => MinusSearch(f, x, h),
                SearchDirection.Plus  => PlusSearch(f, x, h),
                _ => BothSearch(f, x, h),
            };
        }
    }
}
