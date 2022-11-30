namespace MathieuMP {
    public static class RootFinder {
        const double eps = 4e-15;
        const int max_secant_iter = 32;

        private static (double v, bool is_convergence) BothSearch(Func<double, double> f, double x, double h) {
            if (!(h >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h));
            }

            double x0 = x, y0 = f(x0);
            bool is_neighbor_monotone = true, is_convergenced = false;

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

                double dx = (yp2 == yn2) ? 0 : 4 * h / (yp2 - yn2) * y0;
                dx = Math.Max(-h * 2, Math.Min(h * 2, dx));

                x0 -= dx;
                y0 = f(x0);

                if (Math.Abs(dx / (Math.Abs(x0) + double.Epsilon)) <= eps) {
                    if (y0 == 0 || is_convergenced) {
                        return (x0, true);
                    }

                    is_convergenced = true;
                }
            }

            return (x0, false);
        }

        private static (double v, bool is_convergence) MinusSearch(Func<double, double> f, double x, double h) {
            if (!(h >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h));
            }

            double x0 = x, y0 = f(x0);
            bool is_neighbor_monotone = true, is_convergenced = false;

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

                double dx = (y0 == yn4) ? 0 : 4 * h / (y0 - yn4) * y0;
                dx = Math.Max(-h * 2, Math.Min(h * 2, dx));

                if (!is_neighbor_monotone) {
                    dx = Math.Max(-h * 2, Math.Min(h * 2, dx));
                }

                x0 = Math.Min(x, x0 - dx);
                y0 = f(x0);

                if (Math.Abs(dx / (Math.Abs(x0) + double.Epsilon)) <= eps) {
                    if (y0 == 0 || is_convergenced) {
                        return (x0, true);
                    }

                    is_convergenced = true;
                }
                if (x0 == x) {
                    h /= 4;
                }
            }

            return (x0, false);
        }

        private static (double v, bool is_convergence) PlusSearch(Func<double, double> f, double x, double h) {
            if (!(h >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h));
            }

            double x0 = x, y0 = f(x0);
            bool is_neighbor_monotone = true, is_convergenced = false;

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

                double dx = (yp4 == y0) ? 0 : 4 * h / (yp4 - y0) * y0;
                dx = Math.Max(-h * 2, Math.Min(h * 2, dx));

                x0 = Math.Max(x, x0 - dx);
                y0 = f(x0);

                if (Math.Abs(dx / (Math.Abs(x0) + double.Epsilon)) <= eps) {
                    if (y0 == 0 || is_convergenced) {
                        return (x0, true);
                    }

                    is_convergenced = true;
                }
                if (x0 == x) {
                    h /= 4;
                }
            }

            return (x0, false);
        }

        public static (double v, bool is_convergence) AdvancedSearch(Func<double, double> f, double x0, double h0) {
            if (!(h0 >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h0));
            }

            double h = h0, x = x0, y = f(x0);
            bool is_convergenced = false, is_clamp;

            while (Math.Abs(h / (Math.Abs(x) + double.Epsilon)) >= eps) {
                (double yn3, double yn2, double yn1, double yp1, double yp2, double yp3) = 
                    (f(x - h * 3), f(x - h * 2), f(x - h), f(x + h), f(x + h * 2), f(x + h * 3));

                if (SequenceUtil.IsReciprocalCurve(yn3, yn2, yn1, y, yp1, yp2, yp3)) {
                    double dx = (yp2 == yn2) ? 0 : 2 * h / (yp1 - yn1) * y;
                    double dx_clamp = Math.Max(-h * 32, Math.Min(h * 32, dx));

                    (dx, is_clamp) = (dx_clamp, dx != dx_clamp);
                    double x_next = x - dx, y_next = f(x_next);

                    if (yn3 * yp3 <= 0 || SequenceUtil.IsMonotone(y_next, f(x - dx * 3 / 4), f(x - dx / 2), f(x - dx / 4), y)){
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

                        continue;
                    }
                }

                h /= 4;
            }

            return (x, is_convergenced);
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
