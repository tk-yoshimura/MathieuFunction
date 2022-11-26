namespace MathieuMP {
    public static class RootFinder {
        const double eps = 1e-15;
        const int max_secant_iter = 8;

        public static double StepSearch(Func<double, double> f, double x, double h, int max_counts = 1024) {
            if (!(h >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h));
            }

            double y = f(x), xn = double.NaN, xp = double.NaN;
            double xn0, xn1, yn0, yn1, xp0, xp1, yp0, yp1;

            (xn0, xp0) = (x, x);
            (yn0, yp0) = (y, y);

            for (int i = 1; i <= max_counts && !double.IsFinite(xn) && !double.IsFinite(xp); i++) {
                (xn1, xp1) = (xn0 - h, xp0 + h);
                (yn1, yp1) = (f(xn1), f(xp1));

                if (yn0 * yn1 <= 0) {
                    double yh = f(xn0 - h / 2);
                    if ((yn0 <= yh && yh <= yn1) || (yn0 >= yh && yh >= yn1)) { 
                        xn = (yn0 == yn1) ? xn0 : (xn0 + h / (yn1 - yn0) * yn0);
                    }
                }
                if (yp0 * yp1 <= 0) {
                    double yh = f(xp0 + h / 2);
                    if ((yp0 <= yh && yh <= yp1) || (yp0 >= yh && yh >= yp1)) {
                        xp = (yp0 == yp1) ? xp0 : (xp0 - h / (yp1 - yp0) * yp0);
                    }
                }

                (xn0, xp0) = (xn1, xp1);
                (yn0, yp0) = (yn1, yp1);
            }

            double xr = (double.IsNaN(xp) || Math.Abs(x - xn) < Math.Abs(x - xp)) ? xn : xp;

            return xr;
        }

        public static (double v, bool is_convergence) SecantSearch(Func<double, double> f, double x, double h) {
            if (!(h >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h));
            }

            double x0 = StepSearch(f, x, h);
            double yn2, yn1, y0, yp1, yp2;

            if (double.IsNaN(x0)) {
                return (x, is_convergence: false);
            }

            bool is_convergence = false;
            int secant_iter = 0;
            (yn2, y0, yp2) = (f(x0 - h), f(x0), f(x0 + h));

            h /= 2;
            while ((yn2 * yp2 <= 0) && (h / (Math.Abs(x0) + double.Epsilon) >= eps)) {
                (yn1, yp1) = (f(x0 - h), f(x0 + h));

                if (SequenceUtil.IsMonotone(yn2, yn1, y0, yp1, yp2)) {
                    if (yn1 * yp1 <= 0) {
                        double dx = (yp1 == yn1) ? 0 : 2 * h / (yp1 - yn1) * y0;
                        dx = Math.Max(-h, Math.Min(h, dx));

                        x0 -= dx;
                        secant_iter++;

                        if ((secant_iter >= max_secant_iter) || (Math.Abs(dx / (Math.Abs(x0) + double.Epsilon)) <= eps)) {
                            is_convergence = true;
                            break;
                        }

                        y0 = f(x0);
                        if (y0 == 0) {
                            is_convergence = true;
                            break;
                        }

                        (yn2, yp2) = (f(x0 - h), f(x0 + h));
                        h /= 2;

                        continue;
                    }
                    if (yn2 * yn1 <= 0) {
                        x0 -= h;
                        (yn2, y0, yp2) = (f(x0 - h * 2), yn1, yp1);

                        continue;
                    }
                    if (yp2 * yp1 <= 0) {
                        x0 += h;
                        (yn2, y0, yp2) = (yn1, yp1, f(x0 + h * 2));

                        continue;
                    }
                }

                (yn2, yp2) = (yn1, yp1);
                h /= 2;
            }

            return (x0, is_convergence);
        }
    }
}
