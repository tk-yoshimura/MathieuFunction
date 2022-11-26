namespace MathieuMP {
    public static class RootFinder {
        public static double StepSearch(Func<double, double> f, double x, double h, int max_counts = 16) {
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
                    xn = (yn0 == yn1) ? xn0 : (xn0 + h / (yn1 - yn0) * yn0);
                }
                if (yp0 * yp1 <= 0) {
                    xp = (yp0 == yp1) ? xp0 : (xp0 - h / (yp1 - yp0) * yp0);
                }

                (xn0, xp0) = (xn1, xp1);
                (yn0, yp0) = (yn1, yp1);
            }

            double xr = (double.IsNaN(xp) || Math.Abs(x - xn) < Math.Abs(x - xp)) ? xn : xp;

            Console.WriteLine(xr);

            return xr;
        }

        public static (double v, double error) SecantSearch(Func<double, double> f, double x, double h) {
            if (!(h >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(h));
            }

            double x0 = StepSearch(f, x, h), dx = double.NaN;
            double yn2, yn1, y0, yp1, yp2;

            if (double.IsNaN(x0)) {
                Console.WriteLine("not found");

                return (x, dx);
            }

            (yn2, y0, yp2) = (f(x0 - h), f(x0), f(x0 + h));

            h /= 2;
            while (h / (Math.Abs(x0) + double.Epsilon) >= 1e-15) {
                (yn1, yp1) = (f(x0 - h), f(x0 + h));

                Console.WriteLine($"{yn2}, {yn1}, {y0}, {yp1}, {yp2}");

                if (SequenceUtil.IsMonotone(yn2, yn1, y0, yp1, yp2)) {
                    Console.WriteLine("monotone");

                    if (yn1 * yp1 <= 0) {
                        Console.WriteLine("secant");

                        dx = (yp1 == yn1) ? 0 : 2 * h / (yp1 - yn1) * y0;
                        dx = Math.Max(-h, Math.Min(h, dx));

                        x0 -= dx;

                        if (Math.Abs(dx / (Math.Abs(x0) + double.Epsilon)) <= 1e-15) {
                            break;
                        }

                        y0 = f(x0);
                        if (y0 == 0) {
                            dx = 0;
                            break;
                        }

                        (yn2, yp2) = (f(x0 - h), f(x0 + h));
                        h /= 2;

                        continue;
                    }
                    if (yn2 * yn1 <= 0) {
                        Console.WriteLine("sft -");

                        x0 -= h;
                        (yn2, y0, yp2) = (f(x0 - h * 2), yn1, yp1);

                        continue;
                    }
                    if (yp2 * yp1 <= 0) {
                        Console.WriteLine("sft +");

                        x0 += h;
                        (yn2, y0, yp2) = (yn1, yp1, f(x0 + h * 2));

                        continue;
                    }
                }

                Console.WriteLine("shrink");

                (yn2, yp2) = (yn1, yp1);
                h /= 2;
            }

            return (x0, dx);
        }
    }
}
