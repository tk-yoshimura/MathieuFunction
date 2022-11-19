namespace MathieuMP {
    public static class MathieuEigenFP64 {
        public static double FractionA(int n, double q, double a_init, int terms = 32) {
            if (!(q >= 0)) {
                throw new ArgumentException($"{nameof(n)}, {nameof(q)}");
            }

            double a = a_init;
            double h = q * q;
            double u = ((n & 1) == 0)
                ? (n == 0) ? 0d : 2d / (a + checked(n * n))
                : (n == 1) ? 1d / q : 1d / (a + checked(n * n - 1) - q);
            double v = 0d;

            for (int k = 2 + (n & 1); k < n; k += 2) {
                u = 1d / (a + checked(n * n - k * k) - u * h);
            }
            for (int k = checked(terms * 2 + n); k > n; k -= 2) {
                v = 1d / (a + checked(n * n - k * k) - v * h);
            }

            if (n == 0) {
                v *= 2d;
            }

            double y = (u + v) * h - a;

            return y;
        }

        public static double FractionB(int n, double q, double a_init, int terms = 32) {
            if (n < 1 || !(q >= 0)) {
                throw new ArgumentException($"{nameof(n)}, {nameof(q)}");
            }

            double a = a_init;
            double h = q * q;
            double u = ((n & 1) == 0)
                ? 0d
                : (n == 1) ? -1d / q : 1d / (a + checked(n * n - 1) + q);
            double v = 0d;

            for (int k = 2 + (n & 1); k < n; k += 2) {
                u = 1d / (a + checked(n * n - k * k) - u * h);
            }
            for (int k = checked(terms * 2 + n); k > n; k -= 2) {
                v = 1d / (a + checked(n * n - k * k) - v * h);
            }

            double y = (u + v) * h - a;

            return y;
        }
    }
}
