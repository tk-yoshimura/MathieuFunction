namespace MathieuMP {
    public static class SequenceUtil {
        /// <summary>
        /// Determine if the sequence is linear.
        /// </summary>
        public static bool IsLinear(double a, double b, double c, double d, double e) {
            double s = Math.Abs(a - e) / 8, t = Math.Abs(a - e) / 4;
            double p = Math.Abs(b - (a * 3 + e) / 4);
            double q = Math.Abs(c - (a + e) / 2);
            double r = Math.Abs(d - (a + e * 3) / 4);

            return (p <= s) && (q <= t) && (r <= s);
        }

        /// <summary>
        /// Determine if the sequence is monotone.
        /// </summary>
        public static bool IsMonotone(double a, double b, double c, double d, double e) {
            return (a <= b && b <= c && c <= d && d <= e) ||
                   (a >= b && b >= c && c >= d && d >= e);
        }

        /// <summary>
        /// Determine if the sequence is 1 / x like.
        /// </summary>
        public static bool IsReciprocal(double a, double b, double c, double d, double e) {
            return (a <= b && b > c && c <= d && d <= e) ||
                   (a <= b && b <= c && c > d && d <= e) ||
                   (a >= b && b < c && c >= d && d >= e) ||
                   (a >= b && b >= c && c < d && d >= e);
        }
    }
}
