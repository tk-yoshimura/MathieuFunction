namespace MathieuMP {
    public static class SequenceUtil {
        /// <summary>
        /// Determine if the sequence is monotone.
        /// </summary>
        public static bool IsMonotone(double a, double b, double c, double d, double e) {
            if (double.IsFinite(a) && double.IsFinite(e)) {
                return (a <= b && b <= c && c <= d && d <= e) ||
                       (a >= b && b >= c && c >= d && d >= e);
            }
            if (double.IsFinite(a)) { 
                return (a <= b && b <= c && c <= d) ||
                       (a >= b && b >= c && c >= d);
            }
            if (double.IsFinite(e)) {
                return (b <= c && c <= d && d <= e) ||
                       (b >= c && c >= d && d >= e);
            }

            return (b <= c && c <= d) || (b >= c && c >= d);
        }
    }
}
