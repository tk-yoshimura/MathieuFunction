﻿namespace MathieuMP {
    public static class SequenceUtil {
        /// <summary>
        /// Determine if the sequence is monotone.
        /// </summary>
        public static bool IsMonotone(double a, double b, double c, double d, double e) {
            return (a <= b && b <= c && c <= d && d <= e) ||
                   (a >= b && b >= c && c >= d && d >= e);
        }

        /// <summary>
        /// Determine if the sequence is 1/x like curve.
        /// </summary>
        public static bool IsReciprocalCurve(double v0, double v1, double v2, double v3, double v4, double v5, double v6) {

            // is monotone
            if (v0 <= v1 && !(v1 <= v2 && v2 <= v3 && v3 <= v4 && v4 <= v5 && v5 <= v6)) {
                return false;
            }
            if (v0 >= v1 && !(v1 >= v2 && v2 >= v3 && v3 >= v4 && v4 >= v5 && v5 >= v6)) {
                return false;
            }

            // is linear
            if ((Math.Abs((v0 + v2) / 2 - v1) <= Math.Abs((v0 - v2) / 256)) &&
                (Math.Abs((v1 + v3) / 2 - v2) <= Math.Abs((v1 - v3) / 256)) &&
                (Math.Abs((v2 + v4) / 2 - v3) <= Math.Abs((v2 - v4) / 256)) &&
                (Math.Abs((v3 + v5) / 2 - v4) <= Math.Abs((v3 - v5) / 256)) &&
                (Math.Abs((v4 + v6) / 2 - v5) <= Math.Abs((v4 - v6) / 256))) {

                return true;
            }

            // is convex
            if (v0 + v2 <= 2 * v1 && !(v1 + v3 <= 2 * v2 && v2 + v4 <= 2 * v3 && v3 + v5 <= 2 * v4 && v4 + v6 <= 2 * v5)) {
                return false;
            }
            if (v0 + v2 >= 2 * v1 && !(v1 + v3 >= 2 * v2 && v2 + v4 >= 2 * v3 && v3 + v5 >= 2 * v4 && v4 + v6 >= 2 * v5)) {
                return false;
            }

            double d01 = v0 - v1, d12 = v1 - v2, d23 = v2 - v3, d34 = v3 - v4, d45 = v4 - v5, d56 = v5 - v6;

            // is diff monotone
            if (d01 <= d12 && !(d12 <= d23 && d23 <= d34 && d34 <= d45 && d45 <= d56)) {
                return false;
            }
            if (d01 >= d12 && !(d12 >= d23 && d23 >= d34 && d34 >= d45 && d45 >= d56)) {
                return false;
            }

            return true;
        }
    }
}
