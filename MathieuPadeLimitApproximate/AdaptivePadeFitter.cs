using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionCurveFitting;

namespace MathieuPadeApproximate {
    public class AdaptivePadeFitter<N> : PadeFitter<N> where N : struct, IConstant {
        public AdaptivePadeFitter(Vector<N> xs, Vector<N> ys, int numer, int denom, MultiPrecision<N>? intercept = null)
            : base(xs, ys, numer, denom, intercept) { }

        public (Vector<N> parameters, bool success) ExecuteFitting(Vector<N> weights, Func<Vector<N>, Vector<N>, Vector<N>, bool[]> needs_increase_weight, MultiPrecision<N>? norm_cost = null, int iter = 64) {
            Vector<N> parameters = base.ExecuteFitting(weights, norm_cost);

            int min_increase_weight_count = Points;

            for (int j = 0; j < iter; j++) {
                Vector<N> error = Y - base.FittingValue(X, parameters);

                int increase_weight_count = 0;

                bool[] increase_weights = needs_increase_weight(X, Y, error);

                for (int i = 0; i < Points; i++) {
                    if (increase_weights[i]) {
                        weights[i] *= 10;
                        increase_weight_count++;
                    }
                }

                Console.WriteLine($"needs_increase_weight: {increase_weight_count}");

                if (j >= 2 && increase_weight_count > Points * 3 / 4) {
                    break;
                }
                if (increase_weight_count > Points / 16 && increase_weight_count > min_increase_weight_count * 3 / 2) {
                    break;
                }
                if (increase_weight_count <= 0) {
                    return (parameters, true);
                }

                min_increase_weight_count = Math.Min(increase_weight_count, min_increase_weight_count);

                parameters = base.ExecuteFitting(weights, norm_cost);
            }

            return (parameters, false);
        }
    }
}
