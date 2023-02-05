using MultiPrecision;
using MultiPrecisionAlgebra;
using static MultiPrecision.Pow2;

namespace MathieuPadeApproximate {
    class Program {
        static void Main() {
            for (int n = 0; n <= 16; n++) {
                SearchAndPlotA(n);
            }

            for (int n = 1; n <= 16; n++) {
                SearchAndPlotB(n);
            }

            Console.WriteLine("END");
            Console.Read();
        }

        private static void SearchAndPlotA(int n) {
            List<(MultiPrecision<N32> u, MultiPrecision<N32> a, MultiPrecision<N32> a_delta)> expecteds = ReadAExpected<N32>(n);

            List<(MultiPrecision<N32> u, MultiPrecision<N32> a, MultiPrecision<N32> a_delta)> expecteds_transform = new();
            for (int i = 0; i < expecteds.Count; i++) {
                expecteds_transform.Add((1 / expecteds[i].u, expecteds[i].a, expecteds[i].a_delta));
            }

            expecteds_transform.Reverse();

            Vector<N32> parameter, approx;

            for (int k = 4; k <= 1024; k++) {
                int numer = k + 1, denom = k;

                Console.WriteLine($"numer {numer} denom {denom}");

                (parameter, approx, bool success) = PadeApproximate<N32>(expecteds_transform, numer, denom);

                if (success) {
                    using StreamWriter sw = new($"../../../../results/padecoef/eigen_limit_padecoef_precisionbits104_quad_delta_a_n{n}.csv");
                    PlotResult(sw, expecteds_transform, numer, parameter, approx);
                    break;
                }
            }
        }

        private static void SearchAndPlotB(int n) {
            List<(MultiPrecision<N32> u, MultiPrecision<N32> b, MultiPrecision<N32> b_delta)> expecteds = ReadBExpected<N32>(n);

            List<(MultiPrecision<N32> u, MultiPrecision<N32> b, MultiPrecision<N32> b_delta)> expecteds_transform = new();
            for (int i = 0; i < expecteds.Count; i++) {
                expecteds_transform.Add((1 / expecteds[i].u, expecteds[i].b, expecteds[i].b_delta));
            }

            expecteds_transform.Reverse();

            Vector<N32> parameter, approx;
            for (int k = 4; k <= 1024; k++) {
                int numer = k + 1, denom = k;

                Console.WriteLine($"numer {numer} denom {denom}");

                (parameter, approx, bool success) = PadeApproximate<N32>(expecteds_transform, numer, denom);

                if (success) {
                    using StreamWriter sw = new($"../../../../results/padecoef/eigen_limit_padecoef_precisionbits104_quad_delta_b_n{n}.csv");
                    PlotResult(sw, expecteds_transform, numer, parameter, approx);
                    break;
                }
            }
        }

        private static void PlotResult<N>(StreamWriter sw, List<(MultiPrecision<N32> u, MultiPrecision<N32> v, MultiPrecision<N32> v_delta)> expecteds, int numer, Vector<N> parameter, Vector<N> approx) where N : struct, IConstant {
            sw.WriteLine($"numers: {numer}");
            foreach ((_, MultiPrecision<N> v) in parameter[..numer]) {
                sw.WriteLine($"{v:e64}");
            }

            sw.WriteLine($"denoms: {(parameter.Dim - numer)}");
            foreach ((_, MultiPrecision<N> v) in parameter[numer..]) {
                sw.WriteLine($"{v:e64}");
            }

            sw.WriteLine("u,expected,expected_delta,approx_delta,error");
            for (int i = 0; i < expecteds.Count; i++) {
                sw.WriteLine($"{expecteds[i].u},{expecteds[i].v:e32},{expecteds[i].v_delta:e32},{approx[i]:e32},{(expecteds[i].v_delta - approx[i].Convert<N32>()):e10}");
            }
        }

        static List<(MultiPrecision<N> u, MultiPrecision<N> a, MultiPrecision<N> a_delta)> ReadAExpected<N>(int n) where N : struct, IConstant {
            List<(MultiPrecision<N> u, MultiPrecision<N> a, MultiPrecision<N> a_delta)> res = new();

            using StreamReader sr = new($"../../../../results/eigen_limit_r2_precision64_n{n}.csv");
            sr.ReadLine();
            sr.ReadLine();
            sr.ReadLine();

            while (!sr.EndOfStream) {
                string? line = sr.ReadLine();
                if (string.IsNullOrWhiteSpace(line)) {
                    break;
                }

                string[] line_split = line.Split(',');

                MultiPrecision<N> u = line_split[0], a = line_split[1], a_delta = line_split[3];

                res.Add((u, a, a_delta));
            }

            return res;
        }

        static List<(MultiPrecision<N> u, MultiPrecision<N> b, MultiPrecision<N> b_delta)> ReadBExpected<N>(int n) where N : struct, IConstant {
            List<(MultiPrecision<N> u, MultiPrecision<N> b, MultiPrecision<N> a_delta)> res = new();

            using StreamReader sr = new($"../../../../results/eigen_limit_r2_precision64_n{n}.csv");
            sr.ReadLine();
            sr.ReadLine();
            sr.ReadLine();

            while (!sr.EndOfStream) {
                string? line = sr.ReadLine();
                if (string.IsNullOrWhiteSpace(line)) {
                    break;
                }

                string[] line_split = line.Split(',');

                MultiPrecision<N> u = line_split[0], b = line_split[4], b_delta = line_split[6];

                res.Add((u, b, b_delta));
            }

            return res;
        }

        static (Vector<N> parameter, Vector<N> approx, bool success) PadeApproximate<N>(List<(MultiPrecision<N32> u, MultiPrecision<N32> v, MultiPrecision<N32> v_delta)> expecteds, int numer, int denom) where N : struct, IConstant {
            Vector<N> y = ((Vector<N32>)expecteds.Select(expected => expected.v).ToArray()).Convert<N>();

            bool[] needs_increase_weight(Vector<N> x, Vector<N> delta_quad) {
                Vector<N> delta = Vector<N>.Func(delta_quad, v => MultiPrecision<N>.Sqrt(MultiPrecision<N>.Sqrt(v)));

                Vector<N> error = Vector<N>.Func(y, delta, (expected, d) => MultiPrecision<N>.Abs(d / expected));

                return error.Select(e => e.val.Exponent > -104).ToArray();
            }

            Vector<N> xs = expecteds.Select((item) => item.u.Convert<N>()).ToArray();
            Vector<N> ys = expecteds.Select((item) => MultiPrecision<N>.Pow(item.v_delta.Convert<N>(), 4)).ToArray();
            Vector<N> weights = xs.Select(x => (MultiPrecision<N>)1).ToArray();

            AdaptivePadeFitter<N> fitter = new(xs, ys, numer, denom, intercept: 0);

            (Vector<N> parameter, bool success) = fitter.ExecuteFitting(weights, needs_increase_weight);
            Vector<N> approx = fitter.FittingValue(expecteds.Select((item) => item.u.Convert<N>()).ToArray(), parameter);

            return (parameter, approx, success);
        }
    }
}