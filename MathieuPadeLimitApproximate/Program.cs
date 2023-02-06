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

            Vector<N32> parameter, approx;

            for (int k = 4; k <= 1024; k++) {
                int numer = k + 1, denom = k;

                Console.WriteLine($"numer {numer} denom {denom}");

                (parameter, approx, bool success) = PadeApproximate<N32>(expecteds, numer, denom);

                if (success) {
                    using StreamWriter sw = new($"../../../../results/padecoef/eigen_limit_padecoef_precisionbits104_quad_delta_a_n{n}.csv");
                    PlotResult(sw, expecteds, numer, parameter, approx);
                    break;
                }
            }
        }

        private static void SearchAndPlotB(int n) {
            List<(MultiPrecision<N32> u, MultiPrecision<N32> b, MultiPrecision<N32> b_delta)> expecteds = ReadBExpected<N32>(n);

            Vector<N32> parameter, approx;
            for (int k = 4; k <= 1024; k++) {
                int numer = k + 1, denom = k;

                Console.WriteLine($"numer {numer} denom {denom}");

                (parameter, approx, bool success) = PadeApproximate<N32>(expecteds, numer, denom);

                if (success) {
                    using StreamWriter sw = new($"../../../../results/padecoef/eigen_limit_padecoef_precisionbits104_quad_delta_b_n{n}.csv");
                    PlotResult(sw, expecteds, numer, parameter, approx);
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

            sw.WriteLine("1/u,expected,expected_delta,approx_delta,error");
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
            MultiPrecision<N> s = MultiPrecision<N>.Pow(expecteds.Select((item) => item.v_delta).Max().Convert<N>(), 4);
            
            Vector<N> y = ((Vector<N32>)expecteds.Select(expected => expected.v).ToArray()).Convert<N>();

            bool[] needs_increase_weight(Vector<N> x, Vector<N> y, Vector<N> error) {
                Vector<N> relative_error = error / y;

                return relative_error.Select(e => e.val.Exponent > -100).ToArray();
            }

            Vector<N> xs = expecteds.Select((item) => item.u.Convert<N>() * 1024).ToArray();
            Vector<N> ys = expecteds.Select((item) => MultiPrecision<N>.Pow(item.v_delta.Convert<N>(), 4)).ToArray();
            Vector<N> weights = xs.Select(x => 1 / (x.val * x.val)).ToArray();

            AdaptivePadeFitter<N> fitter = new(xs, ys / s, numer, denom, intercept: 0);

            (Vector<N> parameter, bool success) = fitter.ExecuteFitting(weights, needs_increase_weight);
            parameter[..numer] *= s;

            Vector<N> approx = fitter.FittingValue(xs, parameter);
            approx = Vector<N>.Func(approx, v => MultiPrecision<N>.Sqrt(MultiPrecision<N>.Sqrt(v)));

            return (parameter, approx, success);
        }
    }
}