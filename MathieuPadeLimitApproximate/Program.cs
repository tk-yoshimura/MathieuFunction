using MultiPrecision;
using MultiPrecisionAlgebra;
using static MultiPrecision.Pow2;

namespace MathieuPadeApproximate {
    class Program {
        static void Main() {
            {
                Dictionary<(MultiPrecision<N8> umin, MultiPrecision<N8> umax, MultiPrecision<N8> u0), int> ranges = new(){
                    { (0, Math.ScaleB(1, -20), 0) , 4 },
                    { (Math.ScaleB(1, -20), Math.ScaleB(1, -16), Math.ScaleB(81, -24)) , 4 },
                    { (Math.ScaleB(1, -16), Math.ScaleB(1, -12), Math.ScaleB(81, -20)) , 4 },
                };

                for (int n = 0; n <= 16; n++) {
                    foreach ((MultiPrecision<N8> umin, MultiPrecision<N8> umax, MultiPrecision<N8> u0) in ranges.Keys) {
                        SearchAndPlotA(n, ranges, umin, umax, u0);
                    }
                }
            }

            {
                Dictionary<(MultiPrecision<N8> umin, MultiPrecision<N8> umax, MultiPrecision<N8> u0), int> ranges = new(){
                    { (0, Math.ScaleB(1, -20), 0) , 4 },
                    { (Math.ScaleB(1, -20), Math.ScaleB(1, -16), Math.ScaleB(81, -24)) , 4 },
                    { (Math.ScaleB(1, -16), Math.ScaleB(1, -12), Math.ScaleB(81, -20)) , 4 },
                };

                for (int n = 1; n <= 16; n++) {
                    foreach ((MultiPrecision<N8> umin, MultiPrecision<N8> umax, MultiPrecision<N8> u0) in ranges.Keys) {
                        SearchAndPlotB(n, ranges, umin, umax, u0);
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        private static void SearchAndPlotA(int n, Dictionary<(MultiPrecision<N8> umin, MultiPrecision<N8> umax, MultiPrecision<N8> u0), int> ranges, MultiPrecision<N8> umin, MultiPrecision<N8> umax, MultiPrecision<N8> u0) {
            Console.WriteLine($"Plotting n={n} range=[{umin},{umax}]");

            List<(MultiPrecision<N32> u, MultiPrecision<N32> a, MultiPrecision<N32> a_delta)> expecteds = ReadAExpected(n, umin.Convert<N32>(), umax.Convert<N32>());

            Vector<N64> parameter, approx;

            for (ranges[(umin, umax, u0)] = 4; ranges[(umin, umax, u0)] <= 1024; ranges[(umin, umax, u0)]++) {
                int numer = ranges[(umin, umax, u0)] + (u0 == 0 ? 1 : 0), denom = ranges[(umin, umax, u0)];

                Console.WriteLine($"numer {numer} denom {denom}");

                (parameter, approx, bool success) = PadeApproximate<N64>(expecteds, u0, numer, denom);

                if (success) {
                    using StreamWriter sw = new($"../../../../results/padecoef/eigen_limit_padecoef_precisionbits104_delta_range{umin}to{umax}_a_n{n}.csv");
                    PlotResult(sw, expecteds, u0, numer, parameter, approx);
                    break;
                }
            }
        }

        private static void SearchAndPlotB(int n, Dictionary<(MultiPrecision<N8> umin, MultiPrecision<N8> umax, MultiPrecision<N8> u0), int> ranges, MultiPrecision<N8> umin, MultiPrecision<N8> umax, MultiPrecision<N8> u0) {
            Console.WriteLine($"Plotting n={n} range=[{umin},{umax}]");

            List<(MultiPrecision<N32> u, MultiPrecision<N32> b, MultiPrecision<N32> b_delta)> expecteds = ReadBExpected(n, umin.Convert<N32>(), umax.Convert<N32>());

            Vector<N64> parameter, approx;

            for (ranges[(umin, umax, u0)] = 4; ranges[(umin, umax, u0)] <= 1024; ranges[(umin, umax, u0)]++) {
                int numer = ranges[(umin, umax, u0)] + (u0 == 0 ? 1 : 0), denom = ranges[(umin, umax, u0)];

                Console.WriteLine($"numer {numer} denom {denom}");

                (parameter, approx, bool success) = PadeApproximate<N64>(expecteds, u0, numer, denom);

                if (success) {
                    using StreamWriter sw = new($"../../../../results/padecoef/eigen_limit_padecoef_precisionbits104_delta_range{umin}to{umax}_b_n{n}.csv");
                    PlotResult(sw, expecteds, u0, numer, parameter, approx);
                    break;
                }
            }
        }

        private static void PlotResult<N>(StreamWriter sw, List<(MultiPrecision<N32> u, MultiPrecision<N32> v, MultiPrecision<N32> v_delta)> expecteds, MultiPrecision<N8> u0, int numer, Vector<N> parameter, Vector<N> approx) where N : struct, IConstant {
            sw.WriteLine($"u0 = {u0}");
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

        static List<(MultiPrecision<N> u, MultiPrecision<N> a, MultiPrecision<N> a_delta)> ReadAExpected<N>(int n, MultiPrecision<N> umin, MultiPrecision<N> umax) where N : struct, IConstant {
            List<(MultiPrecision<N> u, MultiPrecision<N> a, MultiPrecision<N> a_delta)> res = new();

            using StreamReader sr = new($"../../../../results/eigen_limit_r4_precision64_n{n}.csv");
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

                if (u > umax) {
                    break;
                }

                if (u >= umin) {
                    res.Add((u, a, a_delta));
                }
            }

            return res;
        }

        static List<(MultiPrecision<N> u, MultiPrecision<N> b, MultiPrecision<N> b_delta)> ReadBExpected<N>(int n, MultiPrecision<N> umin, MultiPrecision<N> umax) where N : struct, IConstant {
            List<(MultiPrecision<N> u, MultiPrecision<N> b, MultiPrecision<N> a_delta)> res = new();

            using StreamReader sr = new($"../../../../results/eigen_limit_r4_precision64_n{n}.csv");
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

                if (u > umax) {
                    break;
                }

                if (u >= umin) {
                    res.Add((u, b, b_delta));
                }
            }

            return res;
        }

        static (Vector<N> parameter, Vector<N> approx, bool success) PadeApproximate<N>(List<(MultiPrecision<N32> u, MultiPrecision<N32> v, MultiPrecision<N32> v_delta)> expecteds, MultiPrecision<N8> u0, int numer, int denom) where N : struct, IConstant {
            long[] tol = expecteds.Select(expected => (MultiPrecision<N32>.Ldexp(expected.v.Exponent, -104) / expected.v_delta).Exponent).ToArray();

            bool[] needs_increase_weight(Vector<N> x, Vector<N> y, Vector<N> error) {
                Vector<N> relative_error = error / y;

                return relative_error.Select((e, idx) => e.val.Exponent > tol[idx]).ToArray();
            }

            Vector<N> xs = expecteds.Select((item) => MultiPrecision<N>.Sqrt(MultiPrecision<N>.Sqrt(item.u.Convert<N>()))).ToArray();
            Vector<N> ys = expecteds.Select((item) => item.v_delta.Convert<N>()).ToArray();
            Vector<N> weights = xs.Select(x => 1 / (x.val * x.val)).ToArray();

            xs -= MultiPrecision<N>.Sqrt(MultiPrecision<N>.Sqrt(u0.Convert<N>()));

            AdaptivePadeFitter<N> fitter = new(xs, ys, numer, denom, intercept: u0 == 0 ? 0 : null);

            (Vector<N> parameter, bool success) = fitter.ExecuteFitting(weights, needs_increase_weight);

            Vector<N> approx = fitter.FittingValue(xs, parameter);

            return (parameter, approx, success);
        }
    }
}