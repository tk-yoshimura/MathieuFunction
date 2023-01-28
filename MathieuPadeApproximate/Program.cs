using MultiPrecision;
using MultiPrecisionAlgebra;
using static MultiPrecision.Pow2;

namespace MathieuPadeApproximate {
    class Program {
        static void Main() {
            int k = 16;
            for (int n = 0; n <= 32; n++) {
                if (n < 2) {
                    k = 16;
                }

                Console.WriteLine($"Plotting {n}");

                List<(MultiPrecision<N8> u, MultiPrecision<N8> m, MultiPrecision<N8> d)> expecteds = ReadExpected<N8>(n);

                Vector<N64> parameter, approx;

                for (; k <= 1024; k++) {
                    int numer = k + 1, denom = k;

                    Console.WriteLine($"numer {numer} denom {denom}");

                    (parameter, approx, bool success) = PadeApproximate<N64>(expecteds.Select(v => (v.u, v.m)).ToList(), numer, denom, has_nonzero_root: (n >= 2));

                    if (success) {
                        using StreamWriter sw = new($"../../../../sandbox/eigen_pade_results_m_n{n}.csv");
                        PlotResult(sw, expecteds, numer, parameter, approx);
                        break;
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        private static void PlotResult<N>(StreamWriter sw, List<(MultiPrecision<N8> u, MultiPrecision<N8> m, MultiPrecision<N8> d)> expecteds, int numer, Vector<N> parameter, Vector<N> approx) where N : struct, IConstant {
            sw.WriteLine("numers");
            foreach ((_, MultiPrecision<N> v) in parameter[..numer]) {
                sw.WriteLine($"{v:e64}");
            }

            sw.WriteLine("denoms");
            foreach ((_, MultiPrecision<N> v) in parameter[numer..]) {
                sw.WriteLine($"{v:e64}");
            }

            sw.WriteLine("x,expected,approx,error");
            for (int i = 0; i < expecteds.Count; i++) {
                sw.WriteLine($"{expecteds[i].u},{expecteds[i].m:e64},{approx[i]:e64},{(expecteds[i].m - approx[i].Convert<N8>()):e10}");
            }
        }

        static List<(MultiPrecision<N> u, MultiPrecision<N> m, MultiPrecision<N> d)> ReadExpected<N>(int n) where N : struct, IConstant {
            List<(MultiPrecision<N> u, MultiPrecision<N> m, MultiPrecision<N> d)> res = new();

            using StreamReader sr = new($"../../../../results/eigen_precision64_n{n}.csv");
            sr.ReadLine();
            sr.ReadLine();
            sr.ReadLine();

            while (!sr.EndOfStream) {
                string? line = sr.ReadLine();
                if (string.IsNullOrWhiteSpace(line)) {
                    break;
                }

                string[] line_split = line.Split(',');

                MultiPrecision<N> u = line_split[0], m = line_split[3], d = line_split[4];

                if (u > 16) {
                    break;
                }

                res.Add((u, m, d));
            }

            return res;
        }

        static (Vector<N> parameter, Vector<N> approx, bool success) PadeApproximate<N>(List<(MultiPrecision<N8> u, MultiPrecision<N8> v)> expecteds, int numer, int denom, bool has_nonzero_root) where N : struct, IConstant {
            Func<MultiPrecision<N>, MultiPrecision<N>, MultiPrecision<N>, bool> needs_increase_weight =
                has_nonzero_root
                    ? ((x, y, error) => {
                        if (x.Exponent >= 0 && y.Exponent <= -3) {
                            return error.Exponent > -103;
                        }
                        else {
                            return (error / MultiPrecision<N>.Abs(y)).Exponent > -100;
                        }
                    })
                    : ((x, y, error) => {
                        return (error / MultiPrecision<N>.Abs(y)).Exponent > -100;
                    });

            Vector<N> xs = expecteds.Skip(1).Select((item) => item.u.Convert<N>()).ToArray();
            Vector<N> ys = expecteds.Skip(1).Select((item) => item.v.Convert<N>()).ToArray();
            Vector<N> weights = xs.Select(x => 1 / (x.val * x.val * x.val * x.val)).ToArray();

            AdaptivePadeFitter<N> fitter = new(xs, ys, numer, denom, intercept: 0);

            (Vector<N> parameter, bool success) = fitter.ExecuteFitting(weights, needs_increase_weight);
            Vector<N> approx = fitter.FittingValue(expecteds.Select((item) => item.u.Convert<N>()).ToArray(), parameter);

            return (parameter, approx, success);
        }
    }
}