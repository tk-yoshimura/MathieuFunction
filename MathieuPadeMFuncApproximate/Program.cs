﻿using MultiPrecision;
using MultiPrecisionAlgebra;
using static MultiPrecision.Pow2;

namespace MathieuPadeApproximate {
    class Program {
        static void Main() {
            Dictionary<(MultiPrecision<N8> umin,  MultiPrecision<N8> umax), int> ranges = new(){
                { (0, 1d / 256), 4},
                { (1d / 256, 1d / 16) , 4 },
                { (1d / 16, 0.25) , 4 },
                { (0.25, 1) , 4 },
                { (1, 4) , 4 },
                { (4, 16) , 4 },
                { (16, 64) , 4 },
                { (64, 1024) , 4 },
                //{ (0, 1d / 64), 4},
                //{ (1d / 64, 1d / 16), 4},
                //{ (0, 1d / 16) , 4 },
            };

            for (int n = 0; n <= 64; n++) {
                foreach ((MultiPrecision<N8> umin, MultiPrecision<N8> umax) in ranges.Keys) {
                    if (n < 2) {
                        ranges[(umin, umax)] = 4;
                    }

                    SearchAndPlot(n, ranges, umin, umax);
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        private static void SearchAndPlot(int n, Dictionary<(MultiPrecision<N8> umin, MultiPrecision<N8> umax), int> ranges,  MultiPrecision<N8> umin, MultiPrecision<N8> umax) {
            Console.WriteLine($"Plotting n={n} range=[{umin},{umax}]");

            List<(MultiPrecision<N8> u, MultiPrecision<N8> m)> expecteds = ReadExpected(n, umin, umax);

            Vector<N32> parameter, approx;
            bool success = false;

            for (; ranges[(umin, umax)] <= 1024 && !success; ranges[(umin, umax)]++) {
                int numer = ranges[(umin, umax)], denom = ranges[(umin, umax)];

                Console.WriteLine($"numer {numer} denom {denom}");

                (parameter, approx, success) = PadeApproximate<N32>(expecteds.Select(v => (v.u, v.m)).ToList(), numer, denom, has_nonzero_root: (n >= 1));

                if (success) {
                    using StreamWriter sw = new($"../../../../results/eigen_padecoef_precisionbits104_range{umin}to{umax}_m_n{n}.csv");
                    PlotResult(sw, expecteds, numer, parameter, approx);
                    break;
                }
            }
        }

        private static void PlotResult<N>(StreamWriter sw, List<(MultiPrecision<N8> u, MultiPrecision<N8> m)> expecteds, int numer, Vector<N> parameter, Vector<N> approx) where N : struct, IConstant {
            MultiPrecision<N> u0 = expecteds[0].u.Convert<N>();

            sw.WriteLine($"u0 = {u0}");
            sw.WriteLine($"numers: {numer}");
            foreach ((_, MultiPrecision<N> v) in parameter[..numer]) {
                sw.WriteLine($"{v:e64}");
            }

            sw.WriteLine($"denoms: {(parameter.Dim - numer)}");
            foreach ((_, MultiPrecision<N> v) in parameter[numer..]) {
                sw.WriteLine($"{v:e64}");
            }

            sw.WriteLine("u,expected,approx,error");
            for (int i = 0; i < expecteds.Count; i++) {
                sw.WriteLine($"{expecteds[i].u},{expecteds[i].m:e32},{approx[i]:e32},{(expecteds[i].m - approx[i].Convert<N8>()):e10}");
            }
        }

        static List<(MultiPrecision<N> u, MultiPrecision<N> m)> ReadExpected<N>(int n,  MultiPrecision<N> umin,  MultiPrecision<N> umax) where N : struct, IConstant {
            List<(MultiPrecision<N> u, MultiPrecision<N> m)> res = new();

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

                MultiPrecision<N> u = line_split[0], m = line_split[3];

                if (n >= 1) {
                    m += 1;
                }

                if (u > umax) {
                    break;
                }

                if (u >= umin) {
                    res.Add((u, m));
                }
            }

            return res;
        }

        static (Vector<N> parameter, Vector<N> approx, bool success) PadeApproximate<N>(List<(MultiPrecision<N8> u, MultiPrecision<N8> v)> expecteds, int numer, int denom, bool has_nonzero_root) where N : struct, IConstant {
            MultiPrecision<N> x0 = expecteds[0].u.Convert<N>();

            Func<MultiPrecision<N>, MultiPrecision<N>, MultiPrecision<N>, bool> needs_increase_weight =
                has_nonzero_root
                    ? ((x, y, error) => {
                        if ((x + x0).Exponent >= 0 && y.Exponent <= -3) {
                            return error.Exponent > -107;
                        }
                        else {
                            return (error / MultiPrecision<N>.Abs(y)).Exponent > -104;
                        }
                    })
                    : ((x, y, error) => {
                        return (error / MultiPrecision<N>.Abs(y)).Exponent > -104;
                    });

            Vector<N> xs = expecteds.Select((item) => item.u.Convert<N>()).ToArray();
            Vector<N> ys = expecteds.Select((item) => item.v.Convert<N>()).ToArray();
            Vector<N> weights = xs.Select(x => 1 / (x.val * x.val * x.val * x.val)).ToArray();

            if (x0 == 0) {
                (xs, ys, weights) = (xs[1..], ys[1..], weights[1..]);
            }
            else { 
                xs -= x0;
            }

            AdaptivePadeFitter<N> fitter = new(xs, ys, numer, denom, intercept: (x0 <= 0) ? expecteds[0].v.Convert<N>() : null);

            (Vector<N> parameter, bool success) = fitter.ExecuteFitting(weights, needs_increase_weight);
            Vector<N> approx = fitter.FittingValue(expecteds.Select((item) => item.u.Convert<N>() - x0).ToArray(), parameter);

            return (parameter, approx, success);
        }
    }
}