using DoubleDouble;
using MultiPrecision;

namespace MathieuPadeDDouble {
    class Program {
        static void Main() {
            (ddouble umin, ddouble umax)[] ranges = new (ddouble, ddouble)[]{
                (0, Math.ScaleB(1, -20)),
                (Math.ScaleB(1, -20), Math.ScaleB(1, -16)),
                (Math.ScaleB(1, -16), Math.ScaleB(1, -12)),
            };

            using StreamWriter sw_com = new($"../../../../sandbox/eigen_limit_ddouble_padecoef.txt");

            for (int n = 0; n <= 16; n++) {
                using StreamWriter sw = new($"../../../../results/ddouble/eigen_limit_ddouble_results_a_n{n}.csv");
                sw.WriteLine("x,expected,approx,denom,error,relative_error");

                foreach ((ddouble umin, ddouble umax) in ranges) {
                    List<(ddouble invu, ddouble a, ddouble b)> expected = ReadExpected(n, umin, umax);
                    using StreamReader sr = new($"../../../../results/padecoef/eigen_limit_padecoef_precisionbits104_delta_range{umin}to{umax}_a_n{n}.csv");

                    int s = Math.Max(1, n * n);
                    ddouble[] limits = new ddouble[expected.Count];
                    for (int i = 0; i < limits.Length; i++) {
                        ddouble q = s * ddouble.Sqrt(1 / expected[i].invu);

                        ddouble limit = LimitEigenA(n, q);

                        limits[i] = limit;
                    }

                    ((ddouble[] ms, ddouble[] ns), ddouble u0) = ReadPadecoef(sr);

                    PlotResult(
                        sw,
                        expected.Select(item => ddouble.Sqrt(ddouble.Sqrt(item.invu))).ToArray(),
                        expected.Select(item => item.a).ToArray(),
                        limits,
                        ddouble.Sqrt(ddouble.Sqrt(u0)),
                        ms, ns
                    );
                }
            }

            for (int n = 0; n <= 16; n++) {
                sw_com.WriteLine($"public static ReadOnlyCollection<ReadOnlyCollection<(ddouble c, ddouble d)>> PadeA{n}Tables = Array.AsReadOnly(new ReadOnlyCollection<(ddouble c, ddouble d)>[] {{");

                foreach ((ddouble umin, ddouble umax) in ranges) {
                    using StreamReader sr = new($"../../../../results/padecoef/eigen_limit_padecoef_precisionbits104_delta_range{umin}to{umax}_a_n{n}.csv");

                    ((MultiPrecision<Pow2.N8>[] ms, MultiPrecision<Pow2.N8>[] ns), MultiPrecision<Pow2.N8> u0) = ReadPadecoef<Pow2.N8>(sr);

                    sw_com.WriteLine($"    new ReadOnlyCollection<(ddouble c, ddouble d)>(new (ddouble c, ddouble d)[] {{");

                    if (umin > 0) {
                        for (int i = 0; i < ns.Length; i++) {
                            sw_com.WriteLine($"        ({ToFP128(ms[i])}, {ToFP128(ns[i])}), ");
                        }
                    }
                    else {
                        for (int i = 0; i < ns.Length; i++) {
                            sw_com.WriteLine($"        ({ToFP128(ms[i + 1])}, {ToFP128(ns[i])}), ");
                        }
                    }

                    sw_com.WriteLine("    }),");
                }

                sw_com.WriteLine("});");
            }

            for (int n = 1; n <= 16; n++) {
                using StreamWriter sw = new($"../../../../results/ddouble/eigen_limit_ddouble_results_b_n{n}.csv");
                sw.WriteLine("x,expected,approx,denom,error,relative_error");

                foreach ((ddouble umin, ddouble umax) in ranges) {
                    List<(ddouble invu, ddouble a, ddouble b)> expected = ReadExpected(n, umin, umax);
                    using StreamReader sr = new($"../../../../results/padecoef/eigen_limit_padecoef_precisionbits104_delta_range{umin}to{umax}_b_n{n}.csv");

                    int s = Math.Max(1, n * n);
                    ddouble[] limits = new ddouble[expected.Count];
                    for (int i = 0; i < limits.Length; i++) {
                        ddouble q = s * ddouble.Sqrt(1 / expected[i].invu);

                        ddouble limit = LimitEigenB(n, q);

                        limits[i] = limit;
                    }

                    ((ddouble[] ms, ddouble[] ns), ddouble u0) = ReadPadecoef(sr);

                    PlotResult(
                        sw,
                        expected.Select(item => ddouble.Sqrt(ddouble.Sqrt(item.invu))).ToArray(),
                        expected.Select(item => item.b).ToArray(),
                        limits,
                        ddouble.Sqrt(ddouble.Sqrt(u0)),
                        ms, ns
                    );
                }
            }

            for (int n = 1; n <= 16; n++) {
                sw_com.WriteLine($"public static ReadOnlyCollection<ReadOnlyCollection<(ddouble c, ddouble d)>> PadeB{n}Tables = Array.AsReadOnly(new ReadOnlyCollection<(ddouble c, ddouble d)>[] {{");

                foreach ((ddouble umin, ddouble umax) in ranges) {
                    using StreamReader sr = new($"../../../../results/padecoef/eigen_limit_padecoef_precisionbits104_delta_range{umin}to{umax}_b_n{n}.csv");

                    ((MultiPrecision<Pow2.N8>[] ms, MultiPrecision<Pow2.N8>[] ns), MultiPrecision<Pow2.N8> u0) = ReadPadecoef<Pow2.N8>(sr);

                    sw_com.WriteLine($"    new ReadOnlyCollection<(ddouble c, ddouble d)>(new (ddouble c, ddouble d)[] {{");

                    if (umin > 0) {
                        for (int i = 0; i < ns.Length; i++) {
                            sw_com.WriteLine($"        ({ToFP128(ms[i])}, {ToFP128(ns[i])}), ");
                        }
                    }
                    else {
                        for (int i = 0; i < ns.Length; i++) {
                            sw_com.WriteLine($"        ({ToFP128(ms[i + 1])}, {ToFP128(ns[i])}), ");
                        }
                    }

                    sw_com.WriteLine("    }),");
                }

                sw_com.WriteLine("});");
            }

            Console.WriteLine("END");
            Console.Read();
        }

        static List<(ddouble invu, ddouble a, ddouble b)> ReadExpected(int n, ddouble umin, ddouble umax) {
            List<(ddouble invu, ddouble a, ddouble b)> res = [];

            using StreamReader sr = new($"../../../../results/eigen_limit_precision64_n{n}.csv");
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

                ddouble invu = line_split[0], a = line_split[1], b = line_split[4];

                if (invu > umax) {
                    break;
                }

                if (invu >= umin) {
                    res.Add((invu, a, b));
                }
            }

            return res;
        }

        static ((ddouble[] ms, ddouble[] ns) padecoef, ddouble u0) ReadPadecoef(StreamReader sr) {
            List<ddouble> ms = [], ns = [];

            string line = sr.ReadLine();
            if (!(line.StartsWith("u0 = "))) {
                throw new FormatException();
            }

            ddouble u0 = line[5..];

            if (!(sr.ReadLine().StartsWith("numers"))) {
                throw new FormatException();
            }

            while (!sr.EndOfStream) {
                line = sr.ReadLine();

                if (line.StartsWith("denoms")) {
                    break;
                }

                ms.Add(line);
            }

            while (!sr.EndOfStream) {
                line = sr.ReadLine();

                if (line == "1/u,expected,expected_delta,approx_delta,error") {
                    break;
                }

                ns.Add(line);
            }

            return ((ms.ToArray(), ns.ToArray()), u0);
        }

        static ((MultiPrecision<N>[] ms, MultiPrecision<N>[] ns) padecoef, MultiPrecision<N> u0) ReadPadecoef<N>(StreamReader sr) where N : struct, IConstant {
            List<MultiPrecision<N>> ms = [], ns = [];

            string line = sr.ReadLine();
            if (!(line.StartsWith("u0 = "))) {
                throw new FormatException();
            }

            MultiPrecision<N> u0 = line[5..];

            if (!(sr.ReadLine().StartsWith("numers"))) {
                throw new FormatException();
            }

            while (!sr.EndOfStream) {
                line = sr.ReadLine();

                if (line.StartsWith("denoms")) {
                    break;
                }

                ms.Add(line);
            }

            while (!sr.EndOfStream) {
                line = sr.ReadLine();

                if (line == "1/u,expected,expected_delta,approx_delta,error") {
                    break;
                }

                ns.Add(line);
            }

            return ((ms.ToArray(), ns.ToArray()), u0);
        }

        static void PlotResult(StreamWriter sw, ddouble[] xs, ddouble[] expecteds, ddouble[] limits, ddouble x0, ddouble[] ms, ddouble[] ns) {
            ddouble pade(ddouble x) {
                ddouble p = ms[^1], q = ns[^1];

                for (int i = ms.Length - 2; i >= 0; i--) {
                    p = x * p + ms[i];
                }

                for (int i = ns.Length - 2; i >= 0; i--) {
                    q = x * q + ns[i];
                }

                return p / q;
            }

            ddouble denom(ddouble x) {
                ddouble q = ns[^1];

                for (int i = ns.Length - 2; i >= 0; i--) {
                    q = x * q + ns[i];
                }

                return q;
            }

            for (int i = 0; i < expecteds.Length; i++) {
                ddouble x = xs[i], expected = expecteds[i], approx = limits[i] - pade(x - x0), d = denom(x - x0), error = expected - approx;
                ddouble relative_error = ddouble.Abs(error / expected);

                sw.WriteLine($"{x},{expected},{approx},{d},{error},{relative_error}");
            }
        }


        static ddouble LimitEigenA(int n, ddouble q) {
            int s = 2 * n + 1;
            ddouble h = ddouble.Sqrt(q);

            ddouble y = 2 * (s * h - q) - ddouble.Ldexp(6 * n * n + 2 * n + 1, -2) - DeltaTerm5(s, h);

            return y;
        }

        static ddouble LimitEigenB(int n, ddouble q) {
            int s = 2 * n - 1;
            ddouble h = ddouble.Sqrt(q);

            ddouble y = 2 * (s * h - q) - ddouble.Ldexp(6 * n * n - 2 * n + 1, -2) - DeltaTerm5(s, h);

            return y;
        }

        static ddouble DeltaTerm5(long s, ddouble h) {
            ddouble delta = 0;

            delta += (s * (3 + s * s)) / (128 * h);
            delta += (9 + s * s * (34 + s * s * 5)) / (4096 * h * h);
            delta += (s * (405 + s * s * (410 + s * s * 33))) / (131072 * h * h * h);
            delta += (486 + s * s * (2943 + s * s * (1260 + s * s * 63))) / (1048576 * h * h * h * h);
            delta += (checked(s * (41607 + s * s * (69001 + s * s * (15617 + s * s * 527))))) / (33554432 * h * h * h * h * h);

            return delta;
        }

        public static string ToFP128<N>(MultiPrecision<N> x) where N : struct, IConstant {
            Sign sign = x.Sign;
            long exponent = x.Exponent;
            uint[] mantissa = x.Mantissa.Reverse().ToArray();

            string code = $"({(sign == Sign.Plus ? "+1" : "-1")}, {exponent}, 0x{mantissa[0]:X8}{mantissa[1]:X8}uL, 0x{mantissa[2]:X8}{mantissa[3]:X8}uL)";

            return code;
        }
    }
}