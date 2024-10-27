using DoubleDouble;
using MultiPrecision;

namespace MathieuPadeDDouble {
    class Program {
        static void Main() {
            (ddouble umin, ddouble umax)[] ranges = new (ddouble, ddouble)[]{
                (0, 1d / 4) ,
            };

            using StreamWriter sw_com = new($"../../../../sandbox/eigen_nearzero_ddouble_padecoef.txt");

            for (int n = 1; n <= 16; n++) {
                using StreamWriter sw = new($"../../../../results/ddouble/eigen_nearzero_ddouble_results_az_n{n}.csv");
                sw.WriteLine("x,expected,approx,denom,error,relative_error");

                foreach ((ddouble umin, ddouble umax) in ranges) {
                    List<(ddouble u, ddouble m, ddouble d)> expected = ReadExpected(n, umin, umax);
                    using StreamReader sr = new($"../../../../results/padecoef/eigen_nearzero_padecoef_precisionbits104_range{umin}to{umax}_az_n{n}.csv");

                    ((ddouble[] ms, ddouble[] ns), ddouble u0) = ReadPadecoef(sr);

                    PlotResult(
                        sw,
                        expected.Select(item => item.u).ToArray(),
                        expected.Select(item => item.m).ToArray(),
                        u0,
                        ms, ns
                    );
                }
            }

            for (int n = 1; n <= 16; n++) {
                sw_com.WriteLine($"public static ReadOnlyCollection<ReadOnlyCollection<(ddouble c, ddouble d)>> PadeAz{n}Tables = Array.AsReadOnly(new ReadOnlyCollection<(ddouble c, ddouble d)>[] {{");

                foreach ((ddouble umin, ddouble umax) in ranges) {
                    using StreamReader sr = new($"../../../../results/padecoef/eigen_nearzero_padecoef_precisionbits104_range{umin}to{umax}_az_n{n}.csv");

                    ((MultiPrecision<Pow2.N8>[] ms, MultiPrecision<Pow2.N8>[] ns), MultiPrecision<Pow2.N8> u0) = ReadPadecoef<Pow2.N8>(sr);

                    sw_com.WriteLine($"    new ReadOnlyCollection<(ddouble c, ddouble d)>(new (ddouble c, ddouble d)[] {{");
                    for (int i = 0; i < ms.Length; i++) {
                        sw_com.WriteLine($"        ({ToFP128(ms[i])}, {ToFP128(ns[i])}), ");
                    }
                    sw_com.WriteLine("    }),");
                }

                sw_com.WriteLine("});");
            }

            for (int n = 1; n <= 16; n++) {
                using StreamWriter sw = new($"../../../../results/ddouble/eigen_nearzero_ddouble_results_bz_n{n}.csv");
                sw.WriteLine("x,expected,approx,denom,error,relative_error");

                foreach ((ddouble umin, ddouble umax) in ranges) {
                    List<(ddouble u, ddouble m, ddouble d)> expected = ReadExpected(n, umin, umax);
                    using StreamReader sr = new($"../../../../results/padecoef/eigen_nearzero_padecoef_precisionbits104_range{umin}to{umax}_bz_n{n}.csv");

                    ((ddouble[] ms, ddouble[] ns), ddouble u0) = ReadPadecoef(sr);

                    PlotResult(
                        sw,
                        expected.Select(item => item.u).ToArray(),
                        expected.Select(item => item.d).ToArray(),
                        u0,
                        ms, ns
                    );
                }
            }

            for (int n = 1; n <= 16; n++) {
                sw_com.WriteLine($"public static ReadOnlyCollection<ReadOnlyCollection<(ddouble c, ddouble d)>> PadeBz{n}Tables = Array.AsReadOnly(new ReadOnlyCollection<(ddouble c, ddouble d)>[] {{");

                foreach ((ddouble umin, ddouble umax) in ranges) {
                    using StreamReader sr = new($"../../../../results/padecoef/eigen_nearzero_padecoef_precisionbits104_range{umin}to{umax}_bz_n{n}.csv");

                    ((MultiPrecision<Pow2.N8>[] ms, MultiPrecision<Pow2.N8>[] ns), MultiPrecision<Pow2.N8> u0) = ReadPadecoef<Pow2.N8>(sr);

                    sw_com.WriteLine($"    new ReadOnlyCollection<(ddouble c, ddouble d)>(new (ddouble c, ddouble d)[] {{");
                    for (int i = 0; i < ms.Length; i++) {
                        sw_com.WriteLine($"        ({ToFP128(ms[i])}, {ToFP128(ns[i])}), ");
                    }
                    sw_com.WriteLine("    }),");
                }

                sw_com.WriteLine("});");
            }

            Console.WriteLine("END");
            Console.Read();
        }

        static List<(ddouble u, ddouble az, ddouble bz)> ReadExpected(int n, ddouble umin, ddouble umax) {
            List<(ddouble u, ddouble az, ddouble bz)> res = [];

            using StreamReader sr = new($"../../../../results/eigen_nearzero_precision64_n{n}.csv");
            sr.ReadLine();
            sr.ReadLine();
            sr.ReadLine();

            while (!sr.EndOfStream) {
                string? line = sr.ReadLine();
                if (string.IsNullOrWhiteSpace(line)) {
                    break;
                }

                string[] line_split = line.Split(',');

                ddouble u = line_split[0], az = line_split[5], bz = line_split[6];

                if (u > umax) {
                    break;
                }

                if (u >= umin) {
                    res.Add((u, az, bz));
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

                if (line == "u,expected,approx,error") {
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

                if (line == "u,expected,approx,error") {
                    break;
                }

                ns.Add(line);
            }

            return ((ms.ToArray(), ns.ToArray()), u0);
        }

        static void PlotResult(StreamWriter sw, ddouble[] xs, ddouble[] expecteds, ddouble x0, ddouble[] ms, ddouble[] ns) {
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
                ddouble x = xs[i], expected = expecteds[i], approx = pade(x - x0), d = denom(x - x0), error = expected - approx;
                ddouble relative_error = ((x > 0 && x < 1) || ddouble.Abs(expected) > 0.125) ? ddouble.Abs(error / expected) : ddouble.Abs(error);

                sw.WriteLine($"{x},{expected},{approx},{d},{error},{relative_error}");
            }
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