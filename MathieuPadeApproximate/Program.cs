using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionCurveFitting;

namespace MathieuPadeApproximate {
    class Program {
        static void Main() {
            List<(MultiPrecision<Pow2.N64> umin, MultiPrecision<Pow2.N64> umax, MultiPrecision<Pow2.N64> limit_range)> ranges = [
                (0, 1, 1 / 64d), (1, 2, 1 / 64d), (2, 4, 1 / 64d), (4, 8, 1 / 64d),
                (8, 16, 1 / 64d), (16, 32, 1 / 64d), (32, 64, 1 / 64d), (64, 128, 1 / 64d),
                (128, 256, 1 / 64d), (128, 256, 1 / 64d), (256, 512, 1 / 64d), (512, 1024, 1 / 64d), 
                (512, 1024, 1 / 64d), (1024, 2048, 1 / 64d), (2048, 4096, 1 / 64d)
            ];

            Dictionary<string, List<(MultiPrecision<Pow2.N64> u, MultiPrecision<Pow2.N64> v)>> expecteds = [];

            for (int n = 0; n <= 16; n++) {
                Console.WriteLine($"reading... n{n}");

                using StreamReader sr = new($"../../../../results_disused/eigen_precision145_n{n}.csv");

                List<(MultiPrecision<Pow2.N64> u, MultiPrecision<Pow2.N64> v)> expecteds_m = [], expecteds_d = [];

                sr.ReadLine();
                sr.ReadLine();
                sr.ReadLine();

                int scale = int.Max(1, n * n);

                while (!sr.EndOfStream) {
                    string? line = sr.ReadLine();
                    if (string.IsNullOrWhiteSpace(line)) {
                        break;
                    }

                    string[] line_split = line.Split(',');

                    MultiPrecision<Pow2.N64> u = line_split[0], m = line_split[3];

                    if (n > 0) {
                        m += 1;
                    }

                    expecteds_m.Add((u, m));
                }

                expecteds.Add($"m{n}", expecteds_m);
            }

            Dictionary<string, StreamWriter> sw_list = [];

            foreach (string func in expecteds.Keys) {
                StreamWriter sw = new($"../../../../results_disused/eigen_precision145_{func}_padecoef_2.csv");

                sw_list.Add(func, sw);
            }

            bool approximate(MultiPrecision<Pow2.N64> umin, MultiPrecision<Pow2.N64> umax) {
                Console.WriteLine($"[{umin}, {umax}]");

                Dictionary<string, List<(MultiPrecision<Pow2.N64> u, MultiPrecision<Pow2.N64> v)>> expecteds_range =
                    expecteds.Select(item => (item.Key, item.Value.Where(val => val.u >= umin && val.u <= umax).ToList())).ToDictionary();

                int samples = expecteds_range["m0"].Count;

                if (samples <= 64) {
                    return false;
                }

                foreach (bool forward in new bool[] { true, false }) {
                    if (umin == 0 && !forward) {
                        break;
                    }

                    Dictionary<string, (Vector<Pow2.N64> xs, Vector<Pow2.N64> ys)> value_table =
                        expecteds_range.Select(item => (item.Key, (
                            forward
                            ? new Vector<Pow2.N64>(item.Value.Select(item => item.u - umin).ToArray())
                            : new Vector<Pow2.N64>(item.Value.Select(item => umax - item.u).ToArray()),
                            new Vector<Pow2.N64>(item.Value.Select(item => item.v).ToArray())
                    ))).ToDictionary();

                    Dictionary<string, SumTable<Pow2.N64>> sum_tables =
                        value_table.Select(item => (item.Key, new SumTable<Pow2.N64>(item.Value.xs, item.Value.ys))).ToDictionary();

                    Dictionary<string, (PadeFitter<Pow2.N64> pade, Vector<Pow2.N64> param, MultiPrecision<Pow2.N64> max_err)> pade_tables = [];

                    Console.WriteLine($"expecteds computed : {samples}");

                    foreach ((string func, SumTable<Pow2.N64> sum_table) in sum_tables) {
                        Console.WriteLine(func);
                        Console.WriteLine($"{(forward ? "forward" : "backward")}");

                        Vector<Pow2.N64> xs = value_table[func].xs;
                        Vector<Pow2.N64> ys = value_table[func].ys;

                        bool has_root = umax < 9;

                        bool convergenced = false, over_precision = false;

                        for (int coefs = 5; coefs <= 64 && coefs < samples && !convergenced && !over_precision; coefs++) {
                            foreach ((int m, int n) in CurveFittingUtils.EnumeratePadeDegree(coefs, 4)) {
                                PadeFitter<Pow2.N64> pade = new(sum_table, m, n, intercept: umin == 0 ? ys[0] : null);

                                Vector<Pow2.N64> param = pade.Fit();
                                Vector<Pow2.N64> errs = pade.Error(param);

                                MultiPrecision<Pow2.N64> max_err = has_root
                                    ? CurveFittingUtils.MaxAbsoluteError(ys, pade.Regress(xs, param))
                                    : CurveFittingUtils.MaxRelativeError(ys, pade.Regress(xs, param));

                                Console.WriteLine($"m={m},n={n}");
                                Console.WriteLine($"{max_err:e20}");

                                if (max_err > "1e-22") {
                                    coefs += 4;
                                    break;
                                }

                                if (max_err < "1e-40") {
                                    over_precision = true;
                                    break;
                                }

                                if (has_root) {
                                    if (max_err < "1e-32" &&
                                        !CurveFittingUtils.HasLossDigitsPolynomialCoef(param[m..], 0, umax - umin)) {

                                        pade_tables.Add(func, (pade, param, max_err));
                                        convergenced = true;
                                        break;
                                    }
                                }
                                else {
                                    if (max_err < "1e-32" &&
                                        !CurveFittingUtils.HasLossDigitsPolynomialCoef(param[..m], 0, umax - umin) &&
                                        !CurveFittingUtils.HasLossDigitsPolynomialCoef(param[m..], 0, umax - umin)) {

                                        pade_tables.Add(func, (pade, param, max_err));
                                        convergenced = true;
                                        break;
                                    }
                                }
                            }
                        }

                        if (!convergenced) {
                            break;
                        }
                    }

                    if (pade_tables.Count != sum_tables.Count) {
                        continue;
                    }

                    foreach ((string func, (PadeFitter<Pow2.N64> pade, Vector<Pow2.N64> param, MultiPrecision<Pow2.N64> max_err)) in pade_tables) {
                        StreamWriter sw = sw_list[func];

                        sw.WriteLine($"u=[{umin},{umax}]");
                        sw.WriteLine($"samples={samples}");
                        sw.WriteLine($"{(forward ? "forward" : "backward")}");
                        sw.WriteLine($"m={pade.Numer},n={pade.Denom}");
                        sw.WriteLine("numer");
                        foreach (var (_, val) in param[..pade.Numer]) {
                            sw.WriteLine($"{val:e38}");
                        }
                        sw.WriteLine("denom");
                        foreach (var (_, val) in param[pade.Numer..]) {
                            sw.WriteLine($"{val:e38}");
                        }

                        sw.WriteLine("coef");
                        foreach ((var numer, var denom) in CurveFittingUtils.EnumeratePadeCoef(param, pade.Numer, pade.Denom)) {
                            sw.WriteLine($"({ToFP128(numer)}, {ToFP128(denom)}),");
                        }

                        sw.WriteLine("max err");
                        sw.WriteLine($"{max_err:e20}");
                        sw.Flush();
                    }

                    return true;
                }

                return false;
            }

            Segmenter<Pow2.N64> segmenter = new(ranges, approximate);
            segmenter.Execute();

            Console.WriteLine("END");
            Console.Read();
        }

        public static string ToFP128(MultiPrecision<Pow2.N64> x) {
            Sign sign = x.Sign;
            long exponent = x.Exponent;
            uint[] mantissa = x.Mantissa.Reverse().ToArray();

            string code = $"({(sign == Sign.Plus ? "+1" : "-1")}, {exponent}, 0x{mantissa[0]:X8}{mantissa[1]:X8}uL, 0x{mantissa[2]:X8}{mantissa[3]:X8}uL)";

            return code;
        }
    }
}