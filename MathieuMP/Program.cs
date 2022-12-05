using System.Xml.Serialization;

namespace MathieuMP {
    class Program {
        static void Main() {
             for (double q = 0; q <= 2048; q += 1) {
                double a = EigenFP64.InitialValue(EigenFunc.B, 34, q);
                double b = EigenFP64.Value(EigenFunc.B, 34, q, zero_shift: true).value;
                Console.WriteLine($"{q},{a},{b}");
             }

            //foreach (EigenFunc func in new[] { EigenFunc.A, EigenFunc.B }) {
            //    using StreamWriter sw = new($"../../../../results/eigen_{func}_approx.csv");
            //
            //    SortedDictionary<int, (decimal s, decimal t)> table = new();
            //
            //    for (int n = func == EigenFunc.A ? 0: 1; n <= 32; n++) {
            //        double min_sd = double.PositiveInfinity;
            //        decimal min_s = 0.1m, min_t = 0.5m, d = 1m;
            //        (decimal min, decimal max) range_s = (0.5m, 8.5m), range_t = (0.5m, 16.5m);
            //
            //        while (d >= 1 / 128m) {
            //            for (decimal s = range_s.min; s <= range_s.max; s += d) {
            //                for (decimal t = range_t.min; t <= range_t.max; t += d) {
            //
            //                    if (s > t) {
            //                        continue;
            //                    }
            //
            //                    double sd = 0;
            //                    List<double> xs = new(), ys = new();
            //                    for (double q = 0; q <= 8 * Math.Max(1, n * n); q += Math.Max(1, n * n) / 128d) {
            //                        double x = EigenFP64.InitialValueTest(func, n, q, (double)s, (double)t);
            //                        double y = EigenFP64.SearchFit(func, n, q, x).value;
            //
            //                        xs.Add(x);
            //                        ys.Add(y);
            //
            //                        sd += Math.Abs(x - y);
            //
            //                        if (!(sd < min_sd)) {
            //                            break;
            //                        }
            //                        //Console.WriteLine($"{q},{x}");
            //                    }
            //
            //                    //double sd = 0;
            //                    //for (int i = 1; i < xs.Count - 1; i++) {
            //                    //    double d = Math.Abs(xs[i - 1] + xs[i + 1] - 2 * xs[i]);
            //                    //    sd += d;
            //                    //}
            //
            //                    if (double.IsFinite(sd) && sd < min_sd) {
            //                        (min_sd, min_s, min_t) = (sd, s, t);
            //                    }
            //
            //                    Console.Write('.');
            //                }
            //            }
            //
            //            Console.WriteLine($"\ns={min_s}, t={min_t}: {min_sd}");
            //
            //            range_s = (Math.Max(0, min_s - d), min_s + d);
            //            range_t = (Math.Max(0, min_t - d), min_t + d);
            //            d /= 2;
            //        }
            //
            //        Console.WriteLine($"n={n}, s={min_s}, t={min_t}");
            //        sw.WriteLine($"n={n}, s={min_s}, t={min_t}");
            //
            //        for (double q = 0; q <= 8 * Math.Max(1, n * n); q += Math.Max(1, n * n) / 16d) {
            //            double x = EigenFP64.InitialValueTest(func, n, q, (double)min_s, (double)min_t);
            //            double y = EigenFP64.SearchFit(func, n, q, x).value;
            //
            //            sw.WriteLine($"{q},{x},{y}");
            //        }
            //
            //        sw.Flush();
            //
            //        table.Add(n, (min_s, min_t));
            //    }
            //
            //    sw.WriteLine("n,s,t");
            //    foreach (var item in table) {
            //        sw.WriteLine($"{item.Key},{item.Value.s},{item.Value.t}");
            //    }
            //}
            //
            Console.WriteLine("END");
            Console.Read();
        }
    }
}