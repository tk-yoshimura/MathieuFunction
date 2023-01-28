using MultiPrecision;

namespace MathieuMP {
    class Program {
        static void Main() {
            int degrees = ForwardFiniteDifference<N80, Pow2.N512>.SamplePoints - 1;

            for (int n = 0; n <= 32; n++) {
                Console.WriteLine($"Plotting {n}");

                List<(MultiPrecision<Pow2.N4> u, MultiPrecision<N80> m, MultiPrecision<N80> d)> values = ReadValues(n);

                List<MultiPrecision<N80>[]> mgs = Grads(values.Select(item => (item.u, item.m)).ToList());
                List<MultiPrecision<N80>[]> dgs = Grads(values.Select(item => (item.u, item.d)).ToList());

                using StreamWriter sw = new($"../../../../results/eigen_nearzero_grads_n{n}.csv");

                sw.WriteLine($"# zero shifted mathieu eigen value near zero grads n={n}");

                if (n >= 1) {
                    sw.WriteLine("# u:=q^2/n^4, m:=(a+b-2n^2)/2n^2 (1+u), d:=(q^n/(2^(n-1)(n-1)!)^2 2/(a-b)-1) (1+u)");
                }
                else {
                    sw.WriteLine("# u:=q^2, m:=a (1+u), d:=0");
                }

                sw.WriteLine("degree,m,m_precision,d,d_precision");

                sw.WriteLine("0,0,inf,0,inf");

                for (int deg = 0; deg < degrees; deg++) {
                    List<MultiPrecision<N80>> mg = new(), dg = new();
                    (MultiPrecision<N80> value, int precision) mg_sel = (MultiPrecision<N80>.NaN, 0), dg_sel = (MultiPrecision<N80>.NaN, 0);

                    for (int i = 0; i < mgs.Count; i++) {
                        mg.Add(mgs[i][deg]);
                        dg.Add(dgs[i][deg]);
                    }

                    int mg_precision = MultiPrecision<N80>.DecimalDigits, dg_precision = MultiPrecision<N80>.DecimalDigits;

                    for (; mg_precision >= 0 && mg_sel.precision <= 0; mg_precision--) {
                        string dec_prev = mg_precision >= 2 ? mg[0].ToString($"e{mg_precision}") : ((double)(mg[0])).ToString($"e{mg_precision}");
                        string dec = mg_precision >= 2 ? mg[1].ToString($"e{mg_precision}") : ((double)(mg[1])).ToString($"e{mg_precision}");

                        for (int i = 1; i < mg.Count - 1; i++) {
                            string dec_post = mg_precision >= 2 ? mg[i + 1].ToString($"e{mg_precision}") : ((double)(mg[i + 1])).ToString($"e{mg_precision}");

                            if (dec_prev == dec && dec == dec_post) {
                                mg_sel = (mg[i], mg_precision);
                                break;
                            }

                            (dec_prev, dec) = (dec, dec_post);
                        }
                    }

                    for (; dg_precision >= 0 && dg_sel.precision <= 0; dg_precision--) {
                        string dec_prev = dg_precision >= 2 ? dg[0].ToString($"e{dg_precision}") : ((double)(dg[0])).ToString($"e{dg_precision}");
                        string dec = dg_precision >= 2 ? dg[1].ToString($"e{dg_precision}") : ((double)(dg[1])).ToString($"e{dg_precision}");

                        for (int i = 1; i < dg.Count - 1; i++) {
                            string dec_post = dg_precision >= 2 ? dg[i + 1].ToString($"e{dg_precision}") : ((double)(dg[i + 1])).ToString($"e{dg_precision}");

                            if (dec_prev == dec && dec == dec_post) {
                                dg_sel = (dg[i], dg_precision);
                                break;
                            }

                            (dec_prev, dec) = (dec, dec_post);
                        }
                    }

                    if (n > 0) {
                        sw.WriteLine($"{deg + 1}," +
                            $"{mg_sel.value.ToString($"e{Math.Min(MultiPrecision<N80>.DecimalDigits, mg_sel.precision + 2)}")},{mg_sel.precision}," +
                            $"{dg_sel.value.ToString($"e{Math.Min(MultiPrecision<N80>.DecimalDigits, dg_sel.precision + 2)}")},{dg_sel.precision}"
                        );

                        Console.WriteLine($"{mg_sel.value:e40},{dg_sel.value:e40}");
                    }
                    else {
                        sw.WriteLine($"{deg + 1}," +
                            $"{mg_sel.value.ToString($"e{Math.Min(MultiPrecision<N80>.DecimalDigits, mg_sel.precision + 2)}")},{mg_sel.precision}," +
                            "0,inf"
                        );

                        Console.WriteLine($"{mg_sel.value:e40}");
                    }

                    sw.Flush();
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        static List<(MultiPrecision<Pow2.N4> u, MultiPrecision<N80> m, MultiPrecision<N80> d)> ReadValues(int n) {
            List<(MultiPrecision<Pow2.N4> u, MultiPrecision<N80> m, MultiPrecision<N80> d)> res = new();

            using BinaryReader br = new(File.OpenRead($"../../../../results/eigen_nearzero_log2_precision2560bits_n{n}.bin"));

            for (int i = 0; i <= 80 + 8; i++) {
                MultiPrecision<Pow2.N4> u = br.ReadMultiPrecision<Pow2.N4>();
                MultiPrecision<N80> m = br.ReadMultiPrecision<N80>();
                MultiPrecision<N80> d = br.ReadMultiPrecision<N80>();

                res.Add((u, m, d));
            }

            return res;
        }

        static List<MultiPrecision<N80>[]> Grads(List<(MultiPrecision<Pow2.N4> u, MultiPrecision<N80> x)> values) {
            List<MultiPrecision<N80>[]> res = new();

            int sample_points = ForwardFiniteDifference<N80, Pow2.N512>.SamplePoints;

            for (int i = 0; i < values.Count - sample_points; i++) {
                MultiPrecision<N80> h = values[i].u.Convert<N80>();
                MultiPrecision<N80>[] xs
                    = new (MultiPrecision<Pow2.N4> u, MultiPrecision<N80> x)[] { (0, 0) }
                    .Concat(values.Skip(i).Take(sample_points - 1)).Select(item => item.x).ToArray();

                MultiPrecision<N80>[] gs = ForwardFiniteDifference<N80, Pow2.N512>.Diff(xs, h);

                res.Add(gs);
            }

            return res;
        }
    }
}