using MultiPrecision;
using MultiPrecisionAlgebra;
using static MultiPrecision.Pow2;

namespace MathieuLimitSeries {
    class Program {
        static void Main() {
            const int p = 7, b = 36;
            using StreamWriter sw = new($"../../../../sandbox/eigen_limit_coef{p}_scaleb{b}.csv");

            sw.WriteLine("func,n,c");

            for (int n = 0; n <= 16; n++) {
                List<(MultiPrecision<N8> u, MultiPrecision<N8> a, MultiPrecision<N8> a_delta)> expecteds = ReadAExpected<N8>(n);

                MultiPrecision<N8> c = SearchA(n, expecteds, p, b);

                sw.WriteLine($"A,{n},{(long)c}");
                Console.WriteLine($"A,{n},{(long)c}");
            }

            for (int n = 1; n <= 16; n++) {
                List<(MultiPrecision<N8> u, MultiPrecision<N8> b, MultiPrecision<N8> b_delta)> expecteds = ReadBExpected<N8>(n);
            
                MultiPrecision<N8> c = SearchB(n, expecteds, p, b);

                sw.WriteLine($"B,{n},{(long)c}");
                Console.WriteLine($"B,{n},{(long)c}");
            }

            Console.WriteLine("END");
            Console.Read();
        }

        private static MultiPrecision<N8> SearchA(int n, List<(MultiPrecision<N8> u, MultiPrecision<N8> a, MultiPrecision<N8> a_delta)> expecteds, int p, int b) {
            Vector<N8> invu = expecteds.Select(item => item.u).ToArray();
            Vector<N8> q = expecteds.Select(item => MultiPrecision<N8>.Sqrt(1 / item.u) * Math.Max(1, n * n)).ToArray();

            Vector<N8> delta_expected = expecteds.Select(item => item.a_delta).ToArray();
            Vector<N8> delta_term6 = Vector<N8>.Func(q, v => LimitDeltaTerm6EigenA(n, v));

            Vector<N8> delta_rem = delta_expected - delta_term6;

            Vector<N8> delta_term6_unit = Vector<N8>.Func(q, v => 1 / (Math.ScaleB(1, b) * MultiPrecision<N8>.Pow(v, p * 0.5)));

            return SearchCoef(invu, delta_rem, delta_term6_unit);
        }

        private static MultiPrecision<N8> SearchB(int n, List<(MultiPrecision<N8> u, MultiPrecision<N8> b, MultiPrecision<N8> b_delta)> expecteds, int p, int b) {
            Vector<N8> invu = expecteds.Select(item => item.u).ToArray();
            Vector<N8> q = expecteds.Select(item => MultiPrecision<N8>.Sqrt(1 / item.u) * Math.Max(1, n * n)).ToArray();

            Vector<N8> delta_expected = expecteds.Select(item => item.b_delta).ToArray();
            Vector<N8> delta_term6 = Vector<N8>.Func(q, v => LimitDeltaTerm6EigenB(n, v));

            Vector<N8> delta_rem = delta_expected - delta_term6;

            Vector<N8> delta_term6_unit = Vector<N8>.Func(q, v => 1 / (Math.ScaleB(1, b) * MultiPrecision<N8>.Pow(v, p * 0.5)));

            return SearchCoef(invu, delta_rem, delta_term6_unit);
        }

        private static MultiPrecision<N8> SearchCoef(Vector<N8> invu, Vector<N8> delta_rem, Vector<N8> delta_term_unit) {
            Vector<N8> delta_scale = delta_rem / delta_term_unit;
            long exponent = delta_scale.Select(scale => scale.val.Exponent).Max();

            MultiPrecision<N8> c = 0, dc = MultiPrecision<N8>.Ldexp(1, exponent);
            while (dc >= 1) {
                Vector<N8> delta_term = checked(c + dc) * delta_term_unit;

                Vector<N8> delta_diff = delta_rem - delta_term;

                bool is_ok = true;

                for (int i = 0; i < delta_rem.Dim; i++) {
                    if (delta_diff[i] < 0 || (i > 0 && (delta_diff[i - 1] > delta_diff[i]))) {
                        is_ok = false;
                        break;
                    }
                }

                if(is_ok){
                    Vector<N8> delta_power = (delta_diff[1..] * invu[..^1]) / (delta_diff[..^1] * invu[1..]);

                    for (int i = 1; i < delta_power.Dim; i++) {
                        if (delta_power[i - 1] > delta_power[i] * 1.15) {
                            is_ok = false;
                            break;
                        }
                    }
                }

                if(!is_ok){
                    dc /= 2;
                    continue;
                }

                c += dc;
            }

            return c;
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

        
        static MultiPrecision<N> LimitDeltaTerm5EigenA<N>(int n, MultiPrecision<N> q) where N : struct, IConstant {
            int s = 2 * n + 1;

            MultiPrecision<N> h = MultiPrecision<N>.Sqrt(q);

            MultiPrecision<N> delta = DeltaTerm5(s, h);

            return delta;
        }

        static MultiPrecision<N> LimitDeltaTerm5EigenB<N>(int n, MultiPrecision<N> q) where N : struct, IConstant {
            int s = 2 * n - 1;

            MultiPrecision<N> h = MultiPrecision<N>.Sqrt(q);

            MultiPrecision<N> delta = DeltaTerm5(s, h);

            return delta;
        }

        static MultiPrecision<N> LimitDeltaTerm6EigenA<N>(int n, MultiPrecision<N> q) where N : struct, IConstant {
            int s = 2 * n + 1;

            MultiPrecision<N> h = MultiPrecision<N>.Sqrt(q);

            MultiPrecision<N>[] c6 = new MultiPrecision<N>[] {
                3888840           ,
                306376764         ,
                5813875953        ,
                53467223196       ,
                2964972550241     ,
                11327184240489    ,
                40588474402185    ,
                125453881876663   ,
                333685697794930   ,
                8400175517958246  ,
                11140113047006072 ,
                26316021447278539 ,
                61050830906379734 ,
                108046749099092452,
                188297945113537561,
                313180840876516265,
                520128832895466565,
            };

            MultiPrecision<N> delta = DeltaTerm6(s, h, MultiPrecision<N>.Ldexp(c6[n], -30));

            return delta;
        }

        static MultiPrecision<N> LimitDeltaTerm6EigenB<N>(int n, MultiPrecision<N> q) where N : struct, IConstant {
            int s = 2 * n - 1;

            MultiPrecision<N> h = MultiPrecision<N>.Sqrt(q);

            MultiPrecision<N>[] c6 = new MultiPrecision<N>[] {
                -1,
                2008817           ,
                306245940         ,
                5812787791        ,
                53460959011       ,
                3599731002398     ,
                13226576470455    ,
                46292840099754    ,
                140684377732158   ,
                369241964205361   ,
                10259299831890252 ,
                13261325522499640 ,
                30914292926139497 ,
                70936803722838360 ,
                124102996070214688,
                214191639546155697,
                353223126681787163,
            };

            MultiPrecision<N> delta = DeltaTerm6(s, h, MultiPrecision<N>.Ldexp(c6[n], -30));

            return delta;
        }

        private static MultiPrecision<N> DeltaTerm5<N>(int s, MultiPrecision<N> h) where N : struct, IConstant {
            MultiPrecision<N> delta = 0;
            delta += (s * (3 + s * s)) / (128 * h);
            delta += (9 + s * s * (34 + s * s * 5)) / (4096 * h * h);
            delta += (s * (405 + s * s * (410 + s * s * 33))) / (131072 * h * h * h);
            delta += (486 + s * s * (2943 + s * s * (1260 + s * s * 63))) / (1048576 * h * h * h * h);
            delta += (s * (41607 + s * s * (69001 + s * s * (15617 + s * s * 527)))) / (33554432 * h * h * h * h * h);

            return delta;
        }

        private static MultiPrecision<N> DeltaTerm6<N>(int s, MultiPrecision<N> h, MultiPrecision<N> c6) where N : struct, IConstant {
            MultiPrecision<N> delta = DeltaTerm5(s, h);
            delta += c6 / (h * h * h * h * h * h);

            return delta;
        }
    }
}