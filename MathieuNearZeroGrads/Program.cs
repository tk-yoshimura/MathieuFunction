using MultiPrecision;
using System.Collections.Generic;

namespace MathieuMP {
    class Program {
        static void Main() {
            int degrees = ForwardFiniteDifference<Pow2.N32, Pow2.N256>.SamplePoints - 1;
            
            for (int n = 0; n <= 12; n++) {
                Console.WriteLine($"Plotting {n}");

                List<(MultiPrecision<Pow2.N4> u, MultiPrecision<Pow2.N32> m, MultiPrecision<Pow2.N32> d)> values = ReadValues(n);

                List<MultiPrecision<Pow2.N32>[]> mgs = Grads(values.Select(item => (item.u, item.m)).ToList());
                List<MultiPrecision<Pow2.N32>[]> dgs = Grads(values.Select(item => (item.u, item.d)).ToList());

                using StreamWriter sw = new($"../../../../results/eigen_nearzero_grads_n{n}.csv");

                sw.WriteLine($"# zero shifted mathieu eigen value near zero grads n={n}");
                sw.WriteLine("# u:=q^2/max(1, n^4), m:=mean/max(1, n^2), d:=1/scaled_diff-1");
                sw.WriteLine("degree,h,m,d");

                for (int deg = 0; deg < degrees; deg++) {
                    for (int i = 0; i < mgs.Count; i++) {
                        sw.WriteLine($"{deg},{values[i].u.Convert<Pow2.N16>()},{mgs[i][deg]},{dgs[i][deg]}");
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        static List<(MultiPrecision<Pow2.N4> u, MultiPrecision<Pow2.N32> m, MultiPrecision<Pow2.N32> d)> ReadValues(int n) {
            List<(MultiPrecision<Pow2.N4> u, MultiPrecision<Pow2.N32> m, MultiPrecision<Pow2.N32> d)> res = new();
            
            using BinaryReader br = new(File.OpenRead($"../../../../results/eigen_nearzero_log2_precision1024bits_n{n}.bin"));

            for (int i = 0; i <= 128 + 8; i++) {
                MultiPrecision<Pow2.N4> u = br.ReadMultiPrecision<Pow2.N4>();
                MultiPrecision<Pow2.N32> m = br.ReadMultiPrecision<Pow2.N32>();
                MultiPrecision<Pow2.N32> d = br.ReadMultiPrecision<Pow2.N32>();

                res.Add((u, m, d));
            }

            return res;
        }

        static List<MultiPrecision<Pow2.N32>[]> Grads(List<(MultiPrecision<Pow2.N4> u, MultiPrecision<Pow2.N32> x)> values) {
            List<MultiPrecision<Pow2.N32>[]> res = new();

            int sample_points = ForwardFiniteDifference<Pow2.N32, Pow2.N256>.SamplePoints;

            for (int i = 0; i < values.Count - sample_points; i++) {
                MultiPrecision<Pow2.N32> h = values[i].u.Convert<Pow2.N32>();
                MultiPrecision<Pow2.N32>[] xs 
                    = new (MultiPrecision<Pow2.N4> u, MultiPrecision<Pow2.N32> x)[]{(0, 0)}
                    .Concat(values.Skip(i).Take(sample_points - 1)).Select(item => item.x).ToArray();

                MultiPrecision<Pow2.N32>[] gs = ForwardFiniteDifference<Pow2.N32, Pow2.N256>.Diff(xs, h);

                res.Add(gs);
            }

            return res;
        }
    }
}