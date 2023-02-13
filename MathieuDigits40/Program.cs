using MultiPrecision;
using System.Numerics;
using static MultiPrecision.Pow2;

namespace MathieuMP {
    class Program {
        static void Main() {
            List<MultiPrecision<N8>> us = new();

            for (MultiPrecision<N8> u = MultiPrecision<N8>.Ldexp(1, -20); u <= MultiPrecision<N8>.Ldexp(1, 80); u *= 2) {
                us.Add(u);
            }

            for (int n = 0; n <= 16; n++) {
                Console.WriteLine($"Plotting {n}");

                using StreamWriter sw = new($"../../../../results/eigen_precision40_n{n}.csv");
                sw.WriteLine($"# mathieu eigen value precision_digits=40 n={n}");

                if (n >= 1) {
                    sw.WriteLine("# u:=q^2/n^4");
                }
                else {
                    sw.WriteLine("# u:=q^2");
                }

                sw.WriteLine("u^2,q,a,b");

                int s = Math.Max(1, n * n);

                sw.WriteLine($"0,0,{n * n},{n * n}");

                foreach (MultiPrecision<N8> u in us) {
                    MultiPrecision<N8> q = s * MultiPrecision<N8>.Sqrt(MultiPrecision<N8>.Sqrt(u));

                    (MultiPrecision<N8> a, MultiPrecision<N8> b) = ComputeDigits40(n, q);

                    if (n >= 1) {
                        sw.WriteLine($"{u},{q:e40},{a:e40},{b:e40}");
                        Console.WriteLine($"{u},{a:e20},{b:e20}");
                    }
                    else {
                        sw.WriteLine($"{u},{q:e40},{a:e40},0");
                        Console.WriteLine($"{u},{a:e20}");
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        static (MultiPrecision<N> a, MultiPrecision<N> b) ComputeDigits40<N>(int n, MultiPrecision<N> q) where N : struct, IConstant {
            if (n == 0) { 
                return (EigenMP<N>.Value(EigenFunc.A, n, q, zero_shift: false).value, 0);
            }

            return (EigenMP<N>.Value(EigenFunc.A, n, q, zero_shift: false).value, EigenMP<N>.Value(EigenFunc.B, n, q, zero_shift: false).value);
        }
    }
}