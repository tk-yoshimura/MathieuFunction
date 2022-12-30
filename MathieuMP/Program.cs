using MultiPrecision;

namespace MathieuMP {
    class Program {
        static void Main() {
            foreach (EigenFunc func in new[] { EigenFunc.A, EigenFunc.B }) {
                FracTerms<Pow2.N32>(func);
            }

            Console.WriteLine("END");
            Console.Read();
        }

        private static void FracTerms<N>(EigenFunc func) where N: struct, IConstant {
            int terms_i0 = 0;

            for (int n = func == EigenFunc.A ? 0 : 1; n <= 32; n++) {
                Console.WriteLine($"{func}{n}");

                using StreamWriter sw = new($"../../../../sandbox/needs_frac_log2_mp{MultiPrecision<N>.Length}_{func}_{n}.csv");

                int terms = (terms_i0 <= 0) ? 1 : terms_i0;
                double terms_rate = 1;

                for ((int i, MultiPrecision<N> q) = (0, Math.Max(1, n * n) / 8192d); q <= Math.Max(1, n * n) * 1024d; i++, q *= 2) {
                    int terms_init = Math.Max(1, (int)(terms * terms_rate) - 3);

                    (MultiPrecision<N> a, int terms_res) = EigenMP<N>.ConvergenceFracTerms(func, n, q, terms_init);

                    terms_rate = (i > 0) ? (double)terms_res / terms : 1;
                    terms = terms_res;

                    if (i <= 0) {
                        terms_i0 = terms_res;
                    }

                    Console.WriteLine($"{q},{a},{terms_res},{terms_init},{terms_rate}{((terms_res - terms_init < 1) ? ",*" : "")}");
                    sw.WriteLine($"{q},{a},{terms_res}");
                }
            }
        }
    }
}