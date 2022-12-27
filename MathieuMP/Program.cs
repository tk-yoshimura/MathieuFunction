using MultiPrecision;

namespace MathieuMP {
    class Program {
        static void Main() {
            foreach (EigenFunc func in new[] { EigenFunc.A, EigenFunc.B }) {
                FracTerms<Pow2.N4>(func);
            }
            foreach (EigenFunc func in new[] { EigenFunc.A, EigenFunc.B }) {
                FracTerms<Pow2.N8>(func);
            }
            foreach (EigenFunc func in new[] { EigenFunc.A, EigenFunc.B }) {
                FracTerms<Pow2.N16>(func);
            }
            foreach (EigenFunc func in new[] { EigenFunc.A, EigenFunc.B }) {
                FracTerms<Pow2.N32>(func);
            }

            Console.WriteLine("END");
            Console.Read();
        }

        private static void FracTerms<N>(EigenFunc func) where N: struct, IConstant {
            List<int> terms_list = new();

            for (int n = func == EigenFunc.A ? 0 : 1; n <= 256; n = (n < 64) ? n + 1 : n * 2) {
                Console.WriteLine($"{func}{n}");

                using StreamWriter sw = new($"../../../../sandbox/needs_frac_log2_mp{MultiPrecision<N>.Length}_{func}_{n}.csv");

                int terms = 1, terms_diff = 0;

                for ((int i, MultiPrecision<N> q) = (0, Math.Max(1, n * n) / 4d); q <= Math.Max(1, n * n) * 256d; i++, q *= 2) {
                    terms = (i < terms_list.Count) ? Math.Max(terms, terms_list[i]) : terms;

                    (MultiPrecision<N> a, int terms_res) = EigenMP<N>.ConvergenceFracTerms(func, n, q, Math.Max(1, terms + terms_diff - 2));

                    terms_diff = Math.Max(0, terms_res - terms);
                    terms = terms_res;
 
                    if (i < terms_list.Count) {
                        terms_list[i] = terms_res;
                    }
                    else {
                        terms_list.Add(terms_res);
                    }

                    Console.WriteLine($"{q},{a},{terms_res},{terms_diff}");
                    sw.WriteLine($"{q},{a},{terms_res}");
                }
            }
        }
    }
}