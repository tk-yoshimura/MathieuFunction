using MultiPrecision;

namespace MathieuMP {
    class Program {
        static void Main() {
            List<MultiPrecision<Pow2.N8>> us = new();

            for (MultiPrecision<Pow2.N8> u = MultiPrecision<Pow2.N8>.Ldexp(1, -64); u < 1 / 1024d; u *= 2) {
                us.Add(u);
            }
            for (MultiPrecision<Pow2.N8> u = 1 / 1024d; u < 1; u += 1 / 1024d) {
                us.Add(u);
            }
            for (MultiPrecision<Pow2.N8> u = 1; u < 4; u += 1 / 256d) { 
                us.Add(u);
            }
            for (MultiPrecision<Pow2.N8> u = 4; u < 16; u += 1 / 64d) { 
                us.Add(u);
            }
            for (MultiPrecision<Pow2.N8> u = 16; u < 64; u += 1 / 16d) { 
                us.Add(u);
            }
            for (MultiPrecision<Pow2.N8> u = 64; u < 256; u += 1 / 4d) { 
                us.Add(u);
            }
            for (MultiPrecision<Pow2.N8> u = 256; u <= 1024; u += 1) { 
                us.Add(u);
            }

            for (int n = 0; n <= 64; n++) {
                Console.WriteLine($"Plotting {n}");

                using StreamWriter sw = new($"../../../../results/eigen_precision40_n{n}.csv");
                sw.WriteLine($"# zero shifted mathieu eigen value precision_digits=40 n={n}");
                sw.WriteLine("# u:=q^2/max(1, n^4), mean:=(a+b)/2, diff:=(a-b)/2, scaled_diff:=diff*(2^(n-1)*(n-1)!)^2/q^n");
                sw.WriteLine("u,a,b,mean,diff,scaled_diff,digits_loss(1/0)");

                int s = Math.Max(1, n * n);
                MultiPrecision<Pow2.N8> r = MultiPrecision<Pow2.N8>.Ldexp(MultiPrecision<Pow2.N8>.Square(MultiPrecision<Pow2.N8>.Gamma(n)), 2 * n - 2);

                if (n >= 1) {
                    sw.WriteLine("0,0,0,0,0,1,0");
                }
                else {
                    sw.WriteLine("0,0,0,0,0,0,0");
                }

                foreach (MultiPrecision<Pow2.N8> u in us) {
                    MultiPrecision<Pow2.N8> q = s * MultiPrecision<Pow2.N8>.Sqrt(u);

                    (MultiPrecision<Pow2.N8> a, MultiPrecision<Pow2.N8> b, MultiPrecision<Pow2.N8> m, MultiPrecision<Pow2.N8> d, bool cancellation_digits)
                        = ComputeDigits40(n, q);

                    MultiPrecision<Pow2.N8> sd = r * d / MultiPrecision<Pow2.N8>.Pow(q, n);

                    if (n >= 1) {
                        sw.WriteLine($"{u},{a:e40},{b:e40},{m:e40},{d:e40},{sd:e40},{(cancellation_digits ? "1" : "0")}");
                        Console.WriteLine($"{u},{a:e40},{sd:e40}");
                    }
                    else { 
                        sw.WriteLine($"{u},{a:e40},0,0,0,0,0");
                        Console.WriteLine($"{u},{a:e40}");
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        static (MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> m, MultiPrecision<N> d, bool cancellation_digits) Compute<N, M>(int n, MultiPrecision<M> q, int needs_bits)
            where N: struct, IConstant where M : struct, IConstant {

            if (MultiPrecision<N>.Length > MultiPrecision<M>.Length) {
                throw new ArgumentException("Invalid multi precision length.");
            }
        
            MultiPrecision<M> a = EigenMP<M>.Value(EigenFunc.A, n, q, zero_shift: true).value;
            MultiPrecision<M> b = EigenMP<M>.Value(EigenFunc.B, n, q, zero_shift: true).value;
            MultiPrecision<M> m = (a + b) / 2;
            MultiPrecision<M> d = (a - b) / 2;

            long loss_bits = Math.Max(a.Exponent, b.Exponent) - d.Exponent;
            bool cancellation_digits = loss_bits >= MultiPrecision<M>.Bits - needs_bits;

            return (a.Convert<N>(), b.Convert<N>(), m.Convert<N>(), d.Convert<N>(), cancellation_digits);
        }

        static (MultiPrecision<Pow2.N8> a, MultiPrecision<Pow2.N8> b, MultiPrecision<Pow2.N8> m, MultiPrecision<Pow2.N8> d, bool cancellation_digits) ComputeDigits40(int n, MultiPrecision<Pow2.N8> q) {
            if (n == 0) {
                return (EigenMP<Pow2.N8>.Value(EigenFunc.A, n, q, zero_shift: true).value, 0, 0, 0, cancellation_digits: false);
            }
            
            MultiPrecision<Pow2.N8> a, b, m, d;
            bool cancellation_digits;

            (a, b, m, d, cancellation_digits) = Compute<Pow2.N8, Pow2.N8>(n, q, needs_bits: 144);
            if (!cancellation_digits) {
                return (a, b, m, d, cancellation_digits);
            }

            (a, b, m, d, cancellation_digits) = Compute<Pow2.N8, Pow2.N16>(n, q.Convert<Pow2.N16>(), needs_bits: 144);
            if (!cancellation_digits) {
                return (a, b, m, d, cancellation_digits);
            }

            (a, b, m, d, cancellation_digits) = Compute<Pow2.N8, Pow2.N32>(n, q.Convert<Pow2.N32>(), needs_bits: 144);
            return (a, b, m, d, cancellation_digits);
        }
    }
}