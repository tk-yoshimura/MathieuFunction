using MultiPrecision;
using MathieuNearZero;

namespace MathieuMP {
    class Program {
        static void Main() {
            List<MultiPrecision<Pow2.N4>> us = new();

            MultiPrecision<Pow2.N4> h = MultiPrecision<Pow2.N4>.Ldexp(1, -64);

            for (MultiPrecision<Pow2.N4> u = h; u <= h * 256; u += h) {
                us.Add(u);
            }

            for (int n = 1; n <= 64; n++) {
                Console.WriteLine($"Plotting {n}");

                using StreamWriter sw = new($"../../../../results/eigen_nearzero_precision1024bits_n{n}.csv");
                sw.WriteLine($"# zero shifted mathieu eigen value near zero precision_digits=1024bits n={n}");
                sw.WriteLine("# u:=q^2/max(1, n^4), m:=mean/max(1, n^2), d:=1/scaled_diff-1");
                sw.WriteLine("u,m,d,digits_loss(1/0)");

                using BinaryWriter bw = new(File.OpenWrite($"../../../../results/eigen_nearzero_precision1024bits_n{n}.bin"));
 
                int s = Math.Max(1, n * n);
                
                sw.WriteLine("0,0,0,0");

                foreach (MultiPrecision<Pow2.N4> u in us) {
                    MultiPrecision<Pow2.N32> q = s * MultiPrecision<Pow2.N32>.Sqrt(u.Convert<Pow2.N32>());

                    (MultiPrecision<Pow2.N32> m, MultiPrecision<Pow2.N32> d, bool cancellation_digits)
                        = ComputeDigits1024bits(n, q);

                    sw.WriteLine($"{u.Convert<Pow2.N8>()},{m},{d},{(cancellation_digits ? "1" : "0")}");
                    bw.Write(u);
                    bw.Write(m);
                    bw.Write(d);
                    Console.WriteLine($"{u}, {m:e40}, {d:e40}");
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        static (MultiPrecision<N> m, MultiPrecision<N> d, bool cancellation_digits) Compute<N, M>(int n, MultiPrecision<M> q, int needs_bits)
            where N: struct, IConstant where M : struct, IConstant {

            if (MultiPrecision<N>.Length > MultiPrecision<M>.Length) {
                throw new ArgumentException("Invalid multi precision length.");
            }

            static MultiPrecision<M> frac(MultiPrecision<M> n) {
                MultiPrecision<M> v = 1;

                for (int i = 2; i <= n; i++) {
                    v *= i;
                }

                return v;
            };

            MultiPrecision<M> r = MultiPrecision<M>.Ldexp(MultiPrecision<M>.Square(frac(n - 1)), 2 * n - 2) / MultiPrecision<M>.Pow(q, n);
        
            MultiPrecision<M> a = EigenMP<M>.Value(EigenFunc.A, n, q, zero_shift: true).value;
            MultiPrecision<M> b = EigenMP<M>.Value(EigenFunc.B, n, q, zero_shift: true).value;
            MultiPrecision<M> m = (a + b) / (2 * n * n);

            MultiPrecision<M> e = (a - b) / 2;
            MultiPrecision<M> g = 1 / (r * e);
            MultiPrecision<M> d = g - 1;

            long loss_bits = Math.Max(0, Math.Max(a.Exponent, b.Exponent) - e.Exponent) + Math.Max(0, g.Exponent - d.Exponent);
            bool cancellation_digits = loss_bits >= MultiPrecision<M>.Bits - needs_bits;

            Console.WriteLine($"loss_bits: {loss_bits}");

            return (m.Convert<N>(), d.Convert<N>(), cancellation_digits);
        }

        static (MultiPrecision<Pow2.N32> m, MultiPrecision<Pow2.N32> d, bool cancellation_digits) ComputeDigits1024bits(int n, MultiPrecision<Pow2.N32> q) {
            if (n == 0) {
                return (EigenMP<Pow2.N32>.Value(EigenFunc.A, n, q, zero_shift: true).value, 0, cancellation_digits: false);
            }
            
            MultiPrecision<Pow2.N32> m, d;
            bool cancellation_digits;

            (m, d, cancellation_digits) = Compute<Pow2.N32, Plus4<Pow2.N32>>(n, q.Convert<Plus4<Pow2.N32>>(), needs_bits: 1040);
            if (!cancellation_digits) {
                return (m, d, cancellation_digits);
            }

            (m, d, cancellation_digits) = Compute<Pow2.N32, Plus8<Pow2.N32>>(n, q.Convert<Plus8<Pow2.N32>>(), needs_bits: 1040);
            if (!cancellation_digits) {
                return (m, d, cancellation_digits);
            }

            (m, d, cancellation_digits) = Compute<Pow2.N32, Plus16<Pow2.N32>>(n, q.Convert<Plus16<Pow2.N32>>(), needs_bits: 1040);
            if (!cancellation_digits) {
                return (m, d, cancellation_digits);
            }

            (m, d, cancellation_digits) = Compute<Pow2.N32, Pow2.N64>(n, q.Convert<Pow2.N64>(), needs_bits: 1040);
            if (!cancellation_digits) {
                return (m, d, cancellation_digits);
            }

            (m, d, cancellation_digits) = Compute<Pow2.N32, Plus8<Pow2.N64>>(n, q.Convert<Plus8<Pow2.N64>>(), needs_bits: 1040);
            if (!cancellation_digits) {
                return (m, d, cancellation_digits);
            }

            (m, d, cancellation_digits) = Compute<Pow2.N32, Plus16<Pow2.N64>>(n, q.Convert<Plus16<Pow2.N64>>(), needs_bits: 1040);
            if (!cancellation_digits) {
                return (m, d, cancellation_digits);
            }

            (m, d, cancellation_digits) = Compute<Pow2.N32, Plus32<Pow2.N64>>(n, q.Convert<Plus32<Pow2.N64>>(), needs_bits: 1040);
            if (!cancellation_digits) {
                return (m, d, cancellation_digits);
            }

            (m, d, cancellation_digits) = Compute<Pow2.N32, Pow2.N128>(n, q.Convert<Pow2.N128>(), needs_bits: 1040);
            return (m, d, cancellation_digits);
        }
    }
}