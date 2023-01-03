using MultiPrecision;
using MathieuNearZero;

namespace MathieuMP {
    class Program {
        static void Main() {
            List<MultiPrecision<Pow2.N4>> us = new();

            for (MultiPrecision<Pow2.N4> u = MultiPrecision<Pow2.N4>.Ldexp(1, -32); u <= MultiPrecision<Pow2.N4>.Ldexp(1, 64); u *= 2) {
                us.Add(u);
            }

            for (int n = 0; n <= 64; n++) {
                Console.WriteLine($"Plotting {n}");

                using StreamWriter sw = new($"../../../../results/eigen_nearzero_log2_precision2048bits_n{n}.csv");
                sw.WriteLine($"# zero shifted mathieu eigen value near zero precision_digits=2048bits n={n}");
                sw.WriteLine("# u:=q^2/max(1, n^4), m:=mean/max(1, n^2), d:=1/scaled_diff-1");
                sw.WriteLine("u,m,d,digits_loss(1/0)");

                using BinaryWriter bw = new(File.OpenWrite($"../../../../results/eigen_nearzero_log2_precision2048bits_n{n}.bin"));
 
                int s = Math.Max(1, n * n);
                int i = 0, mp_length = 0;
                
                sw.WriteLine("0,0,0,0");

                foreach (MultiPrecision<Pow2.N4> u in us) {
                    MultiPrecision<Pow2.N64> q = s * MultiPrecision<Pow2.N64>.Sqrt(u.Convert<Pow2.N64>());

                    (MultiPrecision<Pow2.N64> m, MultiPrecision<Pow2.N64> d, bool cancellation_digits, mp_length)
                        = ComputeDigits2048bits(n, q, mp_length);

                    sw.WriteLine($"{u.Convert<Pow2.N16>()},{m},{d},{(cancellation_digits ? "1" : "0")}");
                    bw.Write(u);
                    bw.Write(m);
                    bw.Write(d);
                    Console.WriteLine($"{u}, {m:e40}, {d:e40}");

                    i++;
                    if ((i % 16) == 0) {
                        mp_length = mp_length * 3 / 4;
                    }
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

        static (MultiPrecision<Pow2.N64> m, MultiPrecision<Pow2.N64> d, bool cancellation_digits, int mp_length) ComputeDigits2048bits(int n, MultiPrecision<Pow2.N64> q, int mp_length) {
            if (n == 0) {
                return (EigenMP<Pow2.N64>.Value(EigenFunc.A, n, q, zero_shift: true).value, 0, cancellation_digits: false, mp_length: 0);
            }
            
            MultiPrecision<Pow2.N64> m, d;
            bool cancellation_digits;

            if (mp_length <= 68) {
                (m, d, cancellation_digits) = Compute<Pow2.N64, Plus4<Pow2.N64>>(n, q.Convert<Plus4<Pow2.N64>>(), needs_bits: 2064);
                if (!cancellation_digits) {
                    return (m, d, cancellation_digits, 68);
                }
            }

            if (mp_length <= 72) {
                (m, d, cancellation_digits) = Compute<Pow2.N64, Plus8<Pow2.N64>>(n, q.Convert<Plus8<Pow2.N64>>(), needs_bits: 2064);
                if (!cancellation_digits) {
                    return (m, d, cancellation_digits, 72);
                }
            }

            if (mp_length <= 80) {
                (m, d, cancellation_digits) = Compute<Pow2.N64, Plus16<Pow2.N64>>(n, q.Convert<Plus16<Pow2.N64>>(), needs_bits: 2064);
                if (!cancellation_digits) {
                    return (m, d, cancellation_digits, 80);
                }
            }

            if (mp_length <= 96) {
                (m, d, cancellation_digits) = Compute<Pow2.N64, Plus32<Pow2.N64>>(n, q.Convert<Plus32<Pow2.N64>>(), needs_bits: 2064);
                if (!cancellation_digits) {
                    return (m, d, cancellation_digits, 96);
                }
            }

            (m, d, cancellation_digits) = Compute<Pow2.N64, Pow2.N128>(n, q.Convert<Pow2.N128>(), needs_bits: 2064);
            return (m, d, cancellation_digits, 128);
        }
    }
}