using MultiPrecision;
using System.Numerics;

namespace MathieuMP {
    class Program {
        static void Main() {
            List<MultiPrecision<Pow2.N4>> us = [];

            for (MultiPrecision<Pow2.N4> u = MultiPrecision<Pow2.N4>.Ldexp(1, -80); u <= 256; u *= 2) {
                us.Add(u);
            }

            for (int n = 0; n <= 64; n++) {
                Console.WriteLine($"Plotting {n}");

                using StreamWriter sw = new($"../../../../results/eigen_nearzero_log2_precision2560bits_n{n}.csv");
                sw.WriteLine($"# zero shifted mathieu eigen value near zero precision_digits=2560bits n={n}");

                if (n >= 1) {
                    sw.WriteLine("# u:=q^2/n^4, m:=(a+b-2n^2)/2n^2, d:=q^n/(2^(n-1)(n-1)!)^2 2/(a-b)-1");
                }
                else {
                    sw.WriteLine("# u:=q^2, m:=a, d:=0");
                }

                sw.WriteLine("u,m,d,digits_loss(1/0)");

                using BinaryWriter bw = new(File.OpenWrite($"../../../../results/eigen_nearzero_log2_precision2560bits_n{n}.bin"));

                int s = Math.Max(1, n * n);
                int i = 0, mp_length = 0;

                sw.WriteLine("0,0,0,0");

                foreach (MultiPrecision<Pow2.N4> u in us) {
                    MultiPrecision<Pow2.N128> q = s * MultiPrecision<Pow2.N128>.Sqrt(u.Convert<Pow2.N128>());

                    (MultiPrecision<N80> m, MultiPrecision<N80> d, bool cancellation_digits, mp_length)
                        = ComputeDigits2560bits(n, q, mp_length);

                    sw.WriteLine($"{u.Convert<Pow2.N16>()},{m},{d},{(cancellation_digits ? "1" : "0")}");
                    bw.Write(u);
                    bw.Write(m);
                    bw.Write(d);
                    bw.Flush();

                    Console.WriteLine($"{u}, {m:e40}, {d:e40}");

                    i++;
                    if ((i % 16) == 0) {
                        mp_length -= 8;
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        static (MultiPrecision<N> m, MultiPrecision<N> d, bool cancellation_digits) Compute<N, M>(int n, MultiPrecision<M> q, int needs_bits)
            where N : struct, IConstant where M : struct, IConstant {

            if (MultiPrecision<N>.Length > MultiPrecision<M>.Length) {
                throw new ArgumentException("Invalid multi precision length.");
            }

            static BigInteger frac(long n) {
                BigInteger v = 1;

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

        static (MultiPrecision<N80> m, MultiPrecision<N80> d, bool cancellation_digits, int mp_length) ComputeDigits2560bits(int n, MultiPrecision<Pow2.N128> q, int mp_length) {
            if (n == 0) {
                return (EigenMP<N80>.Value(EigenFunc.A, n, q.Convert<N80>(), zero_shift: true).value, 0, cancellation_digits: false, mp_length: 0);
            }

            int needs_bits = MultiPrecision<N80>.Bits + 16;

            MultiPrecision<N80> m, d;
            bool cancellation_digits;

            if (mp_length <= 84) {
                (m, d, cancellation_digits) = Compute<N80, Plus4<N80>>(n, q.Convert<Plus4<N80>>(), needs_bits);
                if (!cancellation_digits) {
                    return (m, d, cancellation_digits, 84);
                }
            }

            if (mp_length <= 88) {
                (m, d, cancellation_digits) = Compute<N80, Plus8<N80>>(n, q.Convert<Plus8<N80>>(), needs_bits);
                if (!cancellation_digits) {
                    return (m, d, cancellation_digits, 88);
                }
            }

            if (mp_length <= 96) {
                (m, d, cancellation_digits) = Compute<N80, Plus16<N80>>(n, q.Convert<Plus16<N80>>(), needs_bits);
                if (!cancellation_digits) {
                    return (m, d, cancellation_digits, 96);
                }
            }

            if (mp_length <= 104) {
                (m, d, cancellation_digits) = Compute<N80, Plus24<N80>>(n, q.Convert<Plus24<N80>>(), needs_bits);
                if (!cancellation_digits) {
                    return (m, d, cancellation_digits, 104);
                }
            }

            if (mp_length <= 112) {
                (m, d, cancellation_digits) = Compute<N80, Plus32<N80>>(n, q.Convert<Plus32<N80>>(), needs_bits);
                if (!cancellation_digits) {
                    return (m, d, cancellation_digits, 112);
                }
            }

            if (mp_length <= 120) {
                (m, d, cancellation_digits) = Compute<N80, Plus40<N80>>(n, q.Convert<Plus40<N80>>(), needs_bits);
                if (!cancellation_digits) {
                    return (m, d, cancellation_digits, 120);
                }
            }

            (m, d, cancellation_digits) = Compute<N80, Pow2.N128>(n, q.Convert<Pow2.N128>(), needs_bits);
            return (m, d, cancellation_digits, 128);
        }
    }
}