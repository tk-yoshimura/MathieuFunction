using MultiPrecision;
using System.Numerics;
using static MultiPrecision.Pow2;

namespace MathieuMP {
    class Program {
        static void Main() {
            List<MultiPrecision<N8>> us = new();

            for (MultiPrecision<N8> u = MultiPrecision<N8>.Ldexp(1, -64); u < 1 / 1024d; u *= 2) {
                us.Add(u);
            }
            for (MultiPrecision<N8> u = 1 / 1024d; u < 1; u += 1 / 1024d) {
                us.Add(u);
            }
            for (MultiPrecision<N8> u = 1; u < 4; u += 1 / 256d) {
                us.Add(u);
            }
            for (MultiPrecision<N8> u = 4; u < 16; u += 1 / 64d) {
                us.Add(u);
            }
            for (MultiPrecision<N8> u = 16; u < 64; u += 1 / 16d) {
                us.Add(u);
            }
            for (MultiPrecision<N8> u = 64; u < 256; u += 1 / 4d) {
                us.Add(u);
            }
            for (MultiPrecision<N8> u = 256; u <= 1024; u += 1) {
                us.Add(u);
            }

            for (int n = 37; n <= 64; n++) {
                Console.WriteLine($"Plotting {n}");

                using StreamWriter sw = new($"../../../../results/eigen_precision64_n{n}.csv");
                sw.WriteLine($"# zero shifted mathieu eigen value precision_digits=64 n={n}");

                if (n >= 1) {
                    sw.WriteLine("# u:=q^2/n^4, m:=(a+b-2n^2)/2n^2, d:=q^n/(2^(n-1)(n-1)!)^2 2/(a-b)-1");
                }
                else {
                    sw.WriteLine("# u:=q^2, m:=a, d:=0");
                }

                sw.WriteLine("u,a,b,m,d,digits_loss(1/0)");

                int s = Math.Max(1, n * n);
                int i = 0, mp_length = 0;

                sw.WriteLine("0,0,0,0,0,0");

                foreach (MultiPrecision<N8> u in us) {
                    MultiPrecision<N64> q = s * MultiPrecision<N64>.Sqrt(u.Convert<N64>());

                    (MultiPrecision<N8> a, MultiPrecision<N8> b, MultiPrecision<N8> m, MultiPrecision<N8> d, bool cancellation_digits, mp_length)
                        = ComputeDigits64(n, q, mp_length);

                    if (n >= 1) {
                        sw.WriteLine($"{u},{a:e64},{b:e64},{m:e64},{d:e64},{(cancellation_digits ? "1" : "0")}");
                        Console.WriteLine($"{u},{a:e20},{d:e20}");
                    }
                    else {
                        sw.WriteLine($"{u},{a:e64},0,{a:e64},0,0,0");
                        Console.WriteLine($"{u},{a:e20}");
                    }

                    i++;
                    if ((i % 8) == 0) {
                        mp_length -= 4;
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        static (MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> m, MultiPrecision<N> d, bool cancellation_digits) Compute<N, M>(int n, MultiPrecision<M> q, int needs_bits)
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

            return (a.Convert<N>(), b.Convert<N>(), m.Convert<N>(), d.Convert<N>(), cancellation_digits);
        }

        static (MultiPrecision<N8> a, MultiPrecision<N8> b, MultiPrecision<N8> m, MultiPrecision<N8> d, bool cancellation_digits, int mp_length) ComputeDigits64(int n, MultiPrecision<N64> q, int mp_length) {
            if (n == 0) {
                return (EigenMP<N8>.Value(EigenFunc.A, n, q.Convert<N8>(), zero_shift: true).value, 0, 0, 0, cancellation_digits: false, mp_length: 0);
            }

            const int needs_bits = 228;

            MultiPrecision<N8> a, b, m, d;
            bool cancellation_digits;

            if (mp_length <= 8) {
                (a, b, m, d, cancellation_digits) = Compute<N8, N8>(n, q.Convert<N8>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, m, d, cancellation_digits, 8);
                }
            }

            if (mp_length <= 12) {
                (a, b, m, d, cancellation_digits) = Compute<N8, Plus4<N8>>(n, q.Convert<Plus4<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, m, d, cancellation_digits, 12);
                }
            }

            if (mp_length <= 16) {
                (a, b, m, d, cancellation_digits) = Compute<N8, Plus8<N8>>(n, q.Convert<Plus8<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, m, d, cancellation_digits, 16);
                }
            }

            if (mp_length <= 20) {
                (a, b, m, d, cancellation_digits) = Compute<N8, Plus12<N8>>(n, q.Convert<Plus12<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, m, d, cancellation_digits, 20);
                }
            }

            if (mp_length <= 24) {
                (a, b, m, d, cancellation_digits) = Compute<N8, Plus16<N8>>(n, q.Convert<Plus16<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, m, d, cancellation_digits, 24);
                }
            }

            if (mp_length <= 28) {
                (a, b, m, d, cancellation_digits) = Compute<N8, Plus20<N8>>(n, q.Convert<Plus20<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, m, d, cancellation_digits, 28);
                }
            }

            if (mp_length <= 32) {
                (a, b, m, d, cancellation_digits) = Compute<N8, Plus24<N8>>(n, q.Convert<Plus24<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, m, d, cancellation_digits, 32);
                }
            }

            if (mp_length <= 36) {
                (a, b, m, d, cancellation_digits) = Compute<N8, Plus28<N8>>(n, q.Convert<Plus28<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, m, d, cancellation_digits, 36);
                }
            }

            if (mp_length <= 40) {
                (a, b, m, d, cancellation_digits) = Compute<N8, Plus32<N8>>(n, q.Convert<Plus32<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, m, d, cancellation_digits, 40);
                }
            }

            if (mp_length <= 44) {
                (a, b, m, d, cancellation_digits) = Compute<N8, Plus36<N8>>(n, q.Convert<Plus36<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, m, d, cancellation_digits, 44);
                }
            }

            if (mp_length <= 48) {
                (a, b, m, d, cancellation_digits) = Compute<N8, Plus40<N8>>(n, q.Convert<Plus40<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, m, d, cancellation_digits, 48);
                }
            }

            if (mp_length <= 52) {
                (a, b, m, d, cancellation_digits) = Compute<N8, Plus44<N8>>(n, q.Convert<Plus44<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, m, d, cancellation_digits, 52);
                }
            }

            if (mp_length <= 56) {
                (a, b, m, d, cancellation_digits) = Compute<N8, Plus48<N8>>(n, q.Convert<Plus48<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, m, d, cancellation_digits, 56);
                }
            }

            if (mp_length <= 60) {
                (a, b, m, d, cancellation_digits) = Compute<N8, Plus52<N8>>(n, q.Convert<Plus52<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, m, d, cancellation_digits, 60);
                }
            }

            (a, b, m, d, cancellation_digits) = Compute<N8, N64>(n, q.Convert<N64>(), needs_bits);
            if (!cancellation_digits) {
                return (a, b, m, d, cancellation_digits, 64);
            }

            (a, b, m, d, cancellation_digits) = Compute<N8, N128>(n, q.Convert<N128>(), needs_bits);
            return (a, b, m, d, cancellation_digits, 64);
        }
    }
}