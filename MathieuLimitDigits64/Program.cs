using MultiPrecision;
using System.Numerics;
using static MultiPrecision.Pow2;

namespace MathieuMP {
    class Program {
        static void Main() {
            List<MultiPrecision<N8>> us = new();

            for (MultiPrecision<N8> u = 1024; u.Exponent <= 40; u *= 2) {
                for (MultiPrecision<N8> f = 1; f < 2; f += 1d / 16) {
                    us.Add(u * f);
                }
            }

            for (int n = 0; n <= 16; n++) {
                Console.WriteLine($"Plotting {n}");

                using StreamWriter sw = new($"../../../../results/eigen_limit_r2_precision64_n{n}.csv");
                sw.WriteLine($"# zero shifted mathieu eigen limit value precision_digits=64 n={n}");

                if (n >= 1) {
                    sw.WriteLine("# u:=q^2/n^4");
                }
                else {
                    sw.WriteLine("# u:=q^2");
                }

                sw.WriteLine("u,a,a_limit,a_delta,b,b_limit,b_delta,digits_loss(1/0)");

                int s = Math.Max(1, n * n);
                int i = 0, mp_length = 0;

                foreach (MultiPrecision<N8> u in us) {
                    MultiPrecision<N64> q = s * MultiPrecision<N64>.Sqrt(u.Convert<N64>());

                    (MultiPrecision<N8> a, MultiPrecision<N8> a_limit, MultiPrecision<N8> a_delta, MultiPrecision<N8> b, MultiPrecision<N8> b_limit, MultiPrecision<N8> b_delta, bool cancellation_digits, mp_length)
                        = ComputeDigits64(n, q, mp_length);

                    if (n >= 1) {
                        sw.WriteLine($"{u},{a:e64},{a_limit:e64},{a_delta:e64},{b:e64},{b_limit:e64},{b_delta:e64},{(cancellation_digits ? "1" : "0")}");
                        Console.WriteLine($"{u},{a:e20},{a_limit:e20},{b:e20},{b_limit:e20}");
                    }
                    else {
                        sw.WriteLine($"{u},{a:e64},{a_limit:e64},{a_delta:e64},0,0,0,0");
                        Console.WriteLine($"{u},{a:e20},{a_limit:e20}");
                    }

                    sw.Flush();

                    i++;
                    if ((i % 8) == 0) {
                        mp_length -= 4;
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        static (MultiPrecision<N> a, MultiPrecision<N> a_limit, MultiPrecision<N> a_delta, MultiPrecision<N> b, MultiPrecision<N> b_limit, MultiPrecision<N> b_delta, bool cancellation_digits) Compute<N, M>(int n, MultiPrecision<M> q, int needs_bits)
            where N : struct, IConstant where M : struct, IConstant {

            if (MultiPrecision<N>.Length > MultiPrecision<M>.Length) {
                throw new ArgumentException("Invalid multi precision length.");
            }
            MultiPrecision<M> a = EigenMP<M>.Value(EigenFunc.A, n, q, zero_shift: true).value;
            MultiPrecision<M> a_limit = LimitEigenA(n, q);
            MultiPrecision<M> a_delta = a_limit - a;

            MultiPrecision<M> b = (n > 0) ? EigenMP<M>.Value(EigenFunc.B, n, q, zero_shift: true).value : 0;
            MultiPrecision<M> b_limit = (n > 0) ? LimitEigenB(n, q) : 0;
            MultiPrecision<M> b_delta = b_limit - b;

            long loss_bits_a = Math.Max(0, Math.Max(a.Exponent, a_limit.Exponent) - a_delta.Exponent);
            long loss_bits_b = Math.Max(0, Math.Max(b.Exponent, b_limit.Exponent) - b_delta.Exponent);
            long loss_bits = Math.Max(loss_bits_a, loss_bits_b);

            bool cancellation_digits = loss_bits >= MultiPrecision<M>.Bits - needs_bits;

            Console.WriteLine($"loss_bits: {loss_bits}");

            return (a.Convert<N>(), a_limit.Convert<N>(), a_delta.Convert<N>(),
                    b.Convert<N>(), b_limit.Convert<N>(), b_delta.Convert<N>(), cancellation_digits);
        }

        static (MultiPrecision<N8> a, MultiPrecision<N8> a_limit, MultiPrecision<N8> a_delta, MultiPrecision<N8> b, MultiPrecision<N8> b_limit, MultiPrecision<N8> b_delta, bool cancellation_digits, int mp_length) ComputeDigits64(int n, MultiPrecision<N64> q, int mp_length) {
            const int needs_bits = 228;

            MultiPrecision<N8> a, a_limit, a_delta, b, b_limit, b_delta;
            bool cancellation_digits;

            if (mp_length <= 8) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, N8>(n, q.Convert<N8>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 8);
                }
            }

            if (mp_length <= 12) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, Plus4<N8>>(n, q.Convert<Plus4<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 12);
                }
            }

            if (mp_length <= 16) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, Plus8<N8>>(n, q.Convert<Plus8<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 16);
                }
            }

            if (mp_length <= 20) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, Plus12<N8>>(n, q.Convert<Plus12<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 20);
                }
            }

            if (mp_length <= 24) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, Plus16<N8>>(n, q.Convert<Plus16<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 24);
                }
            }

            if (mp_length <= 28) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, Plus20<N8>>(n, q.Convert<Plus20<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 28);
                }
            }

            if (mp_length <= 32) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, Plus24<N8>>(n, q.Convert<Plus24<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 32);
                }
            }

            if (mp_length <= 36) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, Plus28<N8>>(n, q.Convert<Plus28<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 36);
                }
            }

            if (mp_length <= 40) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, Plus32<N8>>(n, q.Convert<Plus32<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 40);
                }
            }

            if (mp_length <= 44) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, Plus36<N8>>(n, q.Convert<Plus36<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 44);
                }
            }

            if (mp_length <= 48) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, Plus40<N8>>(n, q.Convert<Plus40<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 48);
                }
            }

            if (mp_length <= 52) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, Plus44<N8>>(n, q.Convert<Plus44<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 52);
                }
            }

            if (mp_length <= 56) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, Plus48<N8>>(n, q.Convert<Plus48<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 56);
                }
            }

            if (mp_length <= 60) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, Plus52<N8>>(n, q.Convert<Plus52<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 60);
                }
            }

            (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, N64>(n, q.Convert<N64>(), needs_bits);
            if (!cancellation_digits) {
                return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 64);
            }

            (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N8, N128>(n, q.Convert<N128>(), needs_bits);
            return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 64);
        }

        static MultiPrecision<N> LimitEigenA<N>(int n, MultiPrecision<N> q) where N : struct, IConstant {
            int s = 2 * n + 1;
            MultiPrecision<N> y = 2 * (s * MultiPrecision<N>.Sqrt(q) - q) - MultiPrecision<N>.Div(6 * n * n + 2 * n + 1, 4);

            return y;
        }

        static MultiPrecision<N> LimitEigenB<N>(int n, MultiPrecision<N> q) where N : struct, IConstant {
            int s = 2 * n - 1;
            MultiPrecision<N> y = 2 * (s * MultiPrecision<N>.Sqrt(q) - q) - MultiPrecision<N>.Div(6 * n * n - 2 * n + 1, 4);

            return y;
        }
    }
}