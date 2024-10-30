using MultiPrecision;
using static MultiPrecision.Pow2;

namespace MathieuMP {
    class Program {
        static void Main() {
            List<MultiPrecision<N16>> fs = [];

            for (MultiPrecision<N16> f = MultiPrecision<N16>.Ldexp(1, -40); f < MultiPrecision<N16>.Ldexp(1, -20); f *= 2) {
                for (MultiPrecision<N16> r = 1; r < 2; r += 1d / 64) {
                    fs.Add(f * r);
                }
            }
            for (MultiPrecision<N16> f = MultiPrecision<N16>.Ldexp(1, -20); f < MultiPrecision<N16>.Ldexp(1, -16); f *= 2) {
                for (MultiPrecision<N16> r = 1; r < 2; r += 1d / 256) {
                    fs.Add(f * r);
                }
            }
            for (MultiPrecision<N16> f = MultiPrecision<N16>.Ldexp(1, -16); f <= MultiPrecision<N16>.Ldexp(1, -12); f *= 2) {
                for (MultiPrecision<N16> r = 1; r < 2; r += 1d / 256) {
                    fs.Add(f * r);
                }
            }

            for (int n = 14; n <= 15; n++) {
                Console.WriteLine($"Plotting {n}");

                using StreamWriter sw = new($"../../../../results/eigen_limit_precision145_n{n}.csv");
                sw.WriteLine($"# zero shifted mathieu eigen limit value precision_digits=145 n={n}");

                if (n >= 1) {
                    sw.WriteLine("# u:=q^2/n^4");
                }
                else {
                    sw.WriteLine("# u:=q^2");
                }

                sw.WriteLine("1/u,a,a_limit,a_delta,b,b_limit,b_delta,digits_loss(1/0)");

                int s = Math.Max(1, n * n);
                int i = 0, mp_length = 0;

                sw.WriteLine("0,-inf,-inf,0,-inf,-inf,0,0");

                foreach (MultiPrecision<N16> f in fs) {
                    MultiPrecision<N64> u = 1 / f.Convert<N64>();
                    MultiPrecision<N64> q = s * MultiPrecision<N64>.Sqrt(u);

                    (MultiPrecision<N16> a, MultiPrecision<N16> a_limit, MultiPrecision<N16> a_delta, MultiPrecision<N16> b, MultiPrecision<N16> b_limit, MultiPrecision<N16> b_delta, bool cancellation_digits, mp_length)
                        = ComputeDigits64(n, q, mp_length);

                    if (n >= 1) {
                        sw.WriteLine($"{f},{a},{a_limit},{a_delta},{b},{b_limit},{b_delta},{(cancellation_digits ? "1" : "0")}");
                        Console.WriteLine($"{f},{a:e20},{a_limit:e20},{b:e20},{b_limit:e20}");
                    }
                    else {
                        sw.WriteLine($"{f},{a},{a_limit},{a_delta},0,0,0,0");
                        Console.WriteLine($"{f},{a:e20},{a_limit:e20}");
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

        static (MultiPrecision<N16> a, MultiPrecision<N16> a_limit, MultiPrecision<N16> a_delta, MultiPrecision<N16> b, MultiPrecision<N16> b_limit, MultiPrecision<N16> b_delta, bool cancellation_digits, int mp_length) ComputeDigits64(int n, MultiPrecision<N64> q, int mp_length) {
            const int needs_bits = 484;

            MultiPrecision<N16> a, a_limit, a_delta, b, b_limit, b_delta;
            bool cancellation_digits;

            if (mp_length <= 8) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N16, N16>(n, q.Convert<N16>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 8);
                }
            }

            if (mp_length <= 12) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N16, Plus4<N16>>(n, q.Convert<Plus4<N16>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 12);
                }
            }

            if (mp_length <= 16) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N16, Plus8<N16>>(n, q.Convert<Plus8<N16>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 16);
                }
            }

            if (mp_length <= 20) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N16, Plus12<N16>>(n, q.Convert<Plus12<N16>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 20);
                }
            }

            if (mp_length <= 24) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N16, Plus16<N16>>(n, q.Convert<Plus16<N16>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 24);
                }
            }

            if (mp_length <= 28) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N16, Plus20<N16>>(n, q.Convert<Plus20<N16>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 28);
                }
            }

            if (mp_length <= 32) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N16, Plus24<N16>>(n, q.Convert<Plus24<N16>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 32);
                }
            }

            if (mp_length <= 36) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N16, Plus28<N16>>(n, q.Convert<Plus28<N16>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 36);
                }
            }

            if (mp_length <= 40) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N16, Plus32<N16>>(n, q.Convert<Plus32<N16>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 40);
                }
            }

            if (mp_length <= 44) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N16, Plus36<N16>>(n, q.Convert<Plus36<N16>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 44);
                }
            }

            if (mp_length <= 48) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N16, Plus40<N16>>(n, q.Convert<Plus40<N16>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 48);
                }
            }

            if (mp_length <= 52) {
                (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N16, Plus44<N16>>(n, q.Convert<Plus44<N16>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 52);
                }
            }

            (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N16, N64>(n, q.Convert<N64>(), needs_bits);
            if (!cancellation_digits) {
                return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 64);
            }

            (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits) = Compute<N16, N128>(n, q.Convert<N128>(), needs_bits);
            return (a, a_limit, a_delta, b, b_limit, b_delta, cancellation_digits, 64);
        }

        static MultiPrecision<N> LimitEigenA<N>(int n, MultiPrecision<N> q) where N : struct, IConstant {
            int s = 2 * n + 1;
            MultiPrecision<N> h = MultiPrecision<N>.Sqrt(q);

            MultiPrecision<N> y = 2 * (s * h - q) - MultiPrecision<N>.Div(6 * n * n + 2 * n + 1, 4) - DeltaTerm5(s, h);

            return y;
        }

        static MultiPrecision<N> LimitEigenB<N>(int n, MultiPrecision<N> q) where N : struct, IConstant {
            int s = 2 * n - 1;
            MultiPrecision<N> h = MultiPrecision<N>.Sqrt(q);

            MultiPrecision<N> y = 2 * (s * h - q) - MultiPrecision<N>.Div(6 * n * n - 2 * n + 1, 4) - DeltaTerm5(s, h);

            return y;
        }

        static MultiPrecision<N> DeltaTerm5<N>(long s, MultiPrecision<N> h) where N : struct, IConstant {
            MultiPrecision<N> delta = 0;

            delta += (s * (3 + s * s)) / (128 * h);
            delta += (9 + s * s * (34 + s * s * 5)) / (4096 * h * h);
            delta += (s * (405 + s * s * (410 + s * s * 33))) / (131072 * h * h * h);
            delta += (486 + s * s * (2943 + s * s * (1260 + s * s * 63))) / (1048576 * h * h * h * h);
            delta += (checked(s * (41607 + s * s * (69001 + s * s * (15617 + s * s * 527))))) / (33554432 * h * h * h * h * h);

            return delta;
        }
    }
}