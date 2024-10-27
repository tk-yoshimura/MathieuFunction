using MultiPrecision;
using System.Numerics;
using static MultiPrecision.Pow2;

namespace MathieuMP {
    class Program {
        static void Main() {
            List<MultiPrecision<N8>> us = [];

            for (MultiPrecision<N8> u = MultiPrecision<N8>.Ldexp(1, -64); u < 1 / 16384d; u *= 2) {
                us.Add(u);
            }
            for (MultiPrecision<N8> q = 1 / 16384d; q <= 1 / 4d; q += 1 / 16384d) {
                us.Add(q);
            }

            for (int n = 1; n <= 16; n++) {
                Console.WriteLine($"Plotting {n}");

                using StreamWriter sw = new($"../../../../results/eigen_nearzero_r2_precision64_n{n}.csv");
                sw.WriteLine($"# zero shifted mathieu eigen value precision_digits=64 n={n}");
                sw.WriteLine("# u:=q/(n^2)");

                if (n == 1) {
                    sw.WriteLine("# az:=a-(1+q), bz:=b-(1-q)");
                }
                else {
                    sw.WriteLine("# az:=a-n^2, bz:=b-n^2");
                }

                sw.WriteLine("u,q,a,b,az,bz,az/q^2,bz/q^2,digits_loss(1/0)");

                int s = Math.Max(1, n * n);
                int i = 0, mp_length = 0;

                (MultiPrecision<N8> azq2_0, MultiPrecision<N8> bzq2_0) = n switch {
                    1 => (-MultiPrecision<N8>.Rcp(8), -MultiPrecision<N8>.Rcp(8)),
                    2 => (MultiPrecision<N8>.Div(5, 12), -MultiPrecision<N8>.Div(1, 12)),
                    _ => (MultiPrecision<N8>.Rcp(2 * (n * n - 1)), MultiPrecision<N8>.Rcp(2 * (n * n - 1))),
                };

                sw.WriteLine($"0,0,0,0,0,0,{azq2_0:e64},{bzq2_0:e64},0");

                foreach (MultiPrecision<N8> u in us) {
                    MultiPrecision<N64> q = s * u.Convert<N64>();

                    (MultiPrecision<N8> a, MultiPrecision<N8> b, MultiPrecision<N8> az, MultiPrecision<N8> bz, bool cancellation_digits, mp_length)
                        = ComputeDigits64(n, q, mp_length);

                    sw.WriteLine($"{u},{q},{a:e64},{b:e64},{az:e64},{bz:e64},{(az/(q*q).Convert<N8>()):e64},{(bz/(q*q).Convert<N8>()):e64},{(cancellation_digits ? "1" : "0")}");
                    Console.WriteLine($"{u},{q},{a:e20},{az:e20}");

                    i++;
                    if ((i % 8) == 0) {
                        mp_length -= 4;
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        static (MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> az, MultiPrecision<N> bz, bool cancellation_digits) Compute<N, M>(int n, MultiPrecision<M> q, int needs_bits)
            where N : struct, IConstant where M : struct, IConstant {

            if (MultiPrecision<N>.Length > MultiPrecision<M>.Length) {
                throw new ArgumentException("Invalid multi precision length.");
            }

            MultiPrecision<M> a = EigenMP<M>.Value(EigenFunc.A, n, q, zero_shift: true).value;
            MultiPrecision<M> b = EigenMP<M>.Value(EigenFunc.B, n, q, zero_shift: true).value;

            MultiPrecision<M> az = (n != 1) ? a : (a - q);
            MultiPrecision<M> bz = (n != 1) ? b : (b + q);

            long loss_bits_a = Math.Max(0, Math.Max(a.Exponent, q.Exponent) - az.Exponent);
            long loss_bits_b = Math.Max(0, Math.Max(b.Exponent, q.Exponent) - bz.Exponent);
            long loss_bits = Math.Max(loss_bits_a, loss_bits_b);

            bool cancellation_digits = loss_bits >= MultiPrecision<M>.Bits - needs_bits;

            Console.WriteLine($"loss_bits: {loss_bits}");

            return (a.Convert<N>(), b.Convert<N>(), az.Convert<N>(), bz.Convert<N>(), cancellation_digits);
        }

        static (MultiPrecision<N8> a, MultiPrecision<N8> b, MultiPrecision<N8> az, MultiPrecision<N8> bz, bool cancellation_digits, int mp_length) ComputeDigits64(int n, MultiPrecision<N64> q, int mp_length) {
            if (n == 0) {
                return (EigenMP<N8>.Value(EigenFunc.A, n, q.Convert<N8>(), zero_shift: true).value, 0, 0, 0, cancellation_digits: false, mp_length: 0);
            }

            const int needs_bits = 228;

            MultiPrecision<N8> a, b, az, bz;
            bool cancellation_digits;

            if (mp_length <= 8) {
                (a, b, az, bz, cancellation_digits) = Compute<N8, N8>(n, q.Convert<N8>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, az, bz, cancellation_digits, 8);
                }
            }

            if (mp_length <= 12) {
                (a, b, az, bz, cancellation_digits) = Compute<N8, Plus4<N8>>(n, q.Convert<Plus4<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, az, bz, cancellation_digits, 12);
                }
            }

            if (mp_length <= 16) {
                (a, b, az, bz, cancellation_digits) = Compute<N8, Plus8<N8>>(n, q.Convert<Plus8<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, az, bz, cancellation_digits, 16);
                }
            }

            if (mp_length <= 20) {
                (a, b, az, bz, cancellation_digits) = Compute<N8, Plus12<N8>>(n, q.Convert<Plus12<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, az, bz, cancellation_digits, 20);
                }
            }

            if (mp_length <= 24) {
                (a, b, az, bz, cancellation_digits) = Compute<N8, Plus16<N8>>(n, q.Convert<Plus16<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, az, bz, cancellation_digits, 24);
                }
            }

            if (mp_length <= 28) {
                (a, b, az, bz, cancellation_digits) = Compute<N8, Plus20<N8>>(n, q.Convert<Plus20<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, az, bz, cancellation_digits, 28);
                }
            }

            if (mp_length <= 32) {
                (a, b, az, bz, cancellation_digits) = Compute<N8, Plus24<N8>>(n, q.Convert<Plus24<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, az, bz, cancellation_digits, 32);
                }
            }

            if (mp_length <= 36) {
                (a, b, az, bz, cancellation_digits) = Compute<N8, Plus28<N8>>(n, q.Convert<Plus28<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, az, bz, cancellation_digits, 36);
                }
            }

            if (mp_length <= 40) {
                (a, b, az, bz, cancellation_digits) = Compute<N8, Plus32<N8>>(n, q.Convert<Plus32<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, az, bz, cancellation_digits, 40);
                }
            }

            if (mp_length <= 44) {
                (a, b, az, bz, cancellation_digits) = Compute<N8, Plus36<N8>>(n, q.Convert<Plus36<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, az, bz, cancellation_digits, 44);
                }
            }

            if (mp_length <= 48) {
                (a, b, az, bz, cancellation_digits) = Compute<N8, Plus40<N8>>(n, q.Convert<Plus40<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, az, bz, cancellation_digits, 48);
                }
            }

            if (mp_length <= 52) {
                (a, b, az, bz, cancellation_digits) = Compute<N8, Plus44<N8>>(n, q.Convert<Plus44<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, az, bz, cancellation_digits, 52);
                }
            }

            if (mp_length <= 56) {
                (a, b, az, bz, cancellation_digits) = Compute<N8, Plus48<N8>>(n, q.Convert<Plus48<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, az, bz, cancellation_digits, 56);
                }
            }

            if (mp_length <= 60) {
                (a, b, az, bz, cancellation_digits) = Compute<N8, Plus52<N8>>(n, q.Convert<Plus52<N8>>(), needs_bits);
                if (!cancellation_digits) {
                    return (a, b, az, bz, cancellation_digits, 60);
                }
            }

            (a, b, az, bz, cancellation_digits) = Compute<N8, N64>(n, q.Convert<N64>(), needs_bits);
            if (!cancellation_digits) {
                return (a, b, az, bz, cancellation_digits, 64);
            }

            (a, b, az, bz, cancellation_digits) = Compute<N8, N128>(n, q.Convert<N128>(), needs_bits);
            return (a, b, az, bz, cancellation_digits, 64);
        }
    }
}