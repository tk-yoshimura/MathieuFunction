using MultiPrecision;

namespace MathieuMP {
    internal static class ForwardFiniteDifference<N, M> where N : struct, IConstant where M : struct, IConstant {
        public static readonly MultiPrecision<M>[][] table;

        public static int SamplePoints { get; private set; }

        static ForwardFiniteDifference() {
            using StreamReader sr = new("../../../../sandbox/forward_log2pts_acc80.md");
            sr.ReadLine();
            sr.ReadLine();

            string[] lines = sr.ReadToEnd().Split("\n").Where((line) => !string.IsNullOrWhiteSpace(line)).ToArray();

            table = new MultiPrecision<M>[lines.Length][];

            for (int i = 0; i < lines.Length; i++) {
                string[] fracs = lines[i].Split('|').Where((item) => !string.IsNullOrWhiteSpace(item)).Skip(2).ToArray();

                table[i] = new MultiPrecision<M>[fracs.Length];

                for (int j = 0; j < fracs.Length; j++) {
                    Fraction f = fracs[j];

                    table[i][j] = MultiPrecision<M>.Div(f.Numer, f.Denom);
                }

            }

            SamplePoints = table[^1].Length;
        }

        public static MultiPrecision<N>[] Diff(MultiPrecision<N>[] xs, MultiPrecision<N> h) {
            if (xs.Length != SamplePoints) {
                throw new ArgumentException(nameof(xs));
            }

            MultiPrecision<N>[] ds = new MultiPrecision<N>[table.Length];

            for (int i = 0; i < table.Length; i++) {
                MultiPrecision<M> d = 0;
                MultiPrecision<M>[] ts = table[i];

                for (int j = 0; j < ts.Length; j++) {
                    if (ts[j].IsZero || xs[j].IsZero) {
                        continue;
                    }

                    d += ts[j] * xs[j].Convert<M>();
                }

                d /= MultiPrecision<M>.Pow(h.Convert<M>(), i + 1);

                ds[i] = d.Convert<N>();
            }

            return ds;
        }
    }
}
