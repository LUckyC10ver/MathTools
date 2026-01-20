using System;
using System.Collections.Generic;
using System.IO;
using MathTools.Core;

namespace MathTools.Console
{
    internal static class OptimizeTests
    {
        public static Dictionary<string, Action> BuildRegistry()
        {
            return new Dictionary<string, Action>(StringComparer.OrdinalIgnoreCase)
            {
                { "opt1d", TestOpt1DMin },
                { "linprog", TestLinprog },
                { "quadprog", TestQuadprog },
                { "sqp", TestSqp },
                { "bnb", TestBranchAndBound }
            };
        }

        public static void RunAll(Dictionary<string, Action> registry)
        {
            foreach (var item in registry)
            {
                RunSingle(item.Key, item.Value);
            }
        }

        public static void RunSingle(string name, Action test)
        {
            System.Console.WriteLine($"== Running {name} ==");
            try
            {
                test();
                System.Console.WriteLine($"[PASS] {name}");
            }
            catch (Exception ex)
            {
                System.Console.WriteLine($"[FAIL] {name}: {ex.Message}");
            }
            System.Console.WriteLine();
        }

        private static void TestOpt1DMin()
        {
            double a = 0.0 * Math.PI;
            double b = 0.5 * Math.PI;
            double c = 2.0 * Math.PI;
            double fa;
            double fb;
            double fc;
            double xmin;

            Functions.OPTbracketMinimum(ref a, ref b, ref c, out fa, out fb, out fc, Math.Sin);
            double fmin = Functions.OPTbrentsMinimumSearch(a, b, c, Math.Sin, 1e-6, out xmin);

            double expectedX = 1.5 * Math.PI;
            double expectedF = -1.0;
            PrintComparison("OPT 1D min", xmin, expectedX, fmin, expectedF);
            AssertClose("xmin", xmin, expectedX, 1e-3);
            AssertClose("fmin", fmin, expectedF, 1e-3);
        }

        private static void TestLinprog()
        {
            double[] f = { 1.0, 1.0 };
            double[][] Aieq =
            {
                new[] { -1.0, -1.0 },
                new[] { -1.0,  0.0 },
                new[] {  0.0, -1.0 }
            };
            double[] bieq = { -1.0, 0.0, 0.0 };

            Functions.linprog(out double[] x, out int[] actSet, f, Aieq, bieq, null, null, 1e-10);
            double objective = f[0] * x[0] + f[1] * x[1];

            System.Console.WriteLine($"x = [{x[0]:F6}, {x[1]:F6}], objective = {objective:F6}");
            if (objective < 0.999 || objective > 1.001)
            {
                throw new Exception("linprog objective not within expected tolerance");
            }
            if (x[0] < -1e-8 || x[1] < -1e-8 || x[0] + x[1] < 1.0 - 1e-6)
            {
                throw new Exception("linprog constraints violated");
            }
            System.Console.WriteLine($"active set size = {actSet.Length}");
        }

        private static void TestQuadprog()
        {
            string path = ResolveQpDataPath();
            var data = QpData.Load(path);

            Functions.quadprog(
                out double[] x,
                out _,
                out _,
                out _,
                data.H,
                data.f,
                data.Aieq,
                data.bieq,
                data.Aeq,
                data.beq,
                data.lbd,
                data.ubd,
                data.x0,
                0,
                1e-12,
                1e99,
                200);

            double diff = NormDiff(x, data.xExpected);
            System.Console.WriteLine($"x = [{x[0]:F6}, {x[1]:F6}], expected = [{data.xExpected[0]:F6}, {data.xExpected[1]:F6}]");
            System.Console.WriteLine($"difference norm = {diff:E6}");
            if (diff > 1e-6)
            {
                throw new Exception("quadprog result differs from expected solution");
            }
        }

        private static void TestSqp()
        {
            Functions.SqpObjective obj = (out double f, out double[] grad, double[] x, bool calcGrad) =>
            {
                f = Math.Pow(x[0] - 1.0, 2.0) + Math.Pow(x[1] + 2.0, 2.0);
                if (calcGrad)
                {
                    grad = new[] { 2.0 * (x[0] - 1.0), 2.0 * (x[1] + 2.0) };
                }
                else
                {
                    grad = Array.Empty<double>();
                }
            };

            Functions.SqpConstraints con = (out int nequ, out double[] g, out double[][] grad, double[] x, bool calcGrad) =>
            {
                nequ = 0;
                g = Array.Empty<double>();
                grad = new double[0][];
            };

            double[] lb = { -10.0, -10.0 };
            double[] ub = { 10.0, 10.0 };
            double[] x0 = { 0.0, 0.0 };

            Functions.runSqp(
                out double[] xout,
                out SqpInfo info,
                out _,
                out _,
                obj,
                lb,
                ub,
                null,
                null,
                null,
                null,
                x0,
                con,
                new SqpOptions());

            PrintComparison("SQP", xout[0], 1.0, xout[1], -2.0);
            AssertClose("x0", xout[0], 1.0, 1e-4);
            AssertClose("x1", xout[1], -2.0, 1e-4);
            System.Console.WriteLine($"objective = {info.ObjValue:F6}");
        }

        private static void TestBranchAndBound()
        {
            Functions.SqpObjective obj = (out double f, out double[] grad, double[] x, bool calcGrad) =>
            {
                f = Math.Pow(x[0] - 1.0, 2.0) + Math.Pow(x[1] + 2.0, 2.0);
                if (calcGrad)
                {
                    grad = new[] { 2.0 * (x[0] - 1.0), 2.0 * (x[1] + 2.0) };
                }
                else
                {
                    grad = Array.Empty<double>();
                }
            };

            Functions.SqpConstraints con = (out int nequ, out double[] g, out double[][] grad, double[] x, bool calcGrad) =>
            {
                nequ = 0;
                g = Array.Empty<double>();
                grad = new double[0][];
            };

            int[] xVarType = { 1, 0 };
            double[] lb = { 0.0, -10.0 };
            double[] ub = { 3.0, 10.0 };
            double[] xOrig = { 0.0, 0.0 };

            Functions.OPTbranchAndBound(
                out double[] xout,
                out OptBnbInfo info,
                obj,
                xVarType,
                lb,
                ub,
                null,
                null,
                null,
                null,
                xOrig,
                con,
                new OptBnbOptions());

            PrintComparison("BnB", xout[0], 1.0, xout[1], -2.0);
            AssertClose("x0", xout[0], 1.0, 1e-6);
            AssertClose("x1", xout[1], -2.0, 1e-4);
            System.Console.WriteLine($"objective = {info.ObjValue:F6}");
        }

        private static void PrintComparison(string label, double x0, double expected0, double x1, double expected1)
        {
            System.Console.WriteLine($"{label} result: [{x0:F6}, {x1:F6}]");
            System.Console.WriteLine($"{label} expected: [{expected0:F6}, {expected1:F6}]");
        }

        private static void AssertClose(string name, double actual, double expected, double tol)
        {
            if (Math.Abs(actual - expected) > tol)
            {
                throw new Exception($"{name} differs from expected value");
            }
        }

        private static double NormDiff(double[] a, double[] b)
        {
            if (a.Length != b.Length)
            {
                throw new Exception("vector size mismatch");
            }

            double sum = 0.0;
            for (int i = 0; i < a.Length; i++)
            {
                double d = a[i] - b[i];
                sum += d * d;
            }

            return Math.Sqrt(sum);
        }

        private static string ResolveQpDataPath()
        {
            string baseDir = AppDomain.CurrentDomain.BaseDirectory;
            string path = Path.Combine(baseDir, "qpData.txt");
            if (File.Exists(path))
            {
                return path;
            }

            string repoPath = Path.Combine(baseDir, "..", "..", "..", "MathTools.Console", "qpData.txt");
            repoPath = Path.GetFullPath(repoPath);
            if (File.Exists(repoPath))
            {
                return repoPath;
            }

            throw new Exception("qpData.txt not found");
        }

        private sealed class QpData
        {
            public double[][] H { get; private set; }
            public double[] f { get; private set; }
            public double[][] Aeq { get; private set; }
            public double[] beq { get; private set; }
            public double[][] Aieq { get; private set; }
            public double[] bieq { get; private set; }
            public double[] lbd { get; private set; }
            public double[] ubd { get; private set; }
            public double[] x0 { get; private set; }
            public double[] xExpected { get; private set; }

            public static QpData Load(string path)
            {
                var data = new QpData();
                foreach (string raw in File.ReadAllLines(path))
                {
                    string line = raw.Trim();
                    if (line.Length == 0 || line.StartsWith("#", StringComparison.Ordinal))
                    {
                        continue;
                    }

                    string[] parts = line.Split(new[] { '=' }, 2);
                    if (parts.Length != 2)
                    {
                        continue;
                    }

                    string key = parts[0].Trim();
                    string value = parts[1].Trim();
                    switch (key)
                    {
                        case "H":
                            data.H = ParseMatrix(value);
                            break;
                        case "f":
                            data.f = ParseVector(value);
                            break;
                        case "Aeq":
                            data.Aeq = ParseMatrix(value);
                            break;
                        case "beq":
                            data.beq = ParseVector(value);
                            break;
                        case "Aieq":
                            data.Aieq = ParseMatrix(value);
                            break;
                        case "bieq":
                            data.bieq = ParseVector(value);
                            break;
                        case "lbd":
                            data.lbd = ParseVector(value);
                            break;
                        case "ubd":
                            data.ubd = ParseVector(value);
                            break;
                        case "x0":
                            data.x0 = ParseVector(value);
                            break;
                        case "xExpected":
                            data.xExpected = ParseVector(value);
                            break;
                    }
                }

                data.Aeq ??= new double[0][];
                data.Aieq ??= new double[0][];
                data.beq ??= Array.Empty<double>();
                data.bieq ??= Array.Empty<double>();

                if (data.H == null || data.f == null || data.lbd == null || data.ubd == null || data.x0 == null || data.xExpected == null)
                {
                    throw new Exception("qpData.txt missing required data");
                }

                return data;
            }

            private static double[] ParseVector(string value)
            {
                if (string.IsNullOrWhiteSpace(value))
                {
                    return Array.Empty<double>();
                }

                string[] parts = value.Split(',');
                var result = new double[parts.Length];
                for (int i = 0; i < parts.Length; i++)
                {
                    result[i] = double.Parse(parts[i].Trim(), System.Globalization.CultureInfo.InvariantCulture);
                }

                return result;
            }

            private static double[][] ParseMatrix(string value)
            {
                if (string.IsNullOrWhiteSpace(value))
                {
                    return new double[0][];
                }

                string[] rows = value.Split(';');
                var result = new double[rows.Length][];
                for (int i = 0; i < rows.Length; i++)
                {
                    result[i] = ParseVector(rows[i]);
                }

                return result;
            }
        }
    }
}
