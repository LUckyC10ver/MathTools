using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

namespace MathTools.Core
{
    public class OptBnbOptions
    {
        public TextWriter Log { get; set; }
        public SqpOptions SqpOptions { get; set; }
        public int BranchNonConv { get; set; }

        public OptBnbOptions(
            int display = 0,
            TextWriter log = null,
            int branchNonConv = 0,
            double tolArg = 1e-6,
            double tolObj = 1e-4,
            double tolCon = 1e-6,
            double eps = 1e-11,
            bool checkGrad = false,
            int maxFunEvals = 0,
            double diffMinChange = 1e-9,
            double diffMaxChange = 1e-6,
            double maxDouble = 1.7e38,
            int lineSearch = 0,
            int maxItersQp = 1000)
        {
            Log = log;
            BranchNonConv = branchNonConv;
            SqpOptions = new SqpOptions(
                display,
                log,
                tolArg,
                tolObj,
                tolCon,
                eps,
                checkGrad,
                maxFunEvals,
                diffMinChange,
                diffMaxChange,
                maxDouble,
                null,
                lineSearch,
                maxItersQp);
            SqpOptions.SqpLog = null;
        }
    }

    public class OptBnbInfo
    {
        public double ObjValue { get; set; }
        public int Fail { get; set; }
        public double RunTime { get; set; }
        public int NumCycles { get; set; }
        public SqpInfo Sqp { get; set; }

        public OptBnbInfo()
        {
            ObjValue = 0.0;
            Fail = 0;
            RunTime = 0.0;
            NumCycles = 0;
            Sqp = new SqpInfo();
        }
    }

    public partial class Functions
    {
        public static void OPTbranchAndBound(
            out double[] xout,
            out OptBnbInfo info,
            SqpObjective obj,
            int[] xVarType,
            double[] lb,
            double[] ub,
            double[][] Aequ,
            double[] bequ,
            double[][] Aieq,
            double[] bieq,
            double[] xOrig,
            SqpConstraints con,
            OptBnbOptions opt = null)
        {
            if (obj == null)
            {
                throw new Exception("OPTbranchAndBound: objective is required");
            }

            if (con == null)
            {
                throw new Exception("OPTbranchAndBound: constraints are required");
            }

            if (xOrig == null)
            {
                throw new Exception("OPTbranchAndBound: starting point is required");
            }

            if (xVarType == null || xVarType.Length != xOrig.Length)
            {
                throw new Exception("xVarType must have the same size than xOrig");
            }

            opt ??= new OptBnbOptions();
            info = new OptBnbInfo();
            xout = (double[])xOrig.Clone();

            double inf = opt.SqpOptions.MaxDouble;
            double errx = opt.SqpOptions.TolArg;
            int show = opt.Log != null ? opt.SqpOptions.Display : 0;

            int lx = xOrig.Length;
            double[] xlb = new double[lx];
            double[] xub = new double[lx];

            if (lb != null && lb.Length > lx)
            {
                throw new Exception($"size of lb = {lb.Length} > number of variables");
            }

            if (ub != null && ub.Length > lx)
            {
                throw new Exception($"size of ub = {ub.Length} > number of variables");
            }

            for (int i = 0; i < lx; i++)
            {
                double lower = lb != null && i < lb.Length ? lb[i] : -inf;
                double upper = ub != null && i < ub.Length ? ub[i] : inf;

                switch (xVarType[i])
                {
                    case 0:
                    case 1:
                        xlb[i] = lower;
                        xub[i] = upper;
                        break;
                    case 2:
                        xlb[i] = xOrig[i];
                        xub[i] = xOrig[i];
                        break;
                    default:
                        throw new Exception($"xVarType[{i}] must be equal to one of the integers 0, 1, 2.");
                }

                if (xOrig[i] < xlb[i])
                {
                    throw new Exception($"{i}-th component: starting value = {xOrig[i]} < {xlb[i]} = lower bound");
                }

                if (xOrig[i] > xub[i])
                {
                    throw new Exception($"{i}-th component: starting value = {xOrig[i]} > {xub[i]} = upper bound");
                }

                if (xVarType[i] == 1)
                {
                    if (Math.Abs(xlb[i]) >= inf || !IsInteger(xlb[i]))
                    {
                        throw new Exception($"{i}-th component: lower bound must be an integer for an integer variable");
                    }

                    if (Math.Abs(xub[i]) >= inf || !IsInteger(xub[i]))
                    {
                        throw new Exception($"{i}-th component: upper bound must be an integer for an integer variable");
                    }
                }
            }

            Aequ ??= Array.Empty<double[]>();
            Aieq ??= Array.Empty<double[]>();
            bequ ??= Array.Empty<double>();
            bieq ??= Array.Empty<double>();

            if (Aieq.Length > 0)
            {
                if (Aieq[0].Length != lx)
                {
                    throw new Exception($"number of columns of matrix Aieq = {Aieq[0].Length} != {lx} = number of variables");
                }

                if (bieq.Length != Aieq.Length)
                {
                    throw new Exception($"number of rows of matrix Aieq = {Aieq.Length} != {bieq.Length} = size of vector bieq");
                }
            }
            else if (bieq.Length > 0)
            {
                throw new Exception("number of rows of matrix Aieq = 0 != size of vector bieq");
            }

            if (Aequ.Length > 0)
            {
                if (Aequ[0].Length != lx)
                {
                    throw new Exception($"number of columns of matrix Aequ = {Aequ[0].Length} != {lx} = number of variables");
                }

                if (bequ.Length != Aequ.Length)
                {
                    throw new Exception($"number of rows of matrix Aequ = {Aequ.Length} != {bequ.Length} = size of vector bequ");
                }
            }
            else if (bequ.Length > 0)
            {
                throw new Exception("number of rows of matrix Aequ = 0 != size of vector bequ");
            }

            double zIncumbent = inf;
            double[] xIncumbent = new double[lx];
            for (int i = 0; i < lx; i++)
            {
                xIncumbent[i] = inf;
            }

            int stackCapacity = EstimateStackCapacity(xVarType, xlb, xub);
            var stackX0 = new List<double[]>(stackCapacity) { (double[])xOrig.Clone() };
            var stackXlb = new List<double[]>(stackCapacity) { (double[])xlb.Clone() };
            var stackXub = new List<double[]>(stackCapacity) { (double[])xub.Clone() };
            var stackDepth = new List<int>(stackCapacity) { 1 };

            int[] xchoice = new int[lx];
            if (Aequ.Length > 0)
            {
                int choice = 0;
                for (int i = 0; i < Aequ.Length; i++)
                {
                    bool isBinRow = true;
                    for (int k = 0; k < lx; k++)
                    {
                        double val = Aequ[i][k];
                        if (val != 0.0 && val != 1.0)
                        {
                            isBinRow = false;
                            break;
                        }
                    }

                    if (bequ[i] == 1.0 && isBinRow)
                    {
                        var indices = new List<int>();
                        for (int k = 0; k < lx; k++)
                        {
                            if (Aequ[i][k] == 1.0)
                            {
                                indices.Add(k);
                            }
                        }

                        bool isNormAndNotCont = true;
                        bool isNotFixed = true;
                        bool ifFixedThenZero = true;
                        foreach (int idx in indices)
                        {
                            isNormAndNotCont = isNormAndNotCont && xVarType[idx] != 0 && xchoice[idx] == 0 && xlb[idx] == 0 && xub[idx] == 1;
                            isNotFixed = isNotFixed && xVarType[idx] != 2;
                            if (xVarType[idx] == 2)
                            {
                                ifFixedThenZero = ifFixedThenZero && xOrig[idx] == 0.0;
                            }
                        }

                        if (isNormAndNotCont && (isNotFixed || ifFixedThenZero))
                        {
                            choice++;
                            double sumXorig = 0.0;
                            foreach (int idx in indices)
                            {
                                xchoice[idx] = choice;
                                sumXorig += xOrig[idx];
                            }

                            if (sumXorig == 0.0)
                            {
                                throw new Exception("x0 not correct.");
                            }
                        }
                    }
                }
            }

            if (opt.SqpOptions.MaxFunEvals <= 0)
            {
                opt.SqpOptions.MaxFunEvals = 100 * lx;
            }

            if (show > 3)
            {
                opt.SqpOptions.Display = opt.SqpOptions.Display - 3;
                if (opt.SqpOptions.SqpLog == null)
                {
                    opt.SqpOptions.SqpLog = opt.Log;
                }
            }
            else
            {
                opt.SqpOptions.Display = 0;
                opt.SqpOptions.SqpLog = null;
            }

            info.NumCycles = 0;
            info.Fail = 0;

            var stopwatch = Stopwatch.StartNew();

            while (stackX0.Count > 0)
            {
                info.NumCycles++;

                int last = stackX0.Count - 1;
                double[] x0 = stackX0[last];
                xlb = stackXlb[last];
                xub = stackXub[last];
                int depth = stackDepth[last];

                stackX0.RemoveAt(last);
                stackXlb.RemoveAt(last);
                stackXub.RemoveAt(last);
                stackDepth.RemoveAt(last);

                if (show > 1)
                {
                    info.RunTime = stopwatch.Elapsed.TotalSeconds;
                    double seconds = info.RunTime;
                    int hours = (int)(seconds / 3600.0);
                    seconds -= 3600.0 * hours;
                    int minutes = (int)(seconds / 60.0);
                    seconds -= 60.0 * minutes;

                    double sum = 0.0;
                    for (int i = 0; i < stackDepth.Count; i++)
                    {
                        sum += Math.Pow(0.5, stackDepth[i] - 1.0);
                    }
                    double percdone = Math.Truncate(1000.0 * (1.0 - sum)) / 10.0;

                    opt.Log?.WriteLine($"*** searched {percdone,3:0.0} % of tree");
                    opt.Log?.WriteLine($"*** Z          : {zIncumbent,12:0.0000e+00}");
                    opt.Log?.WriteLine($"*** t          : {hours}h,{minutes}m,{seconds:0.0}s");
                    opt.Log?.WriteLine($"*** c          : {info.NumCycles - 1,12} cycles");
                    opt.Log?.WriteLine($"*** fail       : {info.Fail,12} cycles");
                    if (show > 2)
                    {
                        opt.Log?.WriteLine($"*** stackdepth :  {string.Join("  ", stackDepth)}");
                    }
                    opt.Log?.WriteLine();
                }

                double[] x;
                double z;
                int convflag = 1;
                try
                {
                    runSqp(out x, out info.Sqp, out _, out _, obj, xlb, xub, Aequ, bequ, Aieq, bieq, x0, con, opt.SqpOptions, null);
                    z = info.Sqp.ObjValue;
                }
                catch (Exception ex)
                {
                    convflag = 0;
                    z = inf;
                    if (show > 2)
                    {
                        opt.Log?.WriteLine($"sqp did not converge: {ex.Message}");
                    }
                    x = x0;
                }

                var K = new List<int>();
                for (int i = 0; i < lx; i++)
                {
                    if (xVarType[i] == 1 && xlb[i] < xub[i])
                    {
                        K.Add(i);
                    }
                }

                bool separation = true;
                if ((convflag < 0) || (convflag == 0 && opt.BranchNonConv != 0))
                {
                    separation = false;
                    if (show > 2)
                    {
                        opt.Log?.WriteLine("*** branch pruned");
                    }
                    if (convflag == 0)
                    {
                        info.Fail++;
                        if (show > 2)
                        {
                            opt.Log?.WriteLine("*** not convergent");
                        }
                    }
                    else if (show > 2)
                    {
                        opt.Log?.WriteLine("*** not feasible");
                    }
                }
                else if (z >= zIncumbent && convflag > 0)
                {
                    separation = false;
                    if (show > 2)
                    {
                        opt.Log?.WriteLine("*** branch pruned\n*** ghosted");
                    }
                }
                else if (convflag > 0)
                {
                    bool freeIntRealized = true;
                    foreach (int idx in K)
                    {
                        freeIntRealized = freeIntRealized && Math.Abs(Math.Truncate(x[idx]) - x[idx]) < errx;
                    }

                    if (freeIntRealized)
                    {
                        zIncumbent = z;
                        xIncumbent = (double[])x.Clone();
                        separation = false;
                        if (show > 2)
                        {
                            opt.Log?.WriteLine("*** branch pruned\n*** new best solution found");
                        }
                    }
                }

                if (separation && K.Count > 0)
                {
                    double dzsep = -1.0;
                    int ixsep = K[0];
                    foreach (int ki in K)
                    {
                        double dxsepc = Math.Abs(Math.Truncate(x[ki]) - x[ki]);
                        if (dxsepc >= errx || convflag == 0)
                        {
                            double[] xsepc = (double[])x.Clone();
                            xsepc[ki] = Math.Truncate(x[ki]);
                            obj(out double fval, out _, xsepc, false);
                            double dzsepc = Math.Abs(fval - z);
                            if (dzsepc > dzsep)
                            {
                                dzsep = dzsepc;
                                ixsep = ki;
                            }
                        }
                    }

                    if (xchoice[ixsep] == 0)
                    {
                        BranchOnSingle(
                            stackX0,
                            stackXlb,
                            stackXub,
                            stackDepth,
                            x,
                            xlb,
                            xub,
                            ixsep,
                            depth);
                    }
                    else
                    {
                        BranchOnChoice(
                            stackX0,
                            stackXlb,
                            stackXub,
                            stackDepth,
                            x,
                            xlb,
                            xub,
                            depth,
                            xchoice,
                            K,
                            ixsep);
                    }
                }
                else if (separation && K.Count == 0)
                {
                    info.Fail++;
                    if (show > 2)
                    {
                        opt.Log?.WriteLine("*** branch pruned\n*** leaf not convergent");
                    }
                }
            }

            info.RunTime = stopwatch.Elapsed.TotalSeconds;
            info.ObjValue = zIncumbent;
            xout = new double[lx];
            for (int i = 0; i < lx; i++)
            {
                xout[i] = xVarType[i] == 1 ? Math.Truncate(xIncumbent[i]) : xIncumbent[i];
            }

            if (show > 0 && opt.Log != null)
            {
                opt.Log.WriteLine("branch and bound terminated.");
                opt.Log.WriteLine($"incumbent objective value : {info.ObjValue}");
                double seconds = info.RunTime;
                int hours = (int)(seconds / 3600.0);
                seconds -= 3600.0 * hours;
                int minutes = (int)(seconds / 60.0);
                seconds -= 60.0 * minutes;
                opt.Log.WriteLine($"time needed               : {hours}h, {minutes}m, {seconds:0.0}s");
                if (show > 1)
                {
                    opt.Log.WriteLine($"incumbent X               : {string.Join(", ", xout)}");
                }
            }
        }

        private static void BranchOnSingle(
            List<double[]> stackX0,
            List<double[]> stackXlb,
            List<double[]> stackXub,
            List<int> stackDepth,
            double[] x,
            double[] xlb,
            double[] xub,
            int ixsep,
            int depth)
        {
            double low = xlb[ixsep];
            double high = xub[ixsep];
            int sepDepth = depth;
            bool branch = true;
            double[] domain = new[] { low, high };

            while (branch)
            {
                double xboundary = (domain[0] + domain[1]) / 2.0;
                double[] domainA;
                double[] domainB;
                if (x[ixsep] < xboundary)
                {
                    domainA = new[] { domain[0], Math.Truncate(xboundary) };
                    domainB = new[] { Math.Truncate(xboundary + 1.0), domain[1] };
                }
                else
                {
                    domainA = new[] { Math.Truncate(xboundary + 1.0), domain[1] };
                    domainB = new[] { domain[0], Math.Truncate(xboundary) };
                }

                sepDepth++;
                PushNode(stackX0, stackXlb, stackXub, stackDepth, x, xlb, xub, ixsep, domainB, sepDepth);

                if (domainA[0] == domainA[1])
                {
                    PushNode(stackX0, stackXlb, stackXub, stackDepth, x, xlb, xub, ixsep, domainA, sepDepth);
                    branch = false;
                }
                else
                {
                    domain = domainA;
                }
            }
        }

        private static void BranchOnChoice(
            List<double[]> stackX0,
            List<double[]> stackXlb,
            List<double[]> stackXub,
            List<int> stackDepth,
            double[] x,
            double[] xlb,
            double[] xub,
            int depth,
            int[] xchoice,
            List<int> K,
            int ixsep)
        {
            var L = new List<int>();
            for (int i = 0; i < xchoice.Length; i++)
            {
                if (xchoice[i] == xchoice[ixsep])
                {
                    L.Add(i);
                }
            }

            var M = Intersect(K, L);
            var xM = new List<double>();
            foreach (int idx in M)
            {
                xM.Add(x[idx]);
            }

            var N = SortIndicesByValue(xM);
            int halfN = N.Count / 2;
            var part1 = new List<int>();
            var part2 = new List<int>();
            for (int i = 0; i < halfN; i++)
            {
                part1.Add(M[N[i]]);
            }
            for (int i = halfN; i < N.Count; i++)
            {
                part2.Add(M[N[i]]);
            }

            int sepDepth = depth + 1;
            PushChoiceNode(stackX0, stackXlb, stackXub, stackDepth, x, xlb, xub, part1, part2, sepDepth);
            PushChoiceNode(stackX0, stackXlb, stackXub, stackDepth, x, xlb, xub, part2, part1, sepDepth);
        }

        private static void PushChoiceNode(
            List<double[]> stackX0,
            List<double[]> stackXlb,
            List<double[]> stackXub,
            List<int> stackDepth,
            double[] x,
            double[] xlb,
            double[] xub,
            List<int> partKeep,
            List<int> partZero,
            int depth)
        {
            double[] x0 = (double[])x.Clone();
            double sum = 0.0;
            foreach (int idx in partKeep)
            {
                sum += x0[idx];
            }
            double offset = (1.0 - sum) / partKeep.Count;
            foreach (int idx in partKeep)
            {
                x0[idx] += offset;
            }

            double[] newL = (double[])xlb.Clone();
            double[] newU = (double[])xub.Clone();
            foreach (int idx in partZero)
            {
                newU[idx] = 0.0;
            }

            stackX0.Add(x0);
            stackXlb.Add(newL);
            stackXub.Add(newU);
            stackDepth.Add(depth);
        }

        private static void PushNode(
            List<double[]> stackX0,
            List<double[]> stackXlb,
            List<double[]> stackXub,
            List<int> stackDepth,
            double[] x,
            double[] xlb,
            double[] xub,
            int ixsep,
            double[] domain,
            int depth)
        {
            double[] x0 = (double[])x.Clone();
            double[] newL = (double[])xlb.Clone();
            double[] newU = (double[])xub.Clone();
            newL[ixsep] = domain[0];
            newU[ixsep] = domain[1];
            stackX0.Add(x0);
            stackXlb.Add(newL);
            stackXub.Add(newU);
            stackDepth.Add(depth);
        }

        private static int EstimateStackCapacity(int[] xVarType, double[] xlb, double[] xub)
        {
            double sum = 0.0;
            int count = 0;
            for (int i = 0; i < xVarType.Length; i++)
            {
                if (xVarType[i] == 1)
                {
                    double d = xub[i] - xlb[i] + 1.0;
                    double ld = Math.Log(d) / Math.Log(2.0);
                    sum += ld;
                    count++;
                }
            }
            return (int)Math.Ceiling(sum + count + 1.0);
        }

        private static bool IsInteger(double value)
        {
            return Math.Abs(value - Math.Truncate(value)) < 1e-12;
        }

        private static List<int> Intersect(List<int> a, List<int> b)
        {
            var set = new HashSet<int>(b);
            var result = new List<int>();
            foreach (int v in a)
            {
                if (set.Contains(v))
                {
                    result.Add(v);
                }
            }
            return result;
        }

        private static List<int> SortIndicesByValue(List<double> values)
        {
            var indices = new List<int>();
            for (int i = 0; i < values.Count; i++)
            {
                indices.Add(i);
            }
            indices.Sort((a, b) => values[a].CompareTo(values[b]));
            return indices;
        }
    }
}
